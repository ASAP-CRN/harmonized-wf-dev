version 1.0

workflow harmonized_pmdbs_analysis {
	input {
		String project_name
		Array[String] datasets
		File batch_file

		File raw_counts
		File filtered_counts

		Float soup_rate = 0.20

		String container_registry
	}

	scatter (dataset in datasets) {
		call preprocess {
			input:
				dataset = dataset,
				batch_file = batch_file,
				raw_counts = raw_counts,
				filtered_counts = filtered_counts,
				soup_rate = soup_rate,
				container_registry = container_registry
		}
	}

	call doublets {
		input:
			project_name = project_name,
			seurat_objects = preprocess.seurat_object,
			container_registry = container_registry
	}

	call plot_qc {
		input:
			project_name = project_name,
			unfiltered_metadata = doublets.unfiltered_metadata,
			container_registry = container_registry
	}

	scatter (seurat_object in preprocess.seurat_object) {
		call filter {
			input:
				seurat_object = seurat_object,
				unfiltered_metadata = doublets.unfiltered_metadata,
				container_registry = container_registry
		}

		call process {
			input:
				seurat_object = filter.filtered_seurat_object,
				container_registry = container_registry
		}
	}

	output {
		Array[File] seurat_objects = preprocess.seurat_object
		File unfiltered_metadata = doublets.unfiltered_metadata
		File qc_plot1 = plot_qc.plot1
		File qc_plot2 = plot_qc.plot2
		Array[File] filtered_seurat_objects = filter.filtered_seurat_object
		Array[File] normalized_seurat_objects = process.normalized_seurat_object
	}

	meta {
		description: "Harmonized postmortem-derived brain sequencing (PMDBS) workflow"
	}

	parameter_meta {
		project_name: {help: "Name of project"}
		datasets: {help: "Array of dataset names to process"}
		batch_file: {help: "Samples CSV file"}
		raw_counts: {help: "Unfiltered feature-barcode matrices HDF5 output by cellranger"}
		filtered_counts: {help: "Filtered feature-barcode matrices HDF5 output by cellranger"}
		soup_rate: {help: "Dataset contamination rate fraction; used to remove mRNA contamination from the RNAseq data [0.2]"}
		container_registry: {help: "Container registry where Docker images are hosted"}
	}
}

task preprocess {
	input {
		String dataset

		File batch_file
		File raw_counts
		File filtered_counts

		Float soup_rate

		String container_registry
	}

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/preprocess.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--dataset ~{dataset} \
			--batch-file ~{batch_file} \
			--raw-counts ~{raw_counts} \
			--filtered-counts ~{filtered_counts} \
			--soup-rate ~{soup_rate} \
			--output-seurat-object seurat_object_~{dataset}_preprocessed_01.rds
	>>>

	output {
		File seurat_object = "seurat_object_~{dataset}_preprocessed_01.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
	}
}

task doublets {
	input {
		String project_name
		Array[File] seurat_objects

		String container_registry
	}

	Int threads = 2

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/gmm_doublet_calling.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-objects ~{sep=' ' seurat_objects} \
			--project-name ~{project_name} \
			--output-metadata-file ~{project_name}.unfiltered_metadata.csv
	>>>

	output {
		File unfiltered_metadata = "~{project_name}.unfiltered_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: threads
	}
}

task plot_qc {
	input {
		String project_name
		File unfiltered_metadata

		String container_registry
	}

	Int threads = 2

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/plot_qc_metrics.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--metadata ~{unfiltered_metadata} \
			--project-name ~{project_name} \
			--plot1-output-file ~{project_name}.qc_plot1.pdf \
			--plot2-output-file ~{project_name}.qc_plot2.pdf
	>>>

	# TODO come up with better output names for these plots, both here and in the Rscript
	output {
		File plot1 = "~{project_name}.qc_plot1.pdf"
		File plot2 = "~{project_name}.qc_plot2.pdf"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: threads
	}
}

task filter {
	input {
		File seurat_object
		File unfiltered_metadata

		String container_registry
	}

	String seurat_object_basename = basename(seurat_object, "_01.rds")

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/filter.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{seurat_object} \
			--metadata ~{unfiltered_metadata} \
			--output-seurat-object ~{seurat_object_basename}_filtered_02.rds
	>>>

	output {
		File filtered_seurat_object = "~{seurat_object_basename}_filtered_02.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
	}
}

task process {
	input {
		File seurat_object

		String container_registry
	}

	Int threads = 2
	String seurat_object_basename = basename(seurat_object, "_02.rds")

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/filter.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-object ~{seurat_object} \
			--output-seurat-object ~{seurat_object_basename}_filtered_normalized_03.rds
	>>>

	output {
		File normalized_seurat_object = "~{seurat_object_basename}_filtered_normalized_03.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: threads
	}
}
