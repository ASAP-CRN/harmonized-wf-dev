version 1.0

struct Sample {
	String sample_id
	String batch

	File fastq_R1
	File fastq_R2
}

workflow harmonized_pmdbs_analysis {
	input {
		String project_name
		Array[Sample] samples

		Float soup_rate = 0.20

		Int clustering_algorithm = 3
		Float clustering_resolution = 0.3
		File cell_type_markers_list

		Array[String] groups = ["sample", "batch", "seurat_clusters"]
		Array[String] features = ["doublet_scores", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"]

		String container_registry
	}

	scatter (sample in samples) {
		call cellranger {
			input:
				sample_id = sample.sample_id,
				fastq_R1 = sample.fastq_R1,
				fastq_R2 = sample.fastq_R2,
				container_registry = container_registry
		}

		call preprocess {
			input:
				sample_id = sample.sample_id,
				batch = sample.batch,
				raw_counts = cellranger.raw_counts,
				filtered_counts = cellranger.filtered_counts,
				soup_rate = soup_rate,
				container_registry = container_registry
		}
	}

	call doublets {
		input:
			project_name = project_name,
			preprocessed_seurat_objects = preprocess.preprocessed_seurat_object,
			container_registry = container_registry
	}

	call plot_qc {
		input:
			project_name = project_name,
			unfiltered_metadata = doublets.unfiltered_metadata,
			container_registry = container_registry
	}

	scatter (preprocessed_seurat_object in preprocess.preprocessed_seurat_object) {
		call filter {
			input:
				preprocessed_seurat_object = preprocessed_seurat_object,
				unfiltered_metadata = doublets.unfiltered_metadata,
				container_registry = container_registry
		}

		call process {
			input:
				filtered_seurat_object = filter.filtered_seurat_object,
				container_registry = container_registry
		}
	}

	call harmony {
		input:
			project_name = project_name,
			normalized_seurat_objects = process.normalized_seurat_object,
			container_registry = container_registry
	}

	call neighbors {
		input:
			harmony_seurat_object = harmony.harmony_seurat_object,
			container_registry = container_registry
	}

	call umap {
		input:
			neighbors_seurat_object = neighbors.neighbors_seurat_object,
			container_registry = container_registry
	}

	call cluster {
		input:
			project_name = project_name,
			umap_seurat_object = umap.umap_seurat_object,
			clustering_algorithm = clustering_algorithm,
			clustering_resolution = clustering_resolution,
			cell_type_markers_list = cell_type_markers_list,
			container_registry = container_registry
	}

	call sctype {
		input:
			project_name = project_name,
			cluster_seurat_object = cluster.cluster_seurat_object,
			cell_type_markers_list = cell_type_markers_list,
			container_registry = container_registry
	}

	# TODO probably more overhead than needed to scatter this rather than parallelizing within-task
	scatter (group in groups) {
		call plot_groups {
			input:
				project_name = project_name,
				metadata = sctype.metadata,
				group = group,
				container_registry = container_registry
		}
	}

	# TODO same overhead as plot_groups
	scatter (feature in features) {
		call plot_features {
			input:
				project_name = project_name,
				metadata = sctype.metadata,
				feature = feature,
				container_registry = container_registry
		}
	}

	output {
		# QC plots
		File qc_plot1 = plot_qc.plot1
		File qc_plot2 = plot_qc.plot2

		# Final metadata
		File metadata = sctype.metadata

		# Group and feature plots
		Array[File] group_umap_plots = plot_groups.group_umap_plot
		Array[File] feature_umap_plots = plot_features.feature_umap_plot

		# Clustered seurat object
		File cluster_seurat_object = cluster.cluster_seurat_object
	}

	meta {
		description: "Harmonized postmortem-derived brain sequencing (PMDBS) workflow"
	}

	parameter_meta {
		project_name: {help: "Name of project"}
		samples: {help: "The set of samples and their associated reads and metadata"}
		soup_rate: {help: "Dataset contamination rate fraction; used to remove mRNA contamination from the RNAseq data [0.2]"}
		clustering_algorithm: {help: "Clustering algorithm to use. [3]"}
		clustering_resolution: {help: "Clustering resolution to use during clustering. [0.3]"}
		cell_type_markers_list: {help: "Seurat object RDS file containing a list of major cell type markers; used to annotate clusters."}
		groups: {help: "Groups to produce umap plots for. ['sample', 'batch', 'seurat_clusters']"}
		features: {help: "Features to produce umap plots for. ['doublet_scores', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb']"}
		container_registry: {help: "Container registry where Docker images are hosted"}
	}
}

task cellranger {
	input {
		String sample_id

		File fastq_R1
		File fastq_R2

		String container_registry
	}

	command <<<
		set -euo pipefail

		echo "cellranger ~{fastq_R1} ~{fastq_R2}"

		touch ~{sample_id}_count_raw_feature_bc_matrix.h5
		touch ~{sample_id}_count_filtered_feature_bc_matrix.h5
		touch ~{sample_id}_count_molecule_info.h5
	>>>

	output {
		File raw_counts = "~{sample_id}_count_raw_feature_bc_matrix.h5"
		File filtered_counts = "~{sample_id}_count_filtered_feature_bc_matrix.h5"
		File molecule_info = "~{sample_id}_count_molecule_info.h5"
	}

	runtime {
		docker: "~{container_registry}/cellranger" # TODO
	}
}

task preprocess {
	input {
		String sample_id
		String batch

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
			--sample-id ~{sample_id} \
			--batch ~{batch} \
			--raw-counts ~{raw_counts} \
			--filtered-counts ~{filtered_counts} \
			--soup-rate ~{soup_rate} \
			--output-seurat-object seurat_object_~{sample_id}_preprocessed_01.rds
	>>>

	output {
		File preprocessed_seurat_object = "seurat_object_~{sample_id}_preprocessed_01.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
	}
}

task doublets {
	input {
		String project_name
		Array[File] preprocessed_seurat_objects

		String container_registry
	}

	Int threads = 2

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/gmm_doublet_calling.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-objects ~{sep=' ' preprocessed_seurat_objects} \
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
		File preprocessed_seurat_object
		File unfiltered_metadata

		String container_registry
	}

	String seurat_object_basename = basename(preprocessed_seurat_object, "_01.rds")

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/filter.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{preprocessed_seurat_object} \
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
		File filtered_seurat_object

		String container_registry
	}

	Int threads = 2
	String seurat_object_basename = basename(filtered_seurat_object, "_02.rds")

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/filter.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-object ~{filtered_seurat_object} \
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

task harmony {
	input {
		String project_name
		Array[File] normalized_seurat_objects

		String container_registry
	}

	Int threads = 8

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/filter.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-objects ~{sep=' ' normalized_seurat_objects} \
			--output-seurat-object ~{project_name}_seurat_object_harmony_integrated_04.rds
	>>>

	output {
		File harmony_seurat_object = "~{project_name}_seurat_object_harmony_integrated_04.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: threads
	}
}

task neighbors {
	input {
		File harmony_seurat_object

		String container_registry
	}

	String harmony_seurat_object_basename = basename(harmony_seurat_object, "_04.rds")

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/find_neighbors.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{harmony_seurat_object} \
			--output-seurat-object ~{harmony_seurat_object_basename}_neighbors_05.rds
	>>>

	output {
		File neighbors_seurat_object = "~{harmony_seurat_object_basename}_neighbors_05.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
	}
}

task umap {
	input {
		File neighbors_seurat_object

		String container_registry
	}

	String neighbors_seurat_object_basename = basename(neighbors_seurat_object, "_05.rds")

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/find_neighbors.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{neighbors_seurat_object} \
			--output-seurat-object ~{neighbors_seurat_object_basename}_umap_06.rds
	>>>

	output {
		File umap_seurat_object = "~{neighbors_seurat_object_basename}_umap_06.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
	}
}

task cluster {
	input {
		String project_name
		File umap_seurat_object

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		String container_registry
	}

	Int threads = 8
	String umap_seurat_object_basename = basename(umap_seurat_object, "_06.rds")

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/clustering.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-object ~{umap_seurat_object} \
			--clustering-algorithm ~{clustering_algorithm} \
			--clustering-resolution ~{clustering_resolution} \
			--cell-type-markers-list ~{cell_type_markers_list} \
			--output-cell-type-plot ~{project_name}.major_type_module_umap.pdf \
			--output-seurat-object ~{umap_seurat_object_basename}_cluster_07.rds
	>>>

	output {
		File major_cell_type_plot = "~{project_name}.major_type_module_umap.pdf"
		File cluster_seurat_object = "~{umap_seurat_object_basename}_cluster_07.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: threads
	}
}

task sctype {
	input {
		String project_name
		File cluster_seurat_object

		File cell_type_markers_list

		String container_registry
	}

	Int threads = 8

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/annotate_clusters.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-object ~{cluster_seurat_object} \
			--cell-type-markers-list ~{cell_type_markers_list} \
			--output-metadata-file ~{project_name}.final_metadata.csv
	>>>

	output {
		File metadata = "~{project_name}.final_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: threads
	}
}

task plot_groups {
	input {
		String project_name
		File metadata

		String group

		String container_registry
	}

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/plot_groups.R \
			--working-dir "$(pwd)" \
			--metadata ~{metadata} \
			--group ~{group} \
			--output-group-umap-plot ~{project_name}.~{group}_group_umap.pdf
	>>>

	output {
		File group_umap_plot = "~{project_name}.~{group}_group_umap.pdf"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
	}
}

task plot_features {
	input {
		String project_name
		File metadata

		String feature

		String container_registry
	}

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/plot_features.R \
			--working-dir "$(pwd)" \
			--metadata ~{metadata} \
			--feature ~{feature} \
			--output-feature-umap-plot ~{project_name}.~{feature}_feature_umap.pdf
	>>>

	output {
		File feature_umap_plot = "~{project_name}.~{feature}_feature_umap.pdf"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
	}
}
