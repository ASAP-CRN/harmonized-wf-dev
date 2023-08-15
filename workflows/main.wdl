version 1.0

workflow harmonized_pmdbs_analysis {
	input {
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

	output {
		Array[File] seurat_objects = preprocess.seurat_object
	}

	meta {
		description: "Harmonized postmortem-derived brain sequencing (PMDBS) workflow"
	}

	parameter_meta {
		dataset: {help: "Array of dataset names to process"}
		batch_file: {help: "Samples CSV file"}
		raw_counts: {help: "Unfiltered feature-barcode matrices HDF5 output by cellranger"}
		filtered_counts: {help: "Filtered feature-barcode matrices HDF5 output by cellranger"}
		soup_rate: {help: "Dataset contamination rate fraction; used to remove mRNA contamination the RNAseq data. [0.2]"}
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
			--seurat-object seurat_object_~{dataset}_preprocessed_01.rds
	>>>

	output {
		File seurat_object = "seurat_object_~{dataset}_preprocessed_01.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
	}
}
