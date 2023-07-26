version 1.0

workflow harmonized_pmdbs_analysis {
	input {
		String dataset
		File batch_file

		File raw_counts
		File filtered_counts

		Float soup_rate

		String container_registry
	}

	call preprocess {
		input:
			dataset = dataset,
			batch_file = batch_file,
			raw_counts = raw_counts,
			filtered_counts = filtered_counts,
			soup_rate = soup_rate,
			container_registry = container_registry
	}

	output {
		File seurat_object = preprocess.seurat_object
	}

	meta {
		description: "Harmonized postmortem-derived brain sequencing (PMDBS) workflow"
	}

	parameter_meta {
		dataset: {help: ""}
		batch_file: {help: ""}
		raw_counts: {help: ""}
		filtered_counts: {help: ""}
		soup_rate: {help: ""}
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
		docker: "~{container_registry}/multiome:eaaeb73"
	}
}
