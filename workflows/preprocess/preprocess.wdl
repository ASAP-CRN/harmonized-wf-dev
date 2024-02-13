version 1.0

# Generate a preprocessed AnnData object

import "../structs.wdl"

workflow preprocess {
	input {
		String project_id
		Array[Sample] samples

		Array[File] raw_counts

		Float cellbender_fpr

		String run_timestamp
		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	String workflow_name = "preprocess"
	String workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version]]

	String raw_data_path = "~{raw_data_path_prefix}/~{workflow_name}"

	scatter (index in range(length(samples))) {
		Sample sample = samples[index]

		Array[String] project_sample_id = [project_id, sample.sample_id]

		call remove_technical_artifacts {
			input:
				sample_id = sample.sample_id,
				raw_counts = raw_counts[index],
				cellbender_fpr = cellbender_fpr,
				raw_data_path = "~{raw_data_path}/remove_technical_artifacts/~{workflow_version}",
				workflow_info = workflow_info,
				billing_project = billing_project,
				container_registry = container_registry,
				zones = zones
		}

		call counts_to_adata {
			input:
				sample_id = sample.sample_id,
				batch = select_first([sample.batch]),
				project_id = project_id,
				cellbender_counts = remove_technical_artifacts.remove_background_counts, #!FileCoercion
				raw_data_path = "~{raw_data_path}/counts_to_adata/~{workflow_version}",
				workflow_info = workflow_info,
				billing_project = billing_project,
				container_registry = container_registry,
				zones = zones
		}
	}

	output {
		# Sample list
		Array[Array[String]] project_sample_ids = project_sample_id

		# Remove technical artifacts - Cellbender
		Array[File] report_html = remove_technical_artifacts.report_html #!FileCoercion
		Array[File] remove_background_counts = remove_technical_artifacts.remove_background_counts #!FileCoercion
		Array[File] filtered_remove_background_counts = remove_technical_artifacts.filtered_remove_background_counts #!FileCoercion
		Array[File] cell_barcodes_csv = remove_technical_artifacts.cell_barcodes_csv #!FileCoercion
		Array[File] graph_pdf = remove_technical_artifacts.graph_pdf #!FileCoercion
		Array[File] log = remove_technical_artifacts.log #!FileCoercion
		Array[File] metrics_csv = remove_technical_artifacts.metrics_csv #!FileCoercion
		Array[File] checkpoint_tar_gz = remove_technical_artifacts.checkpoint_tar_gz #!FileCoercion
		Array[File] posterior_probability = remove_technical_artifacts.posterior_probability #!FileCoercion

		# AnnData counts
		Array[File] adata_object = counts_to_adata.preprocessed_adata_object #!FileCoercion
	}
}

task remove_technical_artifacts {
	input {
		String sample_id

		File raw_counts

		Float cellbender_fpr

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int disk_size = ceil(size(raw_counts, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/cellbender.py \
			--raw-counts ~{raw_counts} \
			--output-name ~{sample_id}.cellbender. \
			--fpr ~{cellbender_fpr}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.cellbender.output_report.html" \
			-o "~{sample_id}.cellbender.output.h5" \
			-o "~{sample_id}.cellbender.output_filtered.h5" \
			-o "~{sample_id}.cellbender.output_cell_barcodes.csv" \
			-o "~{sample_id}.cellbender.output.pdf" \
			-o "~{sample_id}.cellbender.output.log" \
			-o "~{sample_id}.cellbender.output_metrics.csv" \
			-o "~{sample_id}.cellbender.ckpt.tar.gz" \
			-o "~{sample_id}.cellbender.output_posterior.h5"
	>>>

	output {
		String report_html = "~{raw_data_path}/~{sample_id}.cellbender.output_report.html"
		String remove_background_counts = "~{raw_data_path}/~{sample_id}.cellbender.output.h5"
		String filtered_remove_background_counts = "~{raw_data_path}/~{sample_id}.cellbender.output_filtered.h5"
		String cell_barcodes_csv = "~{raw_data_path}/~{sample_id}.cellbender.output_cell_barcodes.csv"
		String graph_pdf = "~{raw_data_path}/~{sample_id}.cellbender.output.pdf"
		String log = "~{raw_data_path}/~{sample_id}.cellbender.output.log"
		String metrics_csv = "~{raw_data_path}/~{sample_id}.cellbender.output_metrics.csv"
		String checkpoint_tar_gz = "~{raw_data_path}/~{sample_id}.cellbender.ckpt.tar.gz"
		String posterior_probability = "~{raw_data_path}/~{sample_id}.cellbender.output_posterior.h5"
	}

	runtime {
		docker: "~{container_registry}/cellbender:0.3.0"
		cpu: 4
		memory: "16 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
		gpuType: "nvidia-tesla-k80"
		gpuCount: 1
	}
}

task counts_to_adata {
	input {
		String sample_id
		String batch
		String project_id

		File cellbender_counts

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int disk_size = ceil(size(cellbender_counts, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/main/preprocess.py \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--adata-input ~{cellbender_counts} \
			--sample-id ~{sample_id} \
			--batch ~{batch} \
			--project ~{project_id} \
			--adata-output ~{sample_id}.adata_object.h5ad.gz

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.adata_object.h5ad.gz"
	>>>

	output {
		String preprocessed_adata_object = "~{raw_data_path}/~{sample_id}.adata_object.h5ad.gz"
	}

	runtime {
		docker: "~{container_registry}/scvi:1.0.4"
		cpu: 4
		memory: "8 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}
