version 1.0

# Generate a preprocessed AnnData object

import "../structs.wdl"
import "../common/check_output_files_exist.wdl" as CheckOutputFilesExist

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

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{workflow_name}"
	String cellbender_raw_data_path = "~{workflow_raw_data_path_prefix}/remove_technical_artifacts/~{workflow_version}"
	String adata_raw_data_path = "~{workflow_raw_data_path_prefix}/counts_to_adata/~{workflow_version}"

	scatter (sample_object in samples) {
		String cellbender_count_output = "~{cellbender_raw_data_path}/~{sample_object.sample_id}.cellbender.h5"
	}

	call CheckOutputFilesExist.check_output_files_exist as check_cellbender_outputs_exist {
		input:
			output_files = cellbender_count_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (index in range(length(samples))) {
		Sample sample = samples[index]

		Array[String] project_sample_id = [project_id, sample.sample_id]

		String cellbender_remove_background_complete = check_cellbender_outputs_exist.sample_complete[index][0]

		String cellbender_report_html = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender_report.html"
		String cellbender_removed_background_counts = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender.h5"
		String cellbender_filtered_removed_background_counts = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender_filtered.h5"
		String cellbender_cell_barcodes_csv = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender_cell_barcodes.csv"
		String cellbender_graph_pdf = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender.pdf"
		String cellbender_log = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender.log"
		String cellbender_metrics_csv = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender_metrics.csv"
		String cellbender_checkpoint_tar_gz = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbender_ckpt.tar.gz"
		String cellbender_posterior_probability = "~{cellbender_raw_data_path}/~{sample.sample_id}.cellbend_posterior.h5"

		if (cellbender_remove_background_complete == "false") {
			call remove_technical_artifacts {
				input:
					sample_id = sample.sample_id,
					raw_counts = raw_counts[index],
					cellbender_fpr = cellbender_fpr,
					raw_data_path = cellbender_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File report_html_output = select_first([remove_technical_artifacts.report_html, cellbender_report_html]) #!FileCoercion
		File removed_background_counts_output = select_first([remove_technical_artifacts.removed_background_counts, cellbender_removed_background_counts]) #!FileCoercion
		File filtered_removed_background_counts_output = select_first([remove_technical_artifacts.filtered_removed_background_counts, cellbender_filtered_removed_background_counts]) #!FileCoercion
		File cell_barcodes_csv_output = select_first([remove_technical_artifacts.cell_barcodes_csv, cellbender_cell_barcodes_csv]) #!FileCoercion
		File graph_pdf_output = select_first([remove_technical_artifacts.graph_pdf, cellbender_graph_pdf]) #!FileCoercion
		File log_output = select_first([remove_technical_artifacts.log, cellbender_log]) #!FileCoercion
		File metrics_csv_output = select_first([remove_technical_artifacts.metrics_csv, cellbender_metrics_csv]) #!FileCoercion
		File checkpoint_tar_gz_output = select_first([remove_technical_artifacts.checkpoint_tar_gz, cellbender_checkpoint_tar_gz]) #!FileCoercion
		File posterior_probability_output = select_first([remove_technical_artifacts.posterior_probability, cellbender_posterior_probability]) #!FileCoercion

		call counts_to_adata {
			input:
				sample_id = sample.sample_id,
				batch = select_first([sample.batch]),
				project_id = project_id,
				cellbender_counts = removed_background_counts_output,
				raw_data_path = adata_raw_data_path,
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
		Array[File] report_html = report_html_output
		Array[File] removed_background_counts = removed_background_counts_output #!FileCoercion
		Array[File] filtered_removed_background_counts = filtered_removed_background_counts_output #!FileCoercion
		Array[File] cell_barcodes_csv = cell_barcodes_csv_output #!FileCoercion
		Array[File] graph_pdf = graph_pdf_output #!FileCoercion
		Array[File] log = log_output #!FileCoercion
		Array[File] metrics_csv = metrics_csv_output #!FileCoercion
		Array[File] checkpoint_tar_gz = checkpoint_tar_gz_output #!FileCoercion
		Array[File] posterior_probability = posterior_probability_output #!FileCoercion

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

		mv ckpt.tar.gz "~{sample_id}.cellbender_ckpt.tar.gz"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.cellbender_report.html" \
			-o "~{sample_id}.cellbender.h5" \
			-o "~{sample_id}.cellbender_filtered.h5" \
			-o "~{sample_id}.cellbender_cell_barcodes.csv" \
			-o "~{sample_id}.cellbender.pdf" \
			-o "~{sample_id}.cellbender.log" \
			-o "~{sample_id}.cellbender_metrics.csv" \
			-o "~{sample_id}.cellbender_ckpt.tar.gz" \
			-o "~{sample_id}.cellbend_posterior.h5"
	>>>

	output {
		String report_html = "~{raw_data_path}/~{sample_id}.cellbender_report.html"
		String removed_background_counts = "~{raw_data_path}/~{sample_id}.cellbender.h5"
		String filtered_removed_background_counts = "~{raw_data_path}/~{sample_id}.cellbender_filtered.h5"
		String cell_barcodes_csv = "~{raw_data_path}/~{sample_id}.cellbender_cell_barcodes.csv"
		String graph_pdf = "~{raw_data_path}/~{sample_id}.cellbender.pdf"
		String log = "~{raw_data_path}/~{sample_id}.cellbender.log"
		String metrics_csv = "~{raw_data_path}/~{sample_id}.cellbender_metrics.csv"
		String checkpoint_tar_gz = "~{raw_data_path}/~{sample_id}.cellbender_ckpt.tar.gz"
		String posterior_probability = "~{raw_data_path}/~{sample_id}.cellbend_posterior.h5"
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
		gpuCount: 8
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
		docker: "~{container_registry}/scvi:1.1.0"
		cpu: 4
		memory: "8 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}
