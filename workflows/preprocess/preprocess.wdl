version 1.0

# Generate a preprocessed AnnData object

import "../structs.wdl"

workflow preprocess {
	input {
		String project_id
		Array[Sample] samples

		File cellranger_reference_data

		Float cellbender_fpr

		String run_timestamp
		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	# Task and subworkflow versions
	String workflow_name = "preprocess"
	String cellranger_task_version = "1.1.0"
	String workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{workflow_name}"
	String cellranger_raw_data_path = "~{workflow_raw_data_path_prefix}/cellranger/~{cellranger_task_version}"
	String cellbender_raw_data_path = "~{workflow_raw_data_path_prefix}/remove_technical_artifacts/~{workflow_version}"
	String adata_raw_data_path = "~{workflow_raw_data_path_prefix}/counts_to_adata/~{workflow_version}"

	scatter (sample_object in samples) {
		String cellranger_count_output = "~{cellranger_raw_data_path}/~{sample_object.sample_id}.raw_feature_bc_matrix.h5"
		String cellbender_count_output = "~{cellbender_raw_data_path}/~{sample_object.sample_id}.cellbender.h5"
		String adata_object_output = "~{adata_raw_data_path}/~{sample_object.sample_id}.adata_object.h5ad"
	}

	# For each sample, outputs an array of true/false: [cellranger_counts_complete, remove_technical_artifacts_complete]
	call check_output_files_exist {
		input:
			cellranger_count_output_files = cellranger_count_output,
			remove_technical_artifacts_output_files = cellbender_count_output,
			adata_object_output_files = adata_object_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (index in range(length(samples))) {
		Sample sample = samples[index]

		Array[String] project_sample_id = [project_id, sample.sample_id]

		String cellranger_count_complete = check_output_files_exist.sample_preprocessing_complete[index][0]
		String cellbender_remove_background_complete = check_output_files_exist.sample_preprocessing_complete[index][1]
		String adata_object_complete = check_output_files_exist.sample_preprocessing_complete[index][2]

		String cellranger_raw_counts = "~{cellranger_raw_data_path}/~{sample.sample_id}.raw_feature_bc_matrix.h5"
		String cellranger_filtered_counts = "~{cellranger_raw_data_path}/~{sample.sample_id}.filtered_feature_bc_matrix.h5"
		String cellranger_molecule_info = "~{cellranger_raw_data_path}/~{sample.sample_id}.molecule_info.h5"
		String cellranger_metrics_summary_csv = "~{cellranger_raw_data_path}/~{sample.sample_id}.metrics_summary.csv"

		if (cellranger_count_complete == "false") {
			call cellranger_count {
				input:
					sample_id = sample.sample_id,
					fastq_R1s = sample.fastq_R1s,
					fastq_R2s = sample.fastq_R2s,
					fastq_I1s = sample.fastq_I1s,
					fastq_I2s = sample.fastq_I2s,
					cellranger_reference_data = cellranger_reference_data,
					raw_data_path = cellranger_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File raw_counts_output = select_first([cellranger_count.raw_counts, cellranger_raw_counts]) #!FileCoercion
		File filtered_counts_output = select_first([cellranger_count.filtered_counts, cellranger_filtered_counts]) #!FileCoercion
		File molecule_info_output = select_first([cellranger_count.molecule_info, cellranger_molecule_info]) #!FileCoercion
		File metrics_summary_csv_output = select_first([cellranger_count.metrics_summary_csv, cellranger_metrics_summary_csv]) #!FileCoercion

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
					raw_counts = raw_counts_output,
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

		String preprocessed_adata_object = "~{adata_raw_data_path}/~{sample.sample_id}.adata_object.h5ad"

		if (adata_object_complete == "false") {
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

		File preprocessed_adata_object_output = select_first([counts_to_adata.adata_object, preprocessed_adata_object]) #!FileCoercion
	}

	output {
		# Sample list
		Array[Array[String]] project_sample_ids = project_sample_id

		# Cellranger
		Array[File] raw_counts = raw_counts_output #!FileCoercion
		Array[File] filtered_counts = filtered_counts_output #!FileCoercion
		Array[File] molecule_info = molecule_info_output #!FileCoercion
		Array[File] metrics_summary_csv = metrics_summary_csv_output #!FileCoercion

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
		Array[File] adata_object = preprocessed_adata_object_output #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] cellranger_count_output_files
		Array[String] remove_technical_artifacts_output_files
		Array[String] adata_object_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			cellranger_counts_file=$(echo "${output_files}" | cut -f 1)
			cellbender_counts_file=$(echo "${output_files}" | cut -f 2)
			adata_object_file=$(echo "${output_files}" | cut -f 3)

			if gsutil -u ~{billing_project} ls "${cellranger_counts_file}"; then
				if gsutil -u ~{billing_project} ls "${cellbender_counts_file}"; then
					if gsutil -u ~{billing_project} ls "${adata_object_file}"; then
						# If we find all outputs, don't rerun anything
						echo -e "true\ttrue\ttrue" >> sample_preprocessing_complete.tsv
					else
						# If we find cellranger and cellbender outputs, but don't find adata outputs, just rerun counts_to_adata
						echo -e "true\ttrue\tfalse" >> sample_preprocessing_complete.tsv
					fi
				else
					# If we find cellranger, but not cellbender outputs, it does not matter if adata objects exist, so run (or rerun) both cellbender and counts_to_adata
					echo -e "true\tfalse\tfalse" >> sample_preprocessing_complete.tsv
				fi
			else
				# If we don't find cellranger output, we must also need to run (or rerun) preprocessing
				echo -e "false\tfalse\tfalse" >> sample_preprocessing_complete.tsv
			fi
		done < <(paste ~{write_lines(cellranger_count_output_files)} ~{write_lines(remove_technical_artifacts_output_files)} ~{write_lines(adata_object_output_files)})
	>>>

	output {
		Array[Array[String]] sample_preprocessing_complete = read_tsv("sample_preprocessing_complete.tsv")
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
		zones: zones
	}
}

task cellranger_count {
	input {
		String sample_id

		Array[File] fastq_R1s
		Array[File] fastq_R2s
		Array[File] fastq_I1s
		Array[File] fastq_I2s

		File cellranger_reference_data

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 16
	Int mem_gb = 24
	Int disk_size = ceil((size(fastq_R1s, "GB") + size(fastq_R2s, "GB") + size(fastq_I1s, "GB") + size(fastq_I2s, "GB") + size(cellranger_reference_data, "GB")) * 4 + 50)

	command <<<
		set -euo pipefail

		# Unpack refdata
		mkdir cellranger_refdata
		tar \
			-zxvf ~{cellranger_reference_data} \
			-C cellranger_refdata \
			--strip-components 1

		# Ensure fastqs are in the same directory
		mkdir fastqs
		while read -r fastq || [[ -n "${fastq}" ]]; do
			if [[ -n "${fastq}" ]]; then
				validated_fastq_name=$(fix_fastq_names --fastq "${fastq}" --sample-id "~{sample_id}")
				if [[ -e "fastqs/${validated_fastq_name}" ]]; then
					echo "[ERROR] Something's gone wrong with fastq renaming; trying to create fastq [${validated_fastq_name}] but it already exists. Exiting."
					exit 1
				else
					ln -s "${fastq}" "fastqs/${validated_fastq_name}"
				fi
			fi
		done < <(cat \
			~{write_lines(fastq_R1s)} \
			~{write_lines(fastq_R2s)} \
			~{write_lines(fastq_I1s)} \
			~{write_lines(fastq_I2s)})

		cellranger --version

		/usr/bin/time \
		cellranger count \
			--id=~{sample_id} \
			--transcriptome="$(pwd)/cellranger_refdata" \
			--fastqs="$(pwd)/fastqs" \
			--localcores ~{threads} \
			--localmem ~{mem_gb - 4}

		# Rename outputs to include sample ID
		mv ~{sample_id}/outs/raw_feature_bc_matrix.h5 ~{sample_id}.raw_feature_bc_matrix.h5
		mv ~{sample_id}/outs/filtered_feature_bc_matrix.h5 ~{sample_id}.filtered_feature_bc_matrix.h5
		mv ~{sample_id}/outs/molecule_info.h5 ~{sample_id}.molecule_info.h5
		mv ~{sample_id}/outs/metrics_summary.csv ~{sample_id}.metrics_summary.csv

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.raw_feature_bc_matrix.h5" \
			-o "~{sample_id}.filtered_feature_bc_matrix.h5" \
			-o "~{sample_id}.molecule_info.h5" \
			-o "~{sample_id}.metrics_summary.csv"
	>>>

	output {
		String raw_counts = "~{raw_data_path}/~{sample_id}.raw_feature_bc_matrix.h5"
		String filtered_counts = "~{raw_data_path}/~{sample_id}.filtered_feature_bc_matrix.h5"
		String molecule_info = "~{raw_data_path}/~{sample_id}.molecule_info.h5"
		String metrics_summary_csv = "~{raw_data_path}/~{sample_id}.metrics_summary.csv"
	}

	runtime {
		docker: "~{container_registry}/cellranger:7.1.0"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
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

		nvidia-smi

		cellbender remove-background \
			--cuda \
			--input ~{raw_counts} \
			--output ~{sample_id}.cellbender. \
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
		memory: "64 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
		gpuType: "nvidia-tesla-t4"
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
			--adata-output ~{sample_id}.adata_object.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.adata_object.h5ad"
	>>>

	output {
		String adata_object = "~{raw_data_path}/~{sample_id}.adata_object.h5ad"
	}

	runtime {
		docker: "~{container_registry}/scvi:1.1.0"
		cpu: 2
		memory: "16 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}
