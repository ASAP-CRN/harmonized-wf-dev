version 1.0

# Get read counts and generate a preprocessed Seurat object

import "../structs.wdl"
import "../common/upload_final_outputs.wdl" as UploadFinalOutputs

workflow preprocess {
	input {
		String project_id
		Array[Sample] samples

		File cellranger_reference_data

		Float soup_rate

		Boolean regenerate_preprocessed_seurat_objects

		String run_timestamp
		String raw_data_path_prefix
		String staging_data_path_prefix
		String billing_project
		String container_registry
		Int multiome_container_revision
	}

	String workflow_name = "preprocess"

	String cellranger_task_version = "1.0.0"
	String counts_to_seurat_task_version = "1.0.0"

	# Used in the manifest; doesn't influence output locations or whether the task needs to be rerun
	String workflow_version = "~{cellranger_task_version}_~{counts_to_seurat_task_version}"

	String raw_data_path = "~{raw_data_path_prefix}/~{workflow_name}"
	String staging_data_path = "~{staging_data_path_prefix}/~{workflow_name}"

	scatter (sample_object in samples) {
		String cellranger_count_output = "~{raw_data_path}/cellranger/~{cellranger_task_version}/~{sample_object.sample_id}.raw_feature_bc_matrix.h5"
		String counts_to_seurat_output = "~{raw_data_path}/counts_to_seurat/~{counts_to_seurat_task_version}/~{sample_object.sample_id}.seurat_object.preprocessed_01.rds"
	}

	# For each sample, outputs an array of true/false: [cellranger_counts_complete, counts_to_seurat_complete]
	call check_output_files_exist {
		input:
			cellranger_count_output_files = cellranger_count_output,
			counts_to_seurat_output_files = counts_to_seurat_output,
			billing_project = billing_project
	}

	scatter (index in range(length(samples))) {
		Sample sample = samples[index]
		String cellranger_count_complete = check_output_files_exist.sample_preprocessing_complete[index][0]
		String counts_to_seurat_complete = check_output_files_exist.sample_preprocessing_complete[index][1]

		Array[String] project_sample_id = [project_id, sample.sample_id]

		String cellranger_raw_counts = "~{raw_data_path}/cellranger/~{cellranger_task_version}/~{sample.sample_id}.raw_feature_bc_matrix.h5"
		String cellranger_filtered_counts = "~{raw_data_path}/cellranger/~{cellranger_task_version}/~{sample.sample_id}.filtered_feature_bc_matrix.h5"
		String cellranger_molecule_info = "~{raw_data_path}/cellranger/~{cellranger_task_version}/~{sample.sample_id}.molecule_info.h5"
		String cellranger_metrics_csv = "~{raw_data_path}/cellranger/~{cellranger_task_version}/~{sample.sample_id}.metrics_summary.csv"
		String preprocessed_seurat_object = "~{raw_data_path}/counts_to_seurat/~{counts_to_seurat_task_version}/~{sample.sample_id}.seurat_object.preprocessed_01.rds"

		if (cellranger_count_complete == "false") {
			call cellranger_count {
				input:
					sample_id = sample.sample_id,
					fastq_R1 = sample.fastq_R1,
					fastq_R2 = sample.fastq_R2,
					fastq_I1 = sample.fastq_I1,
					fastq_I2 = sample.fastq_I2,
					cellranger_reference_data = cellranger_reference_data,
					raw_data_path = "~{raw_data_path}/cellranger/~{cellranger_task_version}",
					billing_project = billing_project,
					container_registry = container_registry
			}
		}

		File raw_counts_output = select_first([cellranger_count.raw_counts, cellranger_raw_counts]) #!FileCoercion
		File filtered_counts_output = select_first([cellranger_count.filtered_counts, cellranger_filtered_counts]) #!FileCoercion
		File molecule_info_output = select_first([cellranger_count.molecule_info, cellranger_molecule_info]) #!FileCoercion
		File metrics_csv_output = select_first([cellranger_count.metrics_csv, cellranger_metrics_csv]) #!FileCoercion

		if ((counts_to_seurat_complete == "false" && defined(sample.batch)) || regenerate_preprocessed_seurat_objects) {
			# Import counts and convert to a Seurat object
			call counts_to_seurat {
				input:
					sample_id = sample.sample_id,
					batch = select_first([sample.batch]),
					project_id = project_id,
					raw_counts = raw_counts_output, # !FileCoercion
					filtered_counts = filtered_counts_output, # !FileCoercion
					soup_rate = soup_rate,
					raw_data_path = "~{raw_data_path}/counts_to_seurat/~{counts_to_seurat_task_version}",
					billing_project = billing_project,
					container_registry = container_registry,
					multiome_container_revision = multiome_container_revision
			}
		}

		File seurat_object_output = select_first([counts_to_seurat.preprocessed_seurat_object, preprocessed_seurat_object]) #!FileCoercion
	}

	String preprocessing_manifest = "~{staging_data_path}/MANIFEST.tsv"
	Array[String] new_preprocessing_final_outputs = select_all(flatten([
		cellranger_count.raw_counts,
		cellranger_count.filtered_counts,
		cellranger_count.molecule_info,
		cellranger_count.metrics_csv,
		counts_to_seurat.preprocessed_seurat_object
	]))

	if (length(new_preprocessing_final_outputs) > 0) {
		call UploadFinalOutputs.upload_final_outputs {
			input:
				manifest_path = preprocessing_manifest,
				output_file_paths = new_preprocessing_final_outputs,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				run_timestamp = run_timestamp,
				staging_data_path = staging_data_path,
				billing_project = billing_project,
				container_registry = container_registry
		}
	}

	if (length(new_preprocessing_final_outputs) == 0) {
		call UploadFinalOutputs.sync_buckets {
			input:
				source_buckets = [
					"~{raw_data_path}/cellranger/~{cellranger_task_version}",
					"~{raw_data_path}/counts_to_seurat/~{counts_to_seurat_task_version}"
				],
				target_bucket = staging_data_path,
				billing_project = billing_project,
				container_registry = container_registry
		}
	}

	output {
		# Sample list
		Array[Array[String]] project_sample_ids = project_sample_id

		# Cellranger
		Array[File] raw_counts = raw_counts_output
		Array[File] filtered_counts = filtered_counts_output
		Array[File] molecule_info = molecule_info_output
		Array[File] metrics_csv = metrics_csv_output

		# Seurat counts
		Array[File] seurat_object = seurat_object_output

		File preprocessing_manifest_tsv = select_first([upload_final_outputs.updated_manifest, preprocessing_manifest]) #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] cellranger_count_output_files
		Array[String] counts_to_seurat_output_files

		String billing_project
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			counts_file=$(echo "${output_files}" | cut -f 1)
			seurat_object=$(echo "${output_files}" | cut -f 2)

			if gsutil -u ~{billing_project} ls "${seurat_object}"; then
				# If the seurat object exists, assume that the counts file does as well
				echo -e "true\ttrue" >> sample_preprocessing_complete.tsv
			else
				if gsutil -u ~{billing_project} ls "${counts_file}"; then
					echo -e "true\tfalse" >> sample_preprocessing_complete.tsv
				else
					echo -e "false\tfalse" >> sample_preprocessing_complete.tsv
				fi
			fi
		done < <(paste ~{write_lines(cellranger_count_output_files)} ~{write_lines(counts_to_seurat_output_files)})
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
	}
}

task cellranger_count {
	input {
		String sample_id

		File fastq_R1
		File fastq_R2
		File? fastq_I1
		File? fastq_I2

		File cellranger_reference_data

		String raw_data_path
		String billing_project
		String container_registry
	}

	Int threads = 16
	Int disk_size = ceil(size([fastq_R1, fastq_R2, cellranger_reference_data], "GB") * 4 + 50)
	Int mem_gb = 24

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
		ln -s ~{fastq_R1} ~{fastq_R2} ~{fastq_I1} ~{fastq_I2} fastqs/

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

		# Upload outputs
		gsutil -u ~{billing_project} -m cp \
			~{sample_id}.raw_feature_bc_matrix.h5 \
			~{sample_id}.filtered_feature_bc_matrix.h5 \
			~{sample_id}.molecule_info.h5 \
			~{sample_id}.metrics_summary.csv \
			~{raw_data_path}/
	>>>

	output {
		String raw_counts = "~{raw_data_path}/~{sample_id}.raw_feature_bc_matrix.h5"
		String filtered_counts = "~{raw_data_path}/~{sample_id}.filtered_feature_bc_matrix.h5"
		String molecule_info = "~{raw_data_path}/~{sample_id}.molecule_info.h5"
		String metrics_csv = "~{raw_data_path}/~{sample_id}.metrics_summary.csv"
	}

	runtime {
		docker: "~{container_registry}/cellranger:7.1.0"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
	}
}

task counts_to_seurat {
	input {
		String sample_id
		String batch
		String project_id

		File raw_counts
		File filtered_counts

		Float soup_rate

		String raw_data_path
		String billing_project
		String container_registry
		Int multiome_container_revision
	}

	Int disk_size = ceil(size([raw_counts, filtered_counts], "GB") * 2 + 20)
	# Memory scales with filtered_counts size
	Int mem_gb = ceil((size(filtered_counts, "GB") - 0.00132) / 0.001 + 10)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/preprocess.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--sample-id ~{sample_id} \
			--batch ~{batch} \
			--project ~{project_id} \
			--raw-counts ~{raw_counts} \
			--filtered-counts ~{filtered_counts} \
			--soup-rate ~{soup_rate} \
			--output-seurat-object ~{sample_id}.seurat_object.preprocessed_01.rds

		# Upload outputs
		gsutil -u ~{billing_project} -m cp \
			~{sample_id}.seurat_object.preprocessed_01.rds \
			~{raw_data_path}/
	>>>

	output {
		String preprocessed_seurat_object = "~{raw_data_path}/~{sample_id}.seurat_object.preprocessed_01.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_~{multiome_container_revision}"
		cpu: 4
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}
