version 1.0

# Get read counts and generate a preprocessed Seurat object

import "../structs.wdl"

workflow preprocess {
	input {
		String project_id
		Array[Sample] samples

		File cellranger_reference_data

		Float soup_rate

		String run_timestamp
		String raw_data_path_prefix
		String curated_data_path_prefix
		String billing_project
		String container_registry
	}

	String workflow_name = "preprocess"
	String workflow_version = "0.0.1"

	String raw_data_path = "~{raw_data_path_prefix}/~{workflow_name}/~{workflow_version}"
	String curated_data_path = "~{curated_data_path_prefix}/~{workflow_name}"

	scatter (sample_object in samples) {
		String cellranger_count_output = "~{raw_data_path}/~{sample_object.sample_id}.raw_feature_bc_matrix.h5"
		String counts_to_seurat_output = "~{raw_data_path}/~{sample_object.sample_id}.seurat_object.preprocessed_01.rds"
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

		String cellranger_raw_counts = "~{raw_data_path}/~{sample.sample_id}.raw_feature_bc_matrix.h5"
		String cellranger_filtered_counts = "~{raw_data_path}/~{sample.sample_id}.filtered_feature_bc_matrix.h5"
		String cellranger_molecule_info = "~{raw_data_path}/~{sample.sample_id}.molecule_info.h5"
		String cellranger_metrics_csv = "~{raw_data_path}/~{sample.sample_id}.metrics_summary.csv"
		String preprocessed_seurat_object = "~{raw_data_path}/~{sample.sample_id}.seurat_object.preprocessed_01.rds"

		if (cellranger_count_complete == "false") {
			call cellranger_count {
				input:
					sample_id = sample.sample_id,
					fastq_R1 = sample.fastq_R1,
					fastq_R2 = sample.fastq_R2,
					fastq_I1 = sample.fastq_I1,
					fastq_I2 = sample.fastq_I2,
					cellranger_reference_data = cellranger_reference_data,
					raw_data_path = raw_data_path,
					curated_data_path = curated_data_path,
					billing_project = billing_project,
					container_registry = container_registry
			}
		}

		File raw_counts_output = select_first([cellranger_count.raw_counts, cellranger_raw_counts]) #!FileCoercion
		File filtered_counts_output = select_first([cellranger_count.filtered_counts, cellranger_filtered_counts]) #!FileCoercion
		File molecule_info_output = select_first([cellranger_count.molecule_info, cellranger_molecule_info]) #!FileCoercion
		File metrics_csv_output = select_first([cellranger_count.metrics_csv, cellranger_metrics_csv]) #!FileCoercion

		if (counts_to_seurat_complete == "false" && defined(sample.batch)) {
			# Import counts and convert to a Seurat object
			call counts_to_seurat {
				input:
					sample_id = sample.sample_id,
					batch = select_first([sample.batch]),
					raw_counts = raw_counts_output, # !FileCoercion
					filtered_counts = filtered_counts_output, # !FileCoercion
					soup_rate = soup_rate,
					raw_data_path = raw_data_path,
					billing_project = billing_project,
					container_registry = container_registry
			}
		}

		File seurat_object_output = select_first([counts_to_seurat.preprocessed_seurat_object, preprocessed_seurat_object]) #!FileCoercion
	}

	String preprocessing_manifest = "~{curated_data_path}/MANIFEST.tsv"
	Array[String] new_preprocessing_output = select_all(cellranger_count.metrics_csv)
	if (length(new_preprocessing_output) > 0) {
		call upload_final_outputs {
			input:
				manifest_path = preprocessing_manifest,
				output_file_paths = new_preprocessing_output,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				run_timestamp = run_timestamp,
				curated_data_path = curated_data_path,
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
		String curated_data_path
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
			~{raw_data_path}/

		gsutil -u ~{billing_project} -m cp \
			~{sample_id}.metrics_summary.csv \
			~{curated_data_path}/
	>>>

	output {
		String raw_counts = "~{raw_data_path}/~{sample_id}.raw_feature_bc_matrix.h5"
		String filtered_counts = "~{raw_data_path}/~{sample_id}.filtered_feature_bc_matrix.h5"
		String molecule_info = "~{raw_data_path}/~{sample_id}.molecule_info.h5"
		String metrics_csv = "~{curated_data_path}/~{sample_id}.metrics_summary.csv"
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

		File raw_counts
		File filtered_counts

		Float soup_rate

		String raw_data_path
		String billing_project
		String container_registry
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
		docker: "~{container_registry}/multiome:4a7fd84_4"
		cpu: 4
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task upload_final_outputs {
	input {
		String manifest_path
		Array[String] output_file_paths

		String workflow_name
		String workflow_version
		String run_timestamp

		String curated_data_path
		String billing_project
		String container_registry
	}

	command <<<
		set -euo pipefail

		echo -e "filename\tworkflow\tworkflow_version\ttimestamp" > new_files.manifest.tsv
		while read -r output_file || [[ -n "${output_file}" ]]; do
			gsutil -u ~{billing_project} -m cp "${output_file}" ~{curated_data_path}

			echo -e "$(basename "${output_file}")\t~{workflow_name}\t~{workflow_version}\t~{run_timestamp}" >> new_files.manifest.tsv
		done < ~{write_lines(output_file_paths)}

		if gsutil ls ~{manifest_path}; then
			# If a manifest already exists, merge the previous and updated manifests
			#   and replace the existing manifest with the updated one
			gsutil -u ~{billing_project} -m cp ~{manifest_path} previous_manifest.tsv

			update_manifest.py \
				--previous-manifest previous_manifest.tsv \
				--new-files-manifest new_files.manifest.tsv \
				--updated-manifest updated_manifest.tsv

			gsutil -u ~{billing_project} -m cp updated_manifest.tsv ~{manifest_path}
		else
			# If a manifest does not exist, create one (containing info about the files just uploaded)
			gsutil -u ~{billing_project} -m cp new_files_manifest.tsv ~{manifest_path}
		fi
	>>>

	output {
		String updated_manifest = manifest_path
	}

	runtime {
		docker: "~{container_registry}/util:0.0.1"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
	}
}
