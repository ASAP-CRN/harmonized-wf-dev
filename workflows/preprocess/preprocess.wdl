version 1.0

# Get read counts and generate a preprocessed Seurat object

import "../structs.wdl"

workflow preprocess {
	input {
		Sample sample

		File cellranger_reference_data

		Float soup_rate

		String raw_data_path_prefix
		String curated_data_path_prefix
		String container_registry
	}

	String workflow_name = "preprocess"
	String workflow_version = "0.0.1"

	String raw_data_path = "~{raw_data_path_prefix}/~{workflow_name}/~{workflow_version}"
	String curated_data_path = "~{curated_data_path_prefix}/~{workflow_name}/~{workflow_version}"

	String cellranger_raw_counts = "~{raw_data_path}/~{sample.sample_id}.raw_feature_bc_matrix.h5"
	String cellranger_filtered_counts = "~{raw_data_path}/~{sample.sample_id}.filtered_feature_bc_matrix.h5"
	String cellranger_molecule_info = "~{raw_data_path}/~{sample.sample_id}.molecule_info.h5"
	String cellranger_metrics_csv = "~{curated_data_path}/~{sample.sample_id}.metrics_summary.csv"
	String preprocessed_seurat_object = "~{raw_data_path}/~{sample.sample_id}.seurat_object.preprocessed_01.rds"

	# For this version of the pipeline, determine whether the sample has been run
	# (i.e. output files exist for this sample)
	call check_output_files_exist {
		input:
			expected_output_files = [
				cellranger_raw_counts,
				cellranger_filtered_counts,
				cellranger_molecule_info,
				cellranger_metrics_csv,
				preprocessed_seurat_object
			]
	}

	# If preprocessing pipeline outputs do not exist, run preprocessing steps
	if (! check_output_files_exist.outputs_exist) {
		call cellranger_count {
			input:
				sample_id = sample.sample_id,
				fastq_R1 = sample.fastq_R1,
				fastq_R2 = sample.fastq_R2,
				cellranger_reference_data = cellranger_reference_data,
				raw_data_path = raw_data_path,
				curated_data_path = curated_data_path,
				container_registry = container_registry
		}

		# Import counts and convert to a Seurat object
		call counts_to_seurat {
			input:
				sample_id = sample.sample_id,
				batch = sample.batch,
				raw_counts = cellranger_count.raw_counts, # !FileCoercion
				filtered_counts = cellranger_count.filtered_counts, # !FileCoercion
				soup_rate = soup_rate,
				raw_data_path = raw_data_path,
				container_registry = container_registry
		}
	}

	output {
		# Cellranger
		File raw_counts = select_first([cellranger_count.raw_counts, cellranger_raw_counts]) #!FileCoercion
		File filtered_counts = select_first([cellranger_count.filtered_counts, cellranger_filtered_counts]) #!FileCoercion
		File molecule_info = select_first([cellranger_count.molecule_info, cellranger_molecule_info]) #!FileCoercion
		File metrics_csv = select_first([cellranger_count.metrics_csv, cellranger_metrics_csv]) #!FileCoercion

		# Seurat counts
		File seurat_object = select_first([counts_to_seurat.preprocessed_seurat_object, preprocessed_seurat_object]) #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] expected_output_files
	}

	command <<<
		set -euo pipefail

		while read -r output_file || [[ -n "${output_file}" ]]; do
			if gsutil ls "${output_file}"; then
				echo "true" > output_exists.txt
			else
				echo "false" > output_exists.txt
				break
			fi
		done < ~{write_lines(expected_output_files)}
	>>>

	output {
		Boolean outputs_exist = read_boolean("output_exists.txt")
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim"
		cpu: 1
		memory: "1 GB"
		disks: "local-disk 10 HDD"
		preemptible: 3
	}
}

task cellranger_count {
	input {
		String sample_id

		File fastq_R1
		File fastq_R2

		File cellranger_reference_data

		String raw_data_path
		String curated_data_path
		String container_registry = "us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images"
	}

	Int threads = 16
	# TODO not sure this amount of RAM is necessary - cellranger docs claim it is, but test runs using ~15 GB (may be missing part of the process..?)
	Int mem_gb = threads * 8
	Int disk_size = ceil(size([fastq_R1, fastq_R2, cellranger_reference_data], "GB") * 4 + 50)

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
		ln -s ~{fastq_R1} ~{fastq_R2} fastqs/

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
		gsutil -m cp \
			~{sample_id}.raw_feature_bc_matrix.h5 \
			~{sample_id}.filtered_feature_bc_matrix.h5 \
			~{sample_id}.molecule_info.h5 \
			~{raw_data_path}/

		gsutil -m cp \
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
		String container_registry
	}

	Int disk_size = ceil(size([raw_counts, filtered_counts], "GB") * 2 + 20)

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
		gsutil -m cp \
			~{sample_id}.seurat_object.preprocessed_01.rds \
			~{raw_data_path}/
	>>>

	output {
		String preprocessed_seurat_object = "~{raw_data_path}/~{sample_id}.seurat_object.preprocessed_01.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: 8
		memory: "12 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}