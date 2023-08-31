version 1.0

# Get read counts and generate a preprocessed Seurat object

import "../structs.wdl"

workflow preprocess {
	input {
		Sample sample

		File cellranger_reference_data

		Float soup_rate

		String container_registry
	}

	call cellranger {
		input:
			sample_id = sample.sample_id,
			fastq_R1 = sample.fastq_R1,
			fastq_R2 = sample.fastq_R2,
			cellranger_reference_data = cellranger_reference_data,
			container_registry = container_registry
	}

	call counts_to_seurat {
		input:
			sample_id = sample.sample_id,
			batch = sample.batch,
			raw_counts = cellranger.raw_counts,
			filtered_counts = cellranger.filtered_counts,
			soup_rate = soup_rate,
			container_registry = container_registry
	}

	output {
		# Cellranger
		File raw_counts = cellranger.raw_counts
		File filtered_counts = cellranger.filtered_counts
		File molecule_info = cellranger.molecule_info
		File cellranger_metrics_csv = cellranger.metrics_csv

		File preprocessed_seurat_object = counts_to_seurat.preprocessed_seurat_object
	}
}

task cellranger {
	input {
		String sample_id

		File fastq_R1
		File fastq_R2

		File cellranger_reference_data

		String container_registry = "us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images"
	}

	Int threads = 16
	# TODO not sure this amount of RAM is necessary - cellranger docs claim it is, but test runs using ~15 GB (may be missing part of the process..?)
	Int mem_gb = threads * 8
	Int disk_size = ceil(size([fastq_R1, fastq_R2], "GB") * 2 + 30)

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
	>>>

	output {
		File raw_counts = "~{sample_id}.raw_feature_bc_matrix.h5"
		File filtered_counts = "~{sample_id}.filtered_feature_bc_matrix.h5"
		File molecule_info = "~{sample_id}.molecule_info.h5"
		File metrics_csv = "~{sample_id}.metrics_summary.csv"
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

		String container_registry
	}

	Int disk_size = ceil(size([raw_counts, filtered_counts], "GB") * 2 + 20)

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
			--output-seurat-object ~{sample_id}.seurat_object.preprocessed_01.rds
	>>>

	output {
		File preprocessed_seurat_object = "~{sample_id}.seurat_object.preprocessed_01.rds"
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
