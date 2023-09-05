version 1.0

# Perform filtering steps

workflow filter_data {
	input {
		File preprocessed_seurat_object
		File unfiltered_metadata

		String raw_data_path
		String container_registry
	}

	# Filter sample data by nCount_RNA, nFeature_RNA, percent.mt, percent.rb; remove doublets
	call apply_filters {
		input:
			preprocessed_seurat_object = preprocessed_seurat_object,
			unfiltered_metadata = unfiltered_metadata,
			raw_data_path = raw_data_path,
			container_registry = container_registry
	}

	call normalize_scale_data {
		input:
			filtered_seurat_object = apply_filters.filtered_seurat_object, #!FileCoercion
			raw_data_path = raw_data_path,
			container_registry = container_registry
	}

	output {
		File normalized_seurat_object = normalize_scale_data.normalized_seurat_object #!FileCoercion
	}
}

task apply_filters {
	input {
		File preprocessed_seurat_object
		File unfiltered_metadata

		String raw_data_path
		String container_registry
	}

	String seurat_object_basename = basename(preprocessed_seurat_object, "_01.rds")
	Int disk_size = ceil(size(preprocessed_seurat_object, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/filter.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{preprocessed_seurat_object} \
			--metadata ~{unfiltered_metadata} \
			--output-seurat-object ~{seurat_object_basename}_filtered_02.rds

		# Upload outputs
		gsutil -m cp \
			~{seurat_object_basename}_filtered_02.rds \
			~{raw_data_path}/
	>>>

	output {
		String filtered_seurat_object = "~{raw_data_path}/~{seurat_object_basename}_filtered_02.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task normalize_scale_data {
	input {
		File filtered_seurat_object

		String raw_data_path
		String container_registry
	}

	Int threads = 2
	String seurat_object_basename = basename(filtered_seurat_object, "_02.rds")
	Int disk_size = ceil(size(filtered_seurat_object, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/process.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-object ~{filtered_seurat_object} \
			--output-seurat-object ~{seurat_object_basename}_normalized_03.rds

		# Upload outputs
		gsutil -m cp \
			~{seurat_object_basename}_normalized_03.rds \
			~{raw_data_path}/
	>>>

	output {
		String normalized_seurat_object = "~{raw_data_path}/~{seurat_object_basename}_normalized_03.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: threads
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}
