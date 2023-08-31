version 1.0

# Perform filtering steps

workflow filter {
	input {
		File preprocessed_seurat_object
		File unfiltered_metadata

		String container_registry
	}

	call filtering {
		input:
			preprocessed_seurat_object = preprocessed_seurat_object,
			unfiltered_metadata = unfiltered_metadata,
			container_registry = container_registry
	}

	call process {
		input:
			filtered_seurat_object = filtering.filtered_seurat_object,
			container_registry = container_registry
	}

	output {
		File normalized_seurat_object = process.normalized_seurat_object
	}
}

task filtering {
	input {
		File preprocessed_seurat_object
		File unfiltered_metadata

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
	>>>

	output {
		File filtered_seurat_object = "~{seurat_object_basename}_filtered_02.rds"
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

task process {
	input {
		File filtered_seurat_object

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
	>>>

	output {
		File normalized_seurat_object = "~{seurat_object_basename}_normalized_03.rds"
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
