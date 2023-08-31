version 1.0

# Perform harmonization, umapping, and clustering steps

workflow cluster {
	input {
		String project_name
		Array[File] normalized_seurat_objects

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		String container_registry
	}

	call harmony {
		input:
			project_name = project_name,
			normalized_seurat_objects = normalized_seurat_objects,
			container_registry = container_registry
	}

	call neighbors {
		input:
			harmony_seurat_object = harmony.harmony_seurat_object,
			container_registry = container_registry
	}

	call umap {
		input:
			neighbors_seurat_object = neighbors.neighbors_seurat_object,
			container_registry = container_registry
	}

	call clustering {
		input:
			project_name = project_name,
			umap_seurat_object = umap.umap_seurat_object,
			clustering_algorithm = clustering_algorithm,
			clustering_resolution = clustering_resolution,
			cell_type_markers_list = cell_type_markers_list,
			container_registry = container_registry
	}

	call sctype {
		input:
			project_name = project_name,
			cluster_seurat_object = clustering.cluster_seurat_object,
			cell_type_markers_list = cell_type_markers_list,
			container_registry = container_registry
	}

	output {
		File cluster_seurat_object = clustering.cluster_seurat_object
		File metadata = sctype.metadata
	}
}

task harmony {
	input {
		String project_name
		Array[File] normalized_seurat_objects

		String container_registry
	}

	# TODO Seems to only use ~4 threads; following snakemake for now
	Int threads = 8
	Int disk_size = ceil(size(normalized_seurat_objects[0], "GB") * length(normalized_seurat_objects) * 2 + 30)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/harmony.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-objects-fofn ~{write_lines(normalized_seurat_objects)} \
			--output-seurat-object ~{project_name}.seurat_object.harmony_integrated_04.rds
	>>>

	output {
		File harmony_seurat_object = "~{project_name}.seurat_object.harmony_integrated_04.rds"
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

task neighbors {
	input {
		File harmony_seurat_object

		String container_registry
	}

	String harmony_seurat_object_basename = basename(harmony_seurat_object, "_04.rds")
	Int disk_size = ceil(size(harmony_seurat_object, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/find_neighbors.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{harmony_seurat_object} \
			--output-seurat-object ~{harmony_seurat_object_basename}_neighbors_05.rds
	>>>

	output {
		File neighbors_seurat_object = "~{harmony_seurat_object_basename}_neighbors_05.rds"
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

task umap {
	input {
		File neighbors_seurat_object

		String container_registry
	}

	String neighbors_seurat_object_basename = basename(neighbors_seurat_object, "_05.rds")
	Int disk_size = ceil(size(neighbors_seurat_object, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/umap.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{neighbors_seurat_object} \
			--output-seurat-object ~{neighbors_seurat_object_basename}_umap_06.rds
	>>>

	output {
		File umap_seurat_object = "~{neighbors_seurat_object_basename}_umap_06.rds"
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

task clustering {
	input {
		String project_name
		File umap_seurat_object

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		String container_registry
	}

	# TODO only used 1 core
	Int threads = 8
	String umap_seurat_object_basename = basename(umap_seurat_object, "_06.rds")
	Int disk_size = ceil((size(umap_seurat_object, "GB") + size(cell_type_markers_list, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/clustering.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-object ~{umap_seurat_object} \
			--clustering-algorithm ~{clustering_algorithm} \
			--clustering-resolution ~{clustering_resolution} \
			--cell-type-markers-list ~{cell_type_markers_list} \
			--output-cell-type-plot ~{project_name}.major_type_module_umap.pdf \
			--output-seurat-object ~{umap_seurat_object_basename}_cluster_07.rds
	>>>

	output {
		File major_cell_type_plot = "~{project_name}.major_type_module_umap.pdf"
		File cluster_seurat_object = "~{umap_seurat_object_basename}_cluster_07.rds"
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

task sctype {
	input {
		String project_name
		File cluster_seurat_object

		File cell_type_markers_list

		String container_registry
	}

	# TODO uses 2 cores
	Int threads = 8
	Int mem_gb = threads
	Int disk_size = ceil(size(cell_type_markers_list, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/annotate_clusters.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-object ~{cluster_seurat_object} \
			--cell-type-markers-list ~{cell_type_markers_list} \
			--output-metadata-file ~{project_name}.final_metadata.csv
	>>>

	output {
		File metadata = "~{project_name}.final_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}
