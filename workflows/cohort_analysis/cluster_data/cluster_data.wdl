version 1.0

# Perform dataset integration, umap, clustering, and annotation steps

workflow cluster_data {
	input {
		String cohort_id
		Array[File] normalized_seurat_objects

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		String raw_data_path
		String curated_data_path
		String container_registry
	}

	call integrate_sample_data {
		input:
			cohort_id = cohort_id,
			normalized_seurat_objects = normalized_seurat_objects,
			raw_data_path = raw_data_path,
			container_registry = container_registry
	}

	call find_neighbors {
		input:
			integrated_seurat_object = integrate_sample_data.integrated_seurat_object, #!FileCoercion
			raw_data_path = raw_data_path,
			container_registry = container_registry
	}

	call run_umap {
		input:
			neighbors_seurat_object = find_neighbors.neighbors_seurat_object, #!FileCoercion
			raw_data_path = raw_data_path,
			container_registry = container_registry
	}

	call cluster_cells {
		input:
			cohort_id = cohort_id,
			umap_seurat_object = run_umap.umap_seurat_object, #!FileCoercion
			clustering_algorithm = clustering_algorithm,
			clustering_resolution = clustering_resolution,
			cell_type_markers_list = cell_type_markers_list,
			raw_data_path = raw_data_path,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	call annotate_clusters {
		input:
			cohort_id = cohort_id,
			cluster_seurat_object = cluster_cells.cluster_seurat_object, #!FileCoercion
			cell_type_markers_list = cell_type_markers_list,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	output {
		File cluster_seurat_object = cluster_cells.cluster_seurat_object #!FileCoercion
		File metadata = annotate_clusters.metadata #!FileCoercion
	}
}

task integrate_sample_data {
	input {
		String cohort_id
		Array[File] normalized_seurat_objects

		String raw_data_path
		String container_registry
	}

	# TODO Seems to only use ~4 threads; following snakemake for now
	Int threads = 8
	Int disk_size = ceil(size(normalized_seurat_objects[0], "GB") * length(normalized_seurat_objects) * 2 + 30)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/harmony.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-objects-fofn ~{write_lines(normalized_seurat_objects)} \
			--output-seurat-object ~{cohort_id}.seurat_object.harmony_integrated_04.rds

		# Upload outputs
		gsutil -m cp \
			~{cohort_id}.seurat_object.harmony_integrated_04.rds \
			~{raw_data_path}/
	>>>

	output {
		String integrated_seurat_object = "~{raw_data_path}/~{cohort_id}.seurat_object.harmony_integrated_04.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_1"
		cpu: threads
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task find_neighbors {
	input {
		File integrated_seurat_object

		String raw_data_path
		String container_registry
	}

	String integrated_seurat_object_basename = basename(integrated_seurat_object, "_04.rds")
	Int disk_size = ceil(size(integrated_seurat_object, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/find_neighbors.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{integrated_seurat_object} \
			--output-seurat-object ~{integrated_seurat_object_basename}_neighbors_05.rds

		# Upload outputs
		gsutil -m cp \
			~{integrated_seurat_object_basename}_neighbors_05.rds \
			~{raw_data_path}/
	>>>

	output {
		String neighbors_seurat_object = "~{raw_data_path}/~{integrated_seurat_object_basename}_neighbors_05.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_1"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task run_umap {
	input {
		File neighbors_seurat_object

		String raw_data_path
		String container_registry
	}

	String neighbors_seurat_object_basename = basename(neighbors_seurat_object, "_05.rds")
	Int disk_size = ceil(size(neighbors_seurat_object, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/umap.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{neighbors_seurat_object} \
			--output-seurat-object ~{neighbors_seurat_object_basename}_umap_06.rds

		# Upload outputs
		gsutil -m cp \
			~{neighbors_seurat_object_basename}_umap_06.rds \
			~{raw_data_path}/
	>>>

	output {
		String umap_seurat_object = "~{raw_data_path}/~{neighbors_seurat_object_basename}_umap_06.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_1"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task cluster_cells {
	input {
		String cohort_id
		File umap_seurat_object

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		String raw_data_path
		String curated_data_path
		String container_registry
	}

	# TODO only used 1 core
	Int threads = 8
	String umap_seurat_object_basename = basename(umap_seurat_object, "_06.rds")
	Int disk_size = ceil((size(umap_seurat_object, "GB") + size(cell_type_markers_list, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/clustering.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-object ~{umap_seurat_object} \
			--clustering-algorithm ~{clustering_algorithm} \
			--clustering-resolution ~{clustering_resolution} \
			--cell-type-markers-list ~{cell_type_markers_list} \
			--output-cell-type-plot ~{cohort_id}.major_type_module_umap.pdf \
			--output-seurat-object ~{umap_seurat_object_basename}_cluster_07.rds

		# Upload outputs
		gsutil -m cp \
			~{umap_seurat_object_basename}_cluster_07.rds \
			~{raw_data_path}/

		gsutil -m cp \
			~{cohort_id}.major_type_module_umap.pdf \
			~{curated_data_path}/
	>>>

	output {
		String major_cell_type_plot = "~{curated_data_path}/~{cohort_id}.major_type_module_umap.pdf"
		String cluster_seurat_object = "~{raw_data_path}/~{umap_seurat_object_basename}_cluster_07.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_1"
		cpu: threads
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

# TODO output file could be called ~{cohort_id}.annotated.csv ?
task annotate_clusters {
	input {
		String cohort_id
		File cluster_seurat_object

		File cell_type_markers_list

		String curated_data_path
		String container_registry
	}

	# TODO uses 2 cores
	Int threads = 8
	Int mem_gb = threads
	Int disk_size = ceil(size(cell_type_markers_list, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/annotate_clusters.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-object ~{cluster_seurat_object} \
			--cell-type-markers-list ~{cell_type_markers_list} \
			--output-metadata-file ~{cohort_id}.final_metadata.csv

		# Upload outputs
		gsutil -m cp \
			~{cohort_id}.final_metadata.csv \
			~{curated_data_path}/
	>>>

	output {
		String metadata = "~{curated_data_path}/~{cohort_id}.final_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}
