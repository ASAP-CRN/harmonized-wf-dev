version 1.0

# Perform dataset integration, umap, clustering, and annotation steps

workflow cluster_data {
	input {
		String cohort_id
		Array[File] normalized_seurat_objects
		Int n_samples

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
			n_samples = n_samples,
			raw_data_path = raw_data_path,
			container_registry = container_registry
	}

	# Find neighbors, run umap, and cluster cells
	call cluster_cells {
		input:
			cohort_id = cohort_id,
			integrated_seurat_object = integrate_sample_data.integrated_seurat_object, #!FileCoercion
			n_samples = n_samples,
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
			n_samples = n_samples,
			cell_type_markers_list = cell_type_markers_list,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	output {
		File cluster_seurat_object = cluster_cells.cluster_seurat_object #!FileCoercion
		File major_cell_type_plot = cluster_cells.major_cell_type_plot #!FileCoercion
		File metadata = annotate_clusters.metadata #!FileCoercion
	}
}

task integrate_sample_data {
	input {
		String cohort_id
		Array[File] normalized_seurat_objects
		Int n_samples

		String raw_data_path
		String container_registry
	}

	Int threads = 8
	Int disk_size = ceil(size(normalized_seurat_objects[0], "GB") * length(normalized_seurat_objects) * 2 + 30)
	Int mem_gb = ceil(0.4 * n_samples + 10)

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
		docker: "~{container_registry}/multiome:4a7fd84_4"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task cluster_cells {
	input {
		String cohort_id
		File integrated_seurat_object
		Int n_samples

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		String raw_data_path
		String curated_data_path
		String container_registry
	}

	String integrated_seurat_object_basename = basename(integrated_seurat_object, "_04.rds")
	Int threads = 2
	Int disk_size = ceil(size(integrated_seurat_object, "GB") * 6 + 50)
	Int mem_gb = ceil(0.2 * n_samples + 5)

	command <<<
		set -euo pipefail

		# Find neighbors
		/usr/bin/time \
		Rscript /opt/scripts/main/find_neighbors.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{integrated_seurat_object} \
			--output-seurat-object ~{integrated_seurat_object_basename}_neighbors_05.rds

		# Run UMAP
		/usr/bin/time \
		Rscript /opt/scripts/main/umap.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{integrated_seurat_object_basename}_neighbors_05.rds \
			--output-seurat-object ~{integrated_seurat_object_basename}_neighbors_umap_06.rds

		# Cluster cells
		/usr/bin/time \
		Rscript /opt/scripts/main/clustering.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-object ~{integrated_seurat_object_basename}_neighbors_umap_06.rds \
			--clustering-algorithm ~{clustering_algorithm} \
			--clustering-resolution ~{clustering_resolution} \
			--cell-type-markers-list ~{cell_type_markers_list} \
			--output-cell-type-plot ~{cohort_id}.major_type_module_umap.pdf \
			--output-seurat-object ~{integrated_seurat_object_basename}_neighbors_umap_cluster_07.rds

		# Upload outputs
		gsutil -m cp \
			~{integrated_seurat_object_basename}_neighbors_05.rds \
			~{integrated_seurat_object_basename}_neighbors_umap_06.rds \
			~{integrated_seurat_object_basename}_neighbors_umap_cluster_07.rds \
			~{raw_data_path}/

		gsutil -m cp \
			~{cohort_id}.major_type_module_umap.pdf \
			~{curated_data_path}/
	>>>

	output {
		String neighbors_seurat_object = "~{raw_data_path}/~{integrated_seurat_object_basename}_neighbors_05.rds"
		String umap_seurat_object = "~{raw_data_path}/~{integrated_seurat_object_basename}_neighbors_umap_06.rds"
		String cluster_seurat_object = "~{raw_data_path}/~{integrated_seurat_object_basename}_neighbors_umap_cluster_07.rds"
		String major_cell_type_plot = "~{curated_data_path}/~{cohort_id}.major_type_module_umap.pdf"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_4"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task annotate_clusters {
	input {
		String cohort_id
		File cluster_seurat_object
		Int n_samples

		File cell_type_markers_list

		String curated_data_path
		String container_registry
	}

	Int threads = 2
	Int disk_size = ceil(size(cluster_seurat_object, "GB") + size(cell_type_markers_list, "GB") * 2 + 20)
	Int mem_gb = ceil(1.3 * n_samples + 20)

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
		docker: "~{container_registry}/multiome:4a7fd84_4"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}
