version 1.0

# Perform dataset integration, umap, clustering, and annotation steps

workflow cluster_data {
	input {
		String cohort_id
		Array[File] normalized_seurat_objects
		Int n_samples

		Array[String] group_by_vars

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		Int multiome_container_revision
		String zones
	}

	call integrate_sample_data {
		input:
			cohort_id = cohort_id,
			normalized_seurat_objects = normalized_seurat_objects,
			n_samples = n_samples,
			group_by_vars = group_by_vars,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			multiome_container_revision = multiome_container_revision,
			zones = zones
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
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			multiome_container_revision = multiome_container_revision,
			zones = zones
	}

	call annotate_clusters {
		input:
			cohort_id = cohort_id,
			cluster_seurat_object = cluster_cells.cluster_seurat_object, #!FileCoercion
			n_samples = n_samples,
			cell_type_markers_list = cell_type_markers_list,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			multiome_container_revision = multiome_container_revision,
			zones = zones
	}

	output {
		File integrated_seurat_object = integrate_sample_data.integrated_seurat_object #!FileCoercion
		File neighbors_seurat_object = cluster_cells.neighbors_seurat_object #!FileCoercion
		File umap_seurat_object = cluster_cells.umap_seurat_object #!FileCoercion
		File cluster_seurat_object = cluster_cells.cluster_seurat_object #!FileCoercion
		File major_cell_type_plot_pdf = cluster_cells.major_cell_type_plot_pdf #!FileCoercion
		File major_cell_type_plot_png = cluster_cells.major_cell_type_plot_png #!FileCoercion
		File metadata = annotate_clusters.metadata #!FileCoercion
	}
}

task integrate_sample_data {
	input {
		String cohort_id
		Array[File] normalized_seurat_objects
		Int n_samples

		Array[String] group_by_vars

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		Int multiome_container_revision
		String zones
	}

	Int threads = 8
	Int disk_size = ceil(size(normalized_seurat_objects[0], "GB") * length(normalized_seurat_objects) * 2 + 30)
	Int mem_gb = ceil(n_samples * 1.6 + 50)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/harmony.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--group-by-vars ~{sep=' ' group_by_vars} \
			--seurat-objects-fofn ~{write_lines(normalized_seurat_objects)} \
			--output-seurat-object ~{cohort_id}.seurat_object.harmony_integrated_04.rds

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.seurat_object.harmony_integrated_04.rds"
	>>>

	output {
		String integrated_seurat_object = "~{raw_data_path}/~{cohort_id}.seurat_object.harmony_integrated_04.rds"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_~{multiome_container_revision}"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
		zones: zones
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
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		Int multiome_container_revision
		String zones
	}

	String integrated_seurat_object_basename = basename(integrated_seurat_object, "_04.rds")
	Int threads = 2
	Int disk_size = ceil(size(integrated_seurat_object, "GB") * 6 + 50)
	Int mem_gb = ceil(0.3 * n_samples + 15)

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
			--output-cell-type-plot-prefix ~{cohort_id}.major_type_module_umap \
			--output-seurat-object ~{integrated_seurat_object_basename}_neighbors_umap_cluster_07.rds

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{integrated_seurat_object_basename}_neighbors_05.rds" \
			-o "~{integrated_seurat_object_basename}_neighbors_umap_06.rds" \
			-o "~{integrated_seurat_object_basename}_neighbors_umap_cluster_07.rds" \
			-o "~{cohort_id}.major_type_module_umap.pdf" \
			-o "~{cohort_id}.major_type_module_umap.png"
	>>>

	output {
		String neighbors_seurat_object = "~{raw_data_path}/~{integrated_seurat_object_basename}_neighbors_05.rds"
		String umap_seurat_object = "~{raw_data_path}/~{integrated_seurat_object_basename}_neighbors_umap_06.rds"
		String cluster_seurat_object = "~{raw_data_path}/~{integrated_seurat_object_basename}_neighbors_umap_cluster_07.rds"
		String major_cell_type_plot_pdf = "~{raw_data_path}/~{cohort_id}.major_type_module_umap.pdf"
		String major_cell_type_plot_png = "~{raw_data_path}/~{cohort_id}.major_type_module_umap.png"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_~{multiome_container_revision}"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
		zones: zones
	}
}

task annotate_clusters {
	input {
		String cohort_id
		File cluster_seurat_object
		Int n_samples

		File cell_type_markers_list

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		Int multiome_container_revision
		String zones
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

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.final_metadata.csv"
	>>>

	output {
		String metadata = "~{raw_data_path}/~{cohort_id}.final_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_~{multiome_container_revision}"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
		zones: zones
	}
}
