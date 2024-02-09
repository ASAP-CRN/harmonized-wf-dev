version 1.0

# Perform dataset integration, umap, clustering, and annotation steps

workflow cluster_data {
	input {
		String cohort_id
		File normalized_adata_object

		String scvi_latent_key

		Array[String] group_by_vars

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	call integrate_sample_data {
		input:
			cohort_id = cohort_id,
			normalized_adata_object = normalized_adata_object,
			scvi_latent_key = scvi_latent_key,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		File integrated_adata_object = integrate_sample_data.integrated_adata_object #!FileCoercion
		File scvi_model = integrate_sample_data.scvi_model #!FileCoercion
	}
}

task integrate_sample_data {
	input {
		String cohort_id
		File normalized_adata_object

		String scvi_latent_key

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int disk_size = ceil(size(normalized_adata_object, "GB") * 3 + 50)
	Int threads = 8
	Int mem_gb = threads * 2

	command <<<
		set -euo pipefail

		python integrate_scvi.py \
			--latent-key ~{scvi_latent_key} \
			--adata-input ~{normalized_adata_object} \
			--adata-output ~{cohort_id}.adata_object.scvi_integrated.h5ad.gz \
			--output-scvi ~{cohort_id}.scvi_model.pkl

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.adata_object.scvi_integrated.h5ad.gz" \
			-o "~{cohort_id}.scvi_model.pkl"
	>>>

	output {
		String integrated_adata_object = "~{raw_data_path}/~{cohort_id}.adata_object.scvi_integrated.h5ad.gz"
		String scvi_model = "~{raw_data_path}/~{cohort_id}.scvi_model.pkl"
	}

	runtime {
		docker: "~{container_registry}/scvi:1.0.4"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
		zones: zones
		cpuPlatform: "Intel Cascade Lake"
	}
}

task cluster_cells {
	input {
		String cohort_id
		File integrated_adata_object
		Int n_samples

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	String integrated_adata_object_basename = basename(integrated_adata_object, "_04.rds")
	Int threads = 2
	Int disk_size = ceil(size(integrated_adata_object, "GB") * 6 + 50)
	Int mem_gb = ceil(0.8 * n_samples + 50)

	command <<<
		set -euo pipefail

		# Find neighbors
		/usr/bin/time \
		Rscript /opt/scripts/main/find_neighbors.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--adata-object ~{integrated_adata_object} \
			--output-adata-object ~{integrated_adata_object_basename}_neighbors_05.rds

		# Run UMAP
		/usr/bin/time \
		Rscript /opt/scripts/main/umap.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--adata-object ~{integrated_adata_object_basename}_neighbors_05.rds \
			--output-adata-object ~{integrated_adata_object_basename}_neighbors_umap_06.rds

		# Cluster cells
		/usr/bin/time \
		Rscript /opt/scripts/main/clustering.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--adata-object ~{integrated_adata_object_basename}_neighbors_umap_06.rds \
			--clustering-algorithm ~{clustering_algorithm} \
			--clustering-resolution ~{clustering_resolution} \
			--cell-type-markers-list ~{cell_type_markers_list} \
			--output-cell-type-plot-prefix ~{cohort_id}.major_type_module_umap \
			--output-adata-object ~{integrated_adata_object_basename}_neighbors_umap_cluster_07.rds

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{integrated_adata_object_basename}_neighbors_05.rds" \
			-o "~{integrated_adata_object_basename}_neighbors_umap_06.rds" \
			-o "~{integrated_adata_object_basename}_neighbors_umap_cluster_07.rds" \
			-o "~{cohort_id}.major_type_module_umap.pdf" \
			-o "~{cohort_id}.major_type_module_umap.png"
	>>>

	output {
		String neighbors_adata_object = "~{raw_data_path}/~{integrated_adata_object_basename}_neighbors_05.rds"
		String umap_adata_object = "~{raw_data_path}/~{integrated_adata_object_basename}_neighbors_umap_06.rds"
		String cluster_adata_object = "~{raw_data_path}/~{integrated_adata_object_basename}_neighbors_umap_cluster_07.rds"
		String major_cell_type_plot_pdf = "~{raw_data_path}/~{cohort_id}.major_type_module_umap.pdf"
		String major_cell_type_plot_png = "~{raw_data_path}/~{cohort_id}.major_type_module_umap.png"
	}

	runtime {
		docker: "~{container_registry}/scvi:1.0.4"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
		zones: zones
		cpuPlatform: "Intel Cascade Lake"
	}
}

task annotate_clusters {
	input {
		String cohort_id
		File cluster_adata_object
		Int n_samples

		File cell_type_markers_list

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 2
	Int disk_size = ceil(size(cluster_adata_object, "GB") + size(cell_type_markers_list, "GB") * 2 + 20)
	Int mem_gb = ceil(1.3 * n_samples + 20)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/annotate_clusters.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--adata-object ~{cluster_adata_object} \
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
		docker: "~{container_registry}/scvi:1.0.4"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
		zones: zones
	}
}
