version 1.0

# Perform dataset integration, umap, clustering, and annotation steps

workflow cluster_data {
	input {
		String cohort_id
		File normalized_adata_object

		String scvi_latent_key

		Array[String] group_by_vars

		String clustering_method
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

	call cluster_cells {
		input:
			cohort_id = cohort_id,
			integrated_adata_object = integrate_sample_data.integrated_adata_object, #!FileCoercion
			scvi_latent_key =scvi_latent_key,
			clustering_method = clustering_method,
			clustering_algorithm = clustering_algorithm,
			clustering_resolution = clustering_resolution,
			cell_type_markers_list = cell_type_markers_list,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		File integrated_adata_object = integrate_sample_data.integrated_adata_object #!FileCoercion
		File scvi_model = integrate_sample_data.scvi_model #!FileCoercion
		File? umap_cluster_adata_object = cluster_cells.umap_cluster_adata_object #!FileCoercion, #!UnnecessaryQuantifier
		File? mde_cluster_adata_object = cluster_cells.mde_cluster_adata_object #!FileCoercion, #!UnnecessaryQuantifier
		File? major_cell_type_plot_pdf = cluster_cells.major_cell_type_plot_pdf #!FileCoercion, #!UnnecessaryQuantifier
		File? major_cell_type_plot_png = cluster_cells.major_cell_type_plot_png #!FileCoercion, #!UnnecessaryQuantifier
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
	}
}

task cluster_cells {
	input {
		String cohort_id
		File integrated_adata_object

		String scvi_latent_key
		String clustering_method
		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	String integrated_adata_object_basename = basename(integrated_adata_object, ".h5ad.gz")
	Int disk_size = ceil(size(integrated_adata_object, "GB") * 6 + 50)
	Int threads = 2
	Int mem_gb = threads * 2

	# TODO - double check that empty gpuCount and gpuType does not launch GPU for clustering_umap
	String gpu_type = if clustering_method == "mde" then "nvidia-tesla-k80" else ""
	Int gpu_count = if clustering_method == "mde" then 2 else ""

	command <<<
		set -euo pipefail

		if [[ ~{clustering_method} = "umap" ]]; then
			# TODO - script doesn't use cell_type_markers_list (incomplete)
			python clustering_umap.py \
				--working-dir "$(pwd)" \
				--script-dir /opt/scripts \
				--threads ~{threads} \
				--adata-input ~{integrated_adata_object} \
				--adata-output ~{integrated_adata_object_basename}.umap_cluster.h5ad.gz \
				--latent-key ~{scvi_latent_key}

			upload_outputs \
				-b ~{billing_project} \
				-d ~{raw_data_path} \
				-i ~{write_tsv(workflow_info)} \
				-o "~{integrated_adata_object_basename}.umap_cluster.h5ad.gz"

		elif [[ ~{clustering_method} = "mde" ]]; then
			# Note: mde is super fast and efficient on a GPU
			python clustering_mde.py \
				--working-dir "$(pwd)" \
				--script-dir /opt/scripts \
				--threads ~{threads} \
				--adata-object ~{integrated_adata_object} \
				--clustering-algorithm ~{clustering_algorithm} \
				--clustering-resolution ~{clustering_resolution} \
				--cell-type-markers-list ~{cell_type_markers_list} \
				--output-cell-type-plot-prefix ~{cohort_id}.major_type_module_umap \
				--adata-output ~{integrated_adata_object_basename}.mde_cluster.h5ad.gz \
				--latent-key ~{scvi_latent_key}

			# TODO - this will fail because plots are not outputted by script (incomplete)
			#upload_outputs \
			#	-b ~{billing_project} \
			#	-d ~{raw_data_path} \
			#	-i ~{write_tsv(workflow_info)} \
			#	-o "~{integrated_adata_object_basename}.mde_cluster.h5ad.gz" \
			#	-o "~{cohort_id}.major_type_module_umap.pdf" \
			#	-o "~{cohort_id}.major_type_module_umap.png"

		else
			echo "[ERROR] Choose a valid clustering method; options are 'umap' or 'mde'"
			exit 1
		fi
	>>>

	output {
		String umap_cluster_adata_object = "~{raw_data_path}/~{integrated_adata_object_basename}.umap_cluster.h5ad.gz"
		String mde_cluster_adata_object = "~{raw_data_path}/~{integrated_adata_object_basename}.mde_cluster.h5ad.gz"
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
		gpuType: gpu_type
		gpuCount: gpu_count
	}
}
