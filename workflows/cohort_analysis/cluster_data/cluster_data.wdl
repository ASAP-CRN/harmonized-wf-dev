version 1.0

# Perform dataset integration, umap, clustering, and annotation steps

workflow cluster_data {
	input {
		String cohort_id
		File normalized_adata_object

		String scvi_latent_key

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
			integrated_adata_object = integrate_sample_data.integrated_adata_object, #!FileCoercion
			scvi_latent_key =scvi_latent_key,
			cell_type_markers_list = cell_type_markers_list,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call annotate_cells {
		input:
			cohort_id = cohort_id,
			cluster_adata_object = cluster_cells.umap_cluster_adata_object, #!FileCoercion
			cell_type_markers_list = cell_type_markers_list,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		File integrated_adata_object = integrate_sample_data.integrated_adata_object #!FileCoercion
		File scvi_model_tar_gz = integrate_sample_data.scvi_model_tar_gz #!FileCoercion
		File umap_cluster_adata_object = cluster_cells.umap_cluster_adata_object #!FileCoercion
		File cell_types_csv = annotate_cells.cell_types_csv #!FileCoercion
		File cell_annotated_adata_object = annotate_cells.cell_annotated_adata_object #!FileCoercion
		File cell_annotated_metadata = annotate_cells.cell_annotated_metadata #!FileCoercion
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

	command <<<
		set -euo pipefail

		nvidia-smi

		python3 /opt/scripts/main/integrate_scvi.py \
			--latent-key ~{scvi_latent_key} \
			--batch-key "batch_id" \
			--adata-input ~{normalized_adata_object} \
			--adata-output ~{cohort_id}.adata_object.scvi_integrated.h5ad \
			--output-scvi-dir ~{cohort_id}_scvi_model

		# Model name cannot be changed because scvi models serialization expects a path containing a model.pt object
		tar -czvf "~{cohort_id}_scvi_model.tar.gz" "~{cohort_id}_scvi_model"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.adata_object.scvi_integrated.h5ad" \
			-o "~{cohort_id}_scvi_model.tar.gz"
	>>>

	output {
		String integrated_adata_object = "~{raw_data_path}/~{cohort_id}.adata_object.scvi_integrated.h5ad"
		String scvi_model_tar_gz = "~{raw_data_path}/~{cohort_id}_scvi_model.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/scvi:1.1.0"
		cpu: 2
		memory: "100 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
		gpuType: "nvidia-tesla-t4"
		gpuCount: 1
		nvidiaDriverVersion: "530.30.02" #!UnknownRuntimeKey
	}
}

task cluster_cells {
	input {
		File integrated_adata_object

		String scvi_latent_key
		File cell_type_markers_list

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	String integrated_adata_object_basename = basename(integrated_adata_object, ".h5ad")
	Int disk_size = ceil(size([integrated_adata_object, cell_type_markers_list], "GB") * 6 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/main/clustering_umap.py \
			--adata-input ~{integrated_adata_object} \
			--adata-output ~{integrated_adata_object_basename}.umap_cluster.h5ad \
			--latent-key ~{scvi_latent_key}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{integrated_adata_object_basename}.umap_cluster.h5ad"
	>>>

	output {
		String umap_cluster_adata_object = "~{raw_data_path}/~{integrated_adata_object_basename}.umap_cluster.h5ad"
	}

	runtime {
		docker: "~{container_registry}/scvi:1.1.0"
		cpu: 16
		memory: "150 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}

task annotate_cells {
	input {
		String cohort_id
		File cluster_adata_object

		File cell_type_markers_list

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	String cluster_adata_object_basename = basename(cluster_adata_object, ".h5ad")
	Int disk_size = ceil(size([cluster_adata_object, cell_type_markers_list], "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		nvidia-smi

		# This is not annotating the "clusters" but rather, the cells based on marker gene expression
		python3 /opt/scripts/main/annotate_cells.py \
			--adata-input ~{cluster_adata_object} \
			--marker-genes ~{cell_type_markers_list} \
			--batch-key "batch_id" \
			--output-cell-types-file ~{cohort_id}.cell_types.csv \
			--adata-output ~{cluster_adata_object_basename}.annotate_cells.h5ad \
			--output-metadata-file ~{cohort_id}.annotate_cells.metadata.csv

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.cell_types.csv" \
			-o "~{cluster_adata_object_basename}.annotate_cells.h5ad" \
			-o "~{cohort_id}.annotate_cells.metadata.csv"
	>>>

	output {
		String cell_types_csv = "~{raw_data_path}/~{cohort_id}.cell_types.csv"
		String cell_annotated_adata_object = "~{raw_data_path}/~{cluster_adata_object_basename}.annotate_cells.h5ad"
		String cell_annotated_metadata = "~{raw_data_path}/~{cohort_id}.annotate_cells.metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/scvi:1.1.0"
		cpu: 2
		memory: "100 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
		gpuType: "nvidia-tesla-t4"
		gpuCount: 1
		nvidiaDriverVersion: "530.30.02" #!UnknownRuntimeKey
	}
}
