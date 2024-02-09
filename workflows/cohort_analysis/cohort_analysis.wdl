version 1.0

# Run steps in the cohort analysis

import "cluster_data/cluster_data.wdl" as ClusterData
#import "../common/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids
		Array[File] preprocessed_adata_objects

		String scvi_latent_key

		Array[String] group_by_vars

		String clustering_method
		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		Array[String] groups
		Array[String] features

		String run_timestamp
		String raw_data_path_prefix
		Array[String] staging_data_buckets
		String billing_project
		String container_registry
		String zones
	}

	String workflow_name = "cohort_analysis"
	String workflow_version = "2.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version]]

	String raw_data_path = "~{raw_data_path_prefix}/~{workflow_name}/~{workflow_version}/~{run_timestamp}"

	Int n_samples = length(preprocessed_adata_objects)

	call write_cohort_sample_list {
		input:
			cohort_id = cohort_id,
			project_sample_ids = project_sample_ids,
			billing_project = billing_project,
			workflow_info = workflow_info,
			raw_data_path = raw_data_path,
			container_registry = container_registry,
			zones = zones
	}

	call merge_and_plot_qc_metrics {
		input:
			cohort_id = cohort_id,
			preprocessed_adata_objects = preprocessed_adata_objects,
			n_samples = n_samples,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call filter_and_normalize {
		input:
			merged_adata_object = merge_and_plot_qc_metrics.merged_adata_object, #!FileCoercion
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call ClusterData.cluster_data {
		input:
			cohort_id = cohort_id,
			normalized_adata_object = select_first([filter_and_normalize.normalized_adata_object]), #!FileCoercion
			scvi_latent_key = scvi_latent_key,
			group_by_vars = group_by_vars,
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
		File cohort_sample_list = write_cohort_sample_list.cohort_sample_list #!FileCoercion

		# Merged adata objects and QC plots
		File merged_adata_object = merge_and_plot_qc_metrics.merged_adata_object #!FileCoercion
		Array[File] qc_plots_png = merge_and_plot_qc_metrics.qc_plots_png #!FileCoercion

		# Clustering output
		File integrated_adata_object = cluster_data.integrated_adata_object #!FileCoercion
		File scvi_model = cluster_data.scvi_model #!FileCoercion
		File? umap_cluster_adata_object = cluster_data.umap_cluster_adata_object #!FileCoercion
		File? mde_cluster_adata_object = cluster_data.mde_cluster_adata_object #!FileCoercion
		File? major_cell_type_plot_pdf = cluster_data.major_cell_type_plot_pdf #!FileCoercion
		File? major_cell_type_plot_png = cluster_data.major_cell_type_plot_png #!FileCoercion
	}
}

# Upload the list of samples used for this cohort analysis to the output bucket
task write_cohort_sample_list {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	command <<<
		set -euo pipefail

		echo -e "project_id\tsample_id" > ~{cohort_id}.sample_list.tsv
		cat ~{write_tsv(project_sample_ids)} >> ~{cohort_id}.sample_list.tsv

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.sample_list.tsv"
	>>>

	output {
		String cohort_sample_list = "~{raw_data_path}/~{cohort_id}.sample_list.tsv"
	}

	runtime {
		docker: "~{container_registry}/util:1.1.0"
		cpu: 1
		memory: "1 GB"
		disks: "local-disk 10 HDD"
		preemptible: 3
		zones: zones
	}
}

task merge_and_plot_qc_metrics {
	input {
		String cohort_id
		Array[File] preprocessed_adata_objects
		Int n_samples

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 2
	Int disk_size = ceil(size(preprocessed_adata_objects, "GB") * 2 + 20)
	Int mem_gb = ceil(0.02 * n_samples + threads * 2 + 20)

	command <<<
		set -euo pipefail

		while read -r adata_objects || [[ -n "${adata_objects}" ]]; do 
			realpath "${adata_objects}" >> adata_objects_paths.txt
		done < ~{write_lines(preprocessed_adata_objects)}

		python plot_qc_metrics.py \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--adata-objects-fofn adata_objects_paths.txt \
			--project-name ~{cohort_id} \
			--adata-output ~{cohort_id}.merged_adata_object.h5ad.gz

		# TODO - double check file name and type for all plots
		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.merged_adata_object.h5ad.gz" \
			-o violin_n_genes_by_counts.png \
			-o violin_total_counts.png \
			-o violin_pct_counts_mt.png \
			-o violin_pct_counts_rb.png \
			-o violin_doublet_score.png
	>>>

	output {
		String merged_adata_object = "~{raw_data_path}/~{cohort_id}.merged_adata_object.h5ad.gz"

		Array[String] qc_plots_png = [
			"~{raw_data_path}/violin_n_genes_by_counts.png",
			"~{raw_data_path}/violin_total_counts.png",
			"~{raw_data_path}/violin_pct_counts_mt.png",
			"~{raw_data_path}/violin_pct_counts_rb.png",
			"~{raw_data_path}/violin_doublet_score.png"
		]
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

task filter_and_normalize {
	input {
		File merged_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones

		# Purposefully unset
		String? my_none
	}

	Int threads = 2
	String merged_adata_object_basename = basename(merged_adata_object, ".h5ad.gz")
	Int disk_size = ceil(size(merged_adata_object, "GB") * 4 + 20)

	command <<<
		set -euo pipefail

		python filter.py \
			--adata-input ~{merged_adata_object} \
			--adata-output ~{merged_adata_object_basename}_filtered.h5ad.gz

		# If any cells remain after filtering, the data is normalized and variable genes are identified
		if [[ -s "~{merged_adata_object_basename}_filtered.h5ad.gz" ]]; then
			python process.py \
				--working-dir "$(pwd)" \
				--adata-input ~{merged_adata_object_basename}_filtered.h5ad.gz \
				--adata-output ~{merged_adata_object_basename}_filtered_normalized.h5ad.gz

			upload_outputs \
				-b ~{billing_project} \
				-d ~{raw_data_path} \
				-i ~{write_tsv(workflow_info)} \
				-o "~{merged_adata_object_basename}_filtered.h5ad.gz" \
				-o "~{merged_adata_object_basename}_filtered_normalized.h5ad.gz"

			echo true > cells_remaining_post_filter.txt
		else
			echo false > cells_remaining_post_filter.txt
		fi
	>>>

	output {
		String? filtered_adata_object = if read_boolean("cells_remaining_post_filter.txt") then "~{raw_data_path}/~{merged_adata_object_basename}_filtered.h5ad.gz" else my_none
		String? normalized_adata_object = if read_boolean("cells_remaining_post_filter.txt") then "~{raw_data_path}/~{merged_adata_object_basename}_filtered_normalized.h5ad.gz" else my_none
	}

	runtime {
		docker: "~{container_registry}/scvi:1.0.4"
		cpu: threads
		memory: "12 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
		zones: zones
	}
}
