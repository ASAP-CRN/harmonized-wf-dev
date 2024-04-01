version 1.0

# Run steps in the cohort analysis

import "cluster_data/cluster_data.wdl" as ClusterData
import "../common/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids
		Array[File] preprocessed_adata_objects

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] preprocessing_output_file_paths = []

		Int n_top_genes

		String scvi_latent_key

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
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call filter_and_normalize {
		input:
			merged_adata_object = merge_and_plot_qc_metrics.merged_adata_object, #!FileCoercion
			n_top_genes = n_top_genes,
			container_registry = container_registry,
			zones = zones
	}

	call ClusterData.cluster_data {
		input:
			cohort_id = cohort_id,
			normalized_adata_object = select_first([filter_and_normalize.normalized_adata_object]), #!FileCoercion
			scvi_latent_key = scvi_latent_key,
			cell_type_markers_list = cell_type_markers_list,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call plot_groups_and_features {
		input:
			cohort_id = cohort_id,
			cell_annotated_adata_object = cluster_data.cell_annotated_adata_object,
			groups = groups,
			features = features,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	Array[String] cohort_analysis_intermediate_output_paths = flatten([
		select_all([filter_and_normalize.filtered_adata_object]),
		select_all([filter_and_normalize.normalized_adata_object]),
		[
			cluster_data.scvi_model_tar_gz
		]
	]) #!StringCoercion

	call UploadFinalOutputs.upload_final_outputs as upload_preprocess_files {
		input:
			output_file_paths = flatten([preprocessing_output_file_paths, cohort_analysis_intermediate_output_paths]),
			staging_data_buckets = staging_data_buckets,
			staging_data_path = "preprocess",
			billing_project = billing_project,
			zones = zones
	}

	Array[String] cohort_analysis_final_output_paths = flatten([
		[
			write_cohort_sample_list.cohort_sample_list
		],
		merge_and_plot_qc_metrics.qc_plots_png,
		[
			cluster_data.cell_types_csv,
			cluster_data.cell_annotated_adata_object,
			cluster_data.cell_annotated_metadata
		],
		[
			plot_groups_and_features.groups_umap_plot_png,
			plot_groups_and_features.features_umap_plot_png
		]
	]) #!StringCoercion

	call UploadFinalOutputs.upload_final_outputs as upload_cohort_analysis_files {
		input:
			output_file_paths = cohort_analysis_final_output_paths,
			staging_data_buckets = staging_data_buckets,
			staging_data_path = workflow_name,
			billing_project = billing_project,
			zones = zones
	}

	output {
		File cohort_sample_list = write_cohort_sample_list.cohort_sample_list #!FileCoercion

		# Merged adata objects and QC plots
		File merged_adata_object = merge_and_plot_qc_metrics.merged_adata_object
		Array[File] qc_plots_png = merge_and_plot_qc_metrics.qc_plots_png #!FileCoercion

		# Clustering output
		File integrated_adata_object = cluster_data.integrated_adata_object
		File scvi_model_tar_gz = cluster_data.scvi_model_tar_gz
		File umap_cluster_adata_object = cluster_data.umap_cluster_adata_object
		File cell_types_csv = cluster_data.cell_types_csv
		File cell_annotated_adata_object = cluster_data.cell_annotated_adata_object
		File cell_annotated_metadata = cluster_data.cell_annotated_metadata

		# Groups and features plots
		File groups_umap_plot_png = plot_groups_and_features.groups_umap_plot_png #!FileCoercion
		File features_umap_plot_png = plot_groups_and_features.features_umap_plot_png #!FileCoercion

		Array[File] preprocess_manifest_tsvs = upload_preprocess_files.manifests #!FileCoercion
		Array[File] cohort_analysis_manifest_tsvs = upload_cohort_analysis_files.manifests #!FileCoercion
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

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(preprocessed_adata_objects, "GB") * 2.4 + 20)
	Int disk_size = ceil(size(preprocessed_adata_objects, "GB") * 3 + 50)

	command <<<
		set -euo pipefail

		while read -r adata_objects || [[ -n "${adata_objects}" ]]; do 
			adata_path=$(realpath "${adata_objects}")
			sample=$(basename "${adata_path}" ".adata_object.h5ad")
			echo -e "${sample}\t${adata_path}" >> adata_samples_paths.tsv
		done < ~{write_lines(preprocessed_adata_objects)}

		python3 /opt/scripts/main/plot_qc_metrics.py \
			--adata-objects-fofn adata_samples_paths.tsv \
			--adata-output ~{cohort_id}.merged_adata_object.h5ad

		mv "plots/violin_n_genes_by_counts.png" "plots/~{cohort_id}.n_genes_by_counts.violin.png"
		mv "plots/violin_total_counts.png" "plots/~{cohort_id}.total_counts.violin.png"
		mv "plots/violin_pct_counts_mt.png" "plots/~{cohort_id}.pct_counts_mt.violin.png"
		mv "plots/violin_pct_counts_rb.png" "plots/~{cohort_id}.pct_counts_rb.violin.png"
		mv "plots/violin_doublet_score.png" "plots/~{cohort_id}.doublet_score.violin.png"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o plots/"~{cohort_id}.n_genes_by_counts.violin.png" \
			-o plots/"~{cohort_id}.total_counts.violin.png" \
			-o plots/"~{cohort_id}.pct_counts_mt.violin.png" \
			-o plots/"~{cohort_id}.pct_counts_rb.violin.png" \
			-o plots/"~{cohort_id}.doublet_score.violin.png"
	>>>

	output {
		File merged_adata_object = "~{cohort_id}.merged_adata_object.h5ad"

		Array[String] qc_plots_png = [
			"~{raw_data_path}/~{cohort_id}.n_genes_by_counts.violin.png",
			"~{raw_data_path}/~{cohort_id}.total_counts.violin.png",
			"~{raw_data_path}/~{cohort_id}.pct_counts_mt.violin.png",
			"~{raw_data_path}/~{cohort_id}.pct_counts_rb.violin.png",
			"~{raw_data_path}/~{cohort_id}.doublet_score.violin.png"
		]
	}

	runtime {
		docker: "~{container_registry}/scvi:1.1.0_1"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}

task filter_and_normalize {
	input {
		File merged_adata_object

		Int n_top_genes

		String container_registry
		String zones

		# Purposefully unset
		String? my_none
	}

	String merged_adata_object_basename = basename(merged_adata_object, ".h5ad")
	Int mem_gb = ceil(size(merged_adata_object, "GB") * 2.9 + 20)
	Int disk_size = ceil(size(merged_adata_object, "GB") * 4 + 20)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/main/filter.py \
			--adata-input ~{merged_adata_object} \
			--adata-output ~{merged_adata_object_basename}_filtered.h5ad

		# TODO see whether this is still required given the change to python
		# If any cells remain after filtering, the data is normalized and variable genes are identified
		if [[ -s "~{merged_adata_object_basename}_filtered.h5ad" ]]; then
			python3 /opt/scripts/main/process.py \
				--adata-input ~{merged_adata_object_basename}_filtered.h5ad \
				--batch-key "batch_id" \
				--adata-output ~{merged_adata_object_basename}_filtered_normalized.h5ad \
				--n-top-genes ~{n_top_genes}

			echo true > cells_remaining_post_filter.txt
		else
			echo false > cells_remaining_post_filter.txt
		fi
	>>>

	output {
		File? filtered_adata_object = if read_boolean("cells_remaining_post_filter.txt") then "~{merged_adata_object_basename}_filtered.h5ad" else my_none
		File? normalized_adata_object = if read_boolean("cells_remaining_post_filter.txt") then "~{merged_adata_object_basename}_filtered_normalized.h5ad" else my_none
	}

	runtime {
		docker: "~{container_registry}/scvi:1.1.0_1"
		cpu: 4
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}

task plot_groups_and_features {
	input {
		String cohort_id
		File cell_annotated_adata_object

		Array[String] groups
		Array[String] features

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(cell_annotated_adata_object, "GB") * 1.1 + 10)
	Int disk_size = ceil(size(cell_annotated_adata_object, "GB") * 4 + 20)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/main/plot_feats_and_groups.py \
			--adata-input ~{cell_annotated_adata_object} \
			--group ~{sep=',' groups} \
			--output-group-umap-plot-prefix "~{cohort_id}" \
			--feature ~{sep=',' features} \
			--output-feature-umap-plot-prefix "~{cohort_id}"

		mv "plots/umap~{cohort_id}_groups_umap.png" "plots/~{cohort_id}.groups.umap.png"
		mv "plots/umap~{cohort_id}_features_umap.png" "plots/~{cohort_id}.features.umap.png"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o plots/"~{cohort_id}.groups.umap.png" \
			-o plots/"~{cohort_id}.features.umap.png"
	>>>

	output {
		String groups_umap_plot_png = "~{raw_data_path}/~{cohort_id}.groups.umap.png"
		String features_umap_plot_png = "~{raw_data_path}/~{cohort_id}.features.umap.png"
	}

	runtime {
		docker: "~{container_registry}/scvi:1.1.0_1"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}
