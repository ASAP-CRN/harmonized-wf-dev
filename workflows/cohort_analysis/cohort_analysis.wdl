version 1.0

# Run steps in the cohort analysis

import "run_quality_control/run_quality_control.wdl" as RunQualityControl
import "cluster_data/cluster_data.wdl" as ClusterData
import "../common/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids
		Array[File] preprocessed_seurat_objects

		Array[String] group_by_vars

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		Array[String] groups
		Array[String] features

		String run_timestamp
		String raw_data_path_prefix
		String curated_data_path_prefix
		String billing_project
		String container_registry
	}

	String workflow_name = "cohort_analysis"
	String workflow_version = "0.1.0"

	String raw_data_path = "~{raw_data_path_prefix}/~{workflow_name}/~{workflow_version}/~{run_timestamp}"
	String curated_data_path = "~{curated_data_path_prefix}/~{workflow_name}/"

	Int n_samples = length(preprocessed_seurat_objects)

	call write_cohort_sample_list {
		input:
			cohort_id = cohort_id,
			project_sample_ids = project_sample_ids,
			billing_project = billing_project,
			raw_data_path = raw_data_path
	}

	call RunQualityControl.run_quality_control {
		input:
			cohort_id = cohort_id,
			preprocessed_seurat_objects = preprocessed_seurat_objects,
			n_samples = n_samples,
			raw_data_path = raw_data_path,
			billing_project = billing_project,
			container_registry = container_registry
	}

	scatter (preprocessed_seurat_object in preprocessed_seurat_objects) {
		# Filter sample data by nCount_RNA, nFeature_RNA, percent.mt, percent.rb; remove doublets
		# If cells remain after filtering, normalize and scale data
		call filter_and_normalize {
			input:
				preprocessed_seurat_object = preprocessed_seurat_object,
				unfiltered_metadata = run_quality_control.unfiltered_metadata,
				raw_data_path = raw_data_path,
				billing_project = billing_project,
				container_registry = container_registry
		}
	}

	call ClusterData.cluster_data {
		input:
			cohort_id = cohort_id,
			normalized_seurat_objects = select_all(filter_and_normalize.normalized_seurat_object), #!FileCoercion
			n_samples = n_samples,
			group_by_vars = group_by_vars,
			clustering_algorithm = clustering_algorithm,
			clustering_resolution = clustering_resolution,
			cell_type_markers_list = cell_type_markers_list,
			raw_data_path = raw_data_path,
			billing_project = billing_project,
			container_registry = container_registry
	}

	call plot_groups_and_features {
		input:
			cohort_id = cohort_id,
			metadata = cluster_data.metadata,
			groups = groups,
			features = features,
			raw_data_path = raw_data_path,
			billing_project = billing_project,
			container_registry = container_registry
	}

	Array[String] cohort_analysis_final_outputs = [
		write_cohort_sample_list.cohort_sample_list,
		run_quality_control.qc_violin_plots,
		run_quality_control.qc_umis_genes_plot,
		cluster_data.major_cell_type_plot,
		cluster_data.metadata,
		plot_groups_and_features.group_umap_plots,
		plot_groups_and_features.feature_umap_plots
	] #!StringCoercion

	String cohort_analysis_manifest = "~{curated_data_path}/MANIFEST.tsv"
	call UploadFinalOutputs.upload_final_outputs {
		input:
			manifest_path = cohort_analysis_manifest,
			output_file_paths = cohort_analysis_final_outputs,
			workflow_name = workflow_name,
			workflow_version = workflow_version,
			run_timestamp = run_timestamp,
			curated_data_path = curated_data_path,
			billing_project = billing_project,
			container_registry = container_registry
	}

	output {
		File cohort_sample_list = write_cohort_sample_list.cohort_sample_list #!FileCoercion

		# QC plots
		File qc_violin_plots = run_quality_control.qc_violin_plots
		File qc_umis_genes_plot = run_quality_control.qc_umis_genes_plot

		# Clustering and sctyping output
		File cluster_seurat_object = cluster_data.cluster_seurat_object
		File major_cell_type_plot = cluster_data.major_cell_type_plot
		File metadata = cluster_data.metadata

		# Group and feature plots for final metadata
		Array[File] group_umap_plots = plot_groups_and_features.group_umap_plots #!FileCoercion
		Array[File] feature_umap_plots = plot_groups_and_features.feature_umap_plots #!FileCoercion

		File cohort_analysis_manifest_tsv = upload_final_outputs.updated_manifest #!FileCoercion
	}
}

# Upload the list of samples used for this cohort analysis to the output bucket
task write_cohort_sample_list {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids

		String raw_data_path
		String billing_project
	}

	command <<<
		set -euo pipefail

		echo -e "project_id\tsample_id" > ~{cohort_id}.sample_list.tsv
		cat ~{write_tsv(project_sample_ids)} >> ~{cohort_id}.sample_list.tsv

		# Upload outputs
		gsutil -u ~{billing_project} -m cp \
			~{cohort_id}.sample_list.tsv \
			~{raw_data_path}/~{cohort_id}.sample_list.tsv
	>>>

	output {
		String cohort_sample_list = "~{raw_data_path}/~{cohort_id}.sample_list.tsv"
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim"
		cpu: 1
		memory: "1 GB"
		disks: "local-disk 10 HDD"
		preemptible: 3
	}
}

task filter_and_normalize {
	input {
		File preprocessed_seurat_object
		File unfiltered_metadata

		String raw_data_path
		String billing_project
		String container_registry

		# Purposefully unset
		String? my_none
	}

	Int threads = 2
	String seurat_object_basename = basename(preprocessed_seurat_object, "_01.rds")
	Int disk_size = ceil(size(preprocessed_seurat_object, "GB") * 4 + 20)

	command <<<
		set -euo pipefail

		# Filter cells
		/usr/bin/time \
		Rscript /opt/scripts/main/filter.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--seurat-object ~{preprocessed_seurat_object} \
			--metadata ~{unfiltered_metadata} \
			--output-seurat-object ~{seurat_object_basename}_filtered_02.rds

		# If any cells remain after filtering, the data is normalized and scaled
		if [[ -s "~{seurat_object_basename}_filtered_02.rds" ]]; then
			# Normalize and scale filtered cells
			/usr/bin/time \
			Rscript /opt/scripts/main/process.R \
				--working-dir "$(pwd)" \
				--script-dir /opt/scripts \
				--threads ~{threads} \
				--seurat-object ~{seurat_object_basename}_filtered_02.rds \
				--output-seurat-object ~{seurat_object_basename}_filtered_normalized_03.rds

			# Upload outputs
			gsutil -u ~{billing_project} -m cp \
				~{seurat_object_basename}_filtered_02.rds \
				~{seurat_object_basename}_filtered_normalized_03.rds \
				~{raw_data_path}/

			echo true > cells_remaining_post_filter.txt
		else
			echo false > cells_remaining_post_filter.txt
		fi
	>>>

	output {
		String? filtered_seurat_object = if read_boolean("cells_remaining_post_filter.txt") then "~{raw_data_path}/~{seurat_object_basename}_filtered_02.rds" else my_none
		String? normalized_seurat_object = if read_boolean("cells_remaining_post_filter.txt") then "~{raw_data_path}/~{seurat_object_basename}_filtered_normalized_03.rds" else my_none
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_6"
		cpu: threads
		memory: "8 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task plot_groups_and_features {
	input {
		String cohort_id
		File metadata

		Array[String] groups
		Array[String] features

		String raw_data_path
		String billing_project
		String container_registry
	}

	Int disk_size = ceil(size(metadata, "GB") * 4 + 20)

	command <<<
		set -euo pipefail

		# Plot groups
		declare -a group_plots
		while read -r group || [[ -n "${group}" ]]; do
			/usr/bin/time \
			Rscript /opt/scripts/main/plot_groups.R \
				--working-dir "$(pwd)" \
				--metadata ~{metadata} \
				--group "${group}" \
				--output-group-umap-plot "~{cohort_id}.${group}_group_umap.pdf"

			group_plots+=("~{cohort_id}.${group}_group_umap.pdf")
		done < ~{write_lines(groups)}

		# Upload outputs
		gsutil -u ~{billing_project} -m cp \
			"${group_plots[@]}" \
			~{raw_data_path}/

		echo "${group_plots[@]/#/~{raw_data_path}/}" \
		| tr ' ' '\n' \
		> group_plot_locs.txt

		# Plot features
		declare -a feature_plots
		while read -r feature || [[ -n "${feature}" ]]; do
			/usr/bin/time \
			Rscript /opt/scripts/main/plot_features.R \
				--working-dir "$(pwd)" \
				--metadata ~{metadata} \
				--feature "${feature}" \
				--output-feature-umap-plot "~{cohort_id}.${feature}_feature_umap.pdf"

			feature_plots+=("~{cohort_id}.${feature}_feature_umap.pdf")
		done < ~{write_lines(features)}

		# Upload outputs
		gsutil -u ~{billing_project} -m cp \
			"${feature_plots[@]}" \
			~{raw_data_path}/

		echo "${feature_plots[@]/#/~{raw_data_path}/}" \
		| tr ' ' '\n' \
		> feature_plot_locs.txt
	>>>

	output {
		Array[String] group_umap_plots = read_lines("group_plot_locs.txt")
		Array[String] feature_umap_plots = read_lines("feature_plot_locs.txt")
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_6"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}
