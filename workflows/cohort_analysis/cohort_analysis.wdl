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

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] preprocessing_output_file_paths = []

		Array[String] group_by_vars

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		Array[String] groups
		Array[String] features

		String run_timestamp
		String raw_data_path_prefix
		String staging_data_path_prefix
		String billing_project
		String container_registry
		Int multiome_container_revision
		String zones
	}

	String workflow_name = "cohort_analysis"
	String workflow_version = "2.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version]]

	String raw_data_path = "~{raw_data_path_prefix}/~{workflow_name}/~{workflow_version}/~{run_timestamp}"

	Int n_samples = length(preprocessed_seurat_objects)

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

	call RunQualityControl.run_quality_control {
		input:
			cohort_id = cohort_id,
			preprocessed_seurat_objects = preprocessed_seurat_objects,
			n_samples = n_samples,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			multiome_container_revision = multiome_container_revision,
			zones = zones
	}

	scatter (preprocessed_seurat_object in preprocessed_seurat_objects) {
		# Filter sample data by nCount_RNA, nFeature_RNA, percent.mt, percent.rb; remove doublets
		# If cells remain after filtering, normalize and scale data
		call filter_and_normalize {
			input:
				preprocessed_seurat_object = preprocessed_seurat_object,
				unfiltered_metadata = run_quality_control.unfiltered_metadata,
				raw_data_path = raw_data_path,
				workflow_info = workflow_info,
				billing_project = billing_project,
				container_registry = container_registry,
				multiome_container_revision = multiome_container_revision,
				zones = zones
		}
	}

	call ClusterData.cluster_data {
		input:
			cohort_id = cohort_id,
			normalized_seurat_objects = select_all(filter_and_normalize.normalized_seurat_object), #!FileCoercion
			group_by_vars = group_by_vars,
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

	call plot_groups_and_features {
		input:
			cohort_id = cohort_id,
			metadata = cluster_data.metadata,
			groups = groups,
			features = features,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			multiome_container_revision = multiome_container_revision,
			zones = zones
	}

	Array[String] cohort_analysis_intermediate_output_paths = flatten([
		[
			run_quality_control.unfiltered_metadata
		],
		select_all(filter_and_normalize.filtered_seurat_object),
		select_all(filter_and_normalize.normalized_seurat_object),
		[
			cluster_data.integrated_seurat_object,
			cluster_data.neighbors_seurat_object,
			cluster_data.umap_seurat_object
		]
	]) #!StringCoercion

	call UploadFinalOutputs.upload_final_outputs as upload_preprocess_files {
		input:
			output_file_paths = flatten([preprocessing_output_file_paths, cohort_analysis_intermediate_output_paths]),
			staging_data_path = "~{staging_data_path_prefix}/preprocess",
			billing_project = billing_project,
			zones = zones
	}

	Array[String] cohort_analysis_final_output_paths = flatten([
		[
			write_cohort_sample_list.cohort_sample_list,
		],
		run_quality_control.qc_plots_png,
		[
			cluster_data.cluster_seurat_object,
			cluster_data.major_cell_type_plot_png,
			cluster_data.metadata
		],
		plot_groups_and_features.group_umap_plots_png,
		plot_groups_and_features.feature_umap_plots_png
	]) #!StringCoercion

	call UploadFinalOutputs.upload_final_outputs as upload_cohort_analysis_files {
		input:
			output_file_paths = cohort_analysis_final_output_paths,
			staging_data_path = "~{staging_data_path_prefix}/~{workflow_name}",
			billing_project = billing_project,
			zones = zones
	}

	output {
		File cohort_sample_list = write_cohort_sample_list.cohort_sample_list #!FileCoercion

		# QC plots
		Array[File] qc_plots_pdf = run_quality_control.qc_plots_pdf
		Array[File] qc_plots_png = run_quality_control.qc_plots_png

		# Clustering and sctyping output
		File integrated_seurat_object = cluster_data.integrated_seurat_object
		File neighbors_seurat_object = cluster_data.neighbors_seurat_object
		File umap_seurat_object = cluster_data.umap_seurat_object
		File cluster_seurat_object = cluster_data.cluster_seurat_object
		File major_cell_type_plot_pdf = cluster_data.major_cell_type_plot_pdf
		File major_cell_type_plot_png = cluster_data.major_cell_type_plot_png
		File metadata = cluster_data.metadata

		# Group and feature plots for final metadata
		Array[File] group_umap_plots_pdf = plot_groups_and_features.group_umap_plots_pdf #!FileCoercion
		Array[File] group_umap_plots_png = plot_groups_and_features.group_umap_plots_png #!FileCoercion
		Array[File] feature_umap_plots_pdf = plot_groups_and_features.feature_umap_plots_pdf #!FileCoercion
		Array[File] feature_umap_plots_png = plot_groups_and_features.feature_umap_plots_png #!FileCoercion

		File preprocess_manifest_tsv = upload_preprocess_files.manifest #!FileCoercion
		File cohort_analysis_manifest_tsv = upload_cohort_analysis_files.manifest #!FileCoercion
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

task filter_and_normalize {
	input {
		File preprocessed_seurat_object
		File unfiltered_metadata

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		Int multiome_container_revision
		String zones

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

			upload_outputs \
				-b ~{billing_project} \
				-d ~{raw_data_path} \
				-i ~{write_tsv(workflow_info)} \
				-o "~{seurat_object_basename}_filtered_02.rds" \
				-o "~{seurat_object_basename}_filtered_normalized_03.rds"

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
		docker: "~{container_registry}/multiome:4a7fd84_~{multiome_container_revision}"
		cpu: threads
		memory: "12 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
		zones: zones
	}
}

task plot_groups_and_features {
	input {
		String cohort_id
		File metadata

		Array[String] groups
		Array[String] features

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		Int multiome_container_revision
		String zones
	}

	Int disk_size = ceil(size(metadata, "GB") * 4 + 20)

	command <<<
		set -euo pipefail

		# Plot groups
		declare -a group_plots_pdf group_plots_png
		while read -r group || [[ -n "${group}" ]]; do
			/usr/bin/time \
			Rscript /opt/scripts/main/plot_groups.R \
				--working-dir "$(pwd)" \
				--metadata ~{metadata} \
				--group "${group}" \
				--output-group-umap-plot-prefix "~{cohort_id}.${group}_group_umap"

			group_plots_pdf+=("~{cohort_id}.${group}_group_umap.pdf")
			group_plots_png+=("~{cohort_id}.${group}_group_umap.png")
		done < ~{write_lines(groups)}

		# shellcheck disable=SC2068
		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			${group_plots_pdf[@]/#/-o } \
			${group_plots_png[@]/#/-o }

		echo "${group_plots_pdf[@]/#/~{raw_data_path}/}" \
		| tr ' ' '\n' \
		> group_plot_pdf_locs.txt

		echo "${group_plots_png[@]/#/~{raw_data_path}/}" \
		| tr ' ' '\n' \
		> group_plot_png_locs.txt

		# Plot features
		declare -a feature_plots_pdf feature_plots_png
		while read -r feature || [[ -n "${feature}" ]]; do
			/usr/bin/time \
			Rscript /opt/scripts/main/plot_features.R \
				--working-dir "$(pwd)" \
				--metadata ~{metadata} \
				--feature "${feature}" \
				--output-feature-umap-plot-prefix "~{cohort_id}.${feature}_feature_umap"

			feature_plots_pdf+=("~{cohort_id}.${feature}_feature_umap.pdf")
			feature_plots_png+=("~{cohort_id}.${feature}_feature_umap.png")
		done < ~{write_lines(features)}

		# shellcheck disable=SC2068
		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			${feature_plots_pdf[@]/#/-o } \
			${feature_plots_png[@]/#/-o }

		echo "${feature_plots_pdf[@]/#/~{raw_data_path}/}" \
		| tr ' ' '\n' \
		> feature_plot_pdf_locs.txt

		echo "${feature_plots_png[@]/#/~{raw_data_path}/}" \
		| tr ' ' '\n' \
		> feature_plot_png_locs.txt
	>>>

	output {
		Array[String] group_umap_plots_pdf = read_lines("group_plot_pdf_locs.txt")
		Array[String] group_umap_plots_png = read_lines("group_plot_png_locs.txt")
		Array[String] feature_umap_plots_pdf = read_lines("feature_plot_pdf_locs.txt")
		Array[String] feature_umap_plots_png = read_lines("feature_plot_png_locs.txt")
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_~{multiome_container_revision}"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
		zones: zones
	}
}
