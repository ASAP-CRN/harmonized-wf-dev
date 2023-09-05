version 1.0

# Run steps in the cohort analysis

import "run_quality_control/run_quality_control.wdl" as RunQualityControl
import "filter_data/filter_data.wdl" as FilterData
import "cluster_data/cluster_data.wdl" as ClusterData
import "plot_groups_and_features/plot_groups_and_features.wdl" as PlotGroupsAndFeatures

workflow cohort_analysis {
	input {
		String cohort_id
		Array[File] preprocessed_seurat_objects

		Int clustering_algorithm
		Float clustering_resolution
		File cell_type_markers_list

		Array[String] groups
		Array[String] features

		String run_timestamp
		String raw_data_path_prefix
		String curated_data_path_prefix
		String container_registry
	}

	# TODO dummy version
	# TODO this is the version of all underlying workflows, as well
	String workflow_name = "cohort_analysis"
	String workflow_version = "0.0.1"

	String raw_data_path = "~{raw_data_path_prefix}/~{workflow_name}/~{workflow_version}/~{run_timestamp}"
	String curated_data_path = "~{curated_data_path_prefix}/~{workflow_name}/~{workflow_version}/~{run_timestamp}"

	call RunQualityControl.run_quality_control {
		input:
			cohort_id = cohort_id,
			preprocessed_seurat_objects = preprocessed_seurat_objects,
			raw_data_path = raw_data_path,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	scatter (preprocessed_seurat_object in preprocessed_seurat_objects) {
		call FilterData.filter_data {
			input:
				preprocessed_seurat_object = preprocessed_seurat_object,
				unfiltered_metadata = run_quality_control.unfiltered_metadata,
				raw_data_path = raw_data_path,
				container_registry = container_registry
		}
	}

	call ClusterData.cluster_data {
		input:
			cohort_id = cohort_id,
			normalized_seurat_objects = filter_data.normalized_seurat_object,
			clustering_algorithm = clustering_algorithm,
			clustering_resolution = clustering_resolution,
			cell_type_markers_list = cell_type_markers_list,
			raw_data_path = raw_data_path,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	call PlotGroupsAndFeatures.plot_groups_and_features {
		input:
			cohort_id = cohort_id,
			metadata = cluster_data.metadata,
			groups = groups,
			features = features,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	output {
		# QC plots
		File qc_violin_plots = run_quality_control.qc_violin_plots
		File qc_umis_genes_plot = run_quality_control.qc_umis_genes_plot

		# Clustering and sctyping output
		File cluster_seurat_object = cluster_data.cluster_seurat_object
		File metadata = cluster_data.metadata

		# Group and feature plots for final metadata
		Array[File] group_umap_plots = plot_groups_and_features.group_umap_plots
		Array[File] feature_umap_plots = plot_groups_and_features.feature_umap_plots
	}
}
