version 1.0

# Run steps in the cohort analysis

import "quality_control/quality_control.wdl" as QualityControl
import "filter/filter.wdl" as Filter
import "cluster/cluster.wdl" as Cluster
import "plot/plot.wdl" as Plot

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

	call QualityControl.quality_control {
		input:
			cohort_id = cohort_id,
			preprocessed_seurat_objects = preprocessed_seurat_objects,
			raw_data_path = raw_data_path,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	scatter (preprocessed_seurat_object in preprocessed_seurat_objects) {
		call Filter.filter {
			input:
				preprocessed_seurat_object = preprocessed_seurat_object,
				unfiltered_metadata = quality_control.unfiltered_metadata,
				raw_data_path = raw_data_path,
				container_registry = container_registry
		}
	}

	call Cluster.cluster {
		input:
			cohort_id = cohort_id,
			normalized_seurat_objects = filter.normalized_seurat_object,
			clustering_algorithm = clustering_algorithm,
			clustering_resolution = clustering_resolution,
			cell_type_markers_list = cell_type_markers_list,
			raw_data_path = raw_data_path,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	call Plot.plot {
		input:
			cohort_id = cohort_id,
			metadata = cluster.metadata,
			groups = groups,
			features = features,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	output {
		# QC plots
		File qc_violin_plots = quality_control.qc_violin_plots
		File qc_umis_genes_plot = quality_control.qc_umis_genes_plot

		# Clustering and sctyping output
		File cluster_seurat_object = cluster.cluster_seurat_object
		File metadata = cluster.metadata

		# Group and feature plots for final metadata
		Array[File] group_umap_plots = plot.group_umap_plots
		Array[File] feature_umap_plots = plot.feature_umap_plots
	}
}
