version 1.0

import "structs.wdl"
import "preprocess/preprocess.wdl" as Preprocess
import "quality_control/quality_control.wdl" as QualityControl
import "filter/filter.wdl" as Filter
import "cluster/cluster.wdl" as Cluster
import "plot/plot.wdl" as Plot

workflow harmonized_pmdbs_analysis {
	input {
		String project_name
		Array[Sample] samples

		File cellranger_reference_data

		Float soup_rate = 0.20

		Int clustering_algorithm = 3
		Float clustering_resolution = 0.3
		File cell_type_markers_list

		Array[String] groups = ["sample", "batch", "seurat_clusters"]
		Array[String] features = ["doublet_scores", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"]

		String container_registry
	}

	scatter (sample in samples) {
		call Preprocess.preprocess {
			input:
				sample = sample,
				cellranger_reference_data = cellranger_reference_data,
				soup_rate = soup_rate,
				container_registry = container_registry
		}
	}

	call QualityControl.quality_control {
		input:
			project_name = project_name,
			preprocessed_seurat_objects = preprocess.preprocessed_seurat_object,
			container_registry = container_registry
	}

	call Filter.filter {
		input:
			preprocessed_seurat_objects = preprocess.preprocessed_seurat_object,
			unfiltered_metadata = quality_control.unfiltered_metadata,
			container_registry = container_registry
	}

	call Cluster.cluster {
		input:
			project_name = project_name,
			normalized_seurat_objects = filter.normalized_seurat_objects,
			clustering_algorithm = clustering_algorithm,
			clustering_resolution = clustering_resolution,
			cell_type_markers_list = cell_type_markers_list,
			container_registry = container_registry
	}

	call Plot.plot {
		input:
			project_name = project_name,
			metadata = cluster.metadata,
			groups = groups,
			features = features,
			container_registry = container_registry
	}

	output {
		# Cellranger
		Array[File] raw_counts = preprocess.raw_counts
		Array[File] filtered_counts = preprocess.filtered_counts
		Array[File] molecule_info = preprocess.molecule_info
		Array[File] cellranger_metrics_csv = preprocess.cellranger_metrics_csv

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

	meta {
		description: "Harmonized postmortem-derived brain sequencing (PMDBS) workflow"
	}

	parameter_meta {
		project_name: {help: "Name of project"}
		samples: {help: "The set of samples and their associated reads and metadata"}
		cellranger_reference_data: {help: "Cellranger transcriptome reference data; see https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest."}
		soup_rate: {help: "Dataset contamination rate fraction; used to remove mRNA contamination from the RNAseq data [0.2]"}
		clustering_algorithm: {help: "Clustering algorithm to use. [3]"}
		clustering_resolution: {help: "Clustering resolution to use during clustering. [0.3]"}
		cell_type_markers_list: {help: "Seurat object RDS file containing a list of major cell type markers; used to annotate clusters."}
		groups: {help: "Groups to produce umap plots for. ['sample', 'batch', 'seurat_clusters']"}
		features: {help: "Features to produce umap plots for. ['doublet_scores', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb']"}
		container_registry: {help: "Container registry where Docker images are hosted"}
	}
}
