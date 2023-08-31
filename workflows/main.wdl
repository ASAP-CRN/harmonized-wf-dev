version 1.0

# Harmonized workflow entrypoint

import "structs.wdl"
import "preprocess/preprocess.wdl" as Preprocess
import "cohort_analysis/cohort_analysis.wdl" as CohortAnalysis

workflow harmonized_pmdbs_analysis {
	input {
		String cohort_id
		Array[Project] projects

		File cellranger_reference_data

		Float soup_rate = 0.20

		Boolean run_cohort_analysis = false

		Int clustering_algorithm = 3
		Float clustering_resolution = 0.3
		File cell_type_markers_list

		Array[String] groups = ["sample", "batch", "seurat_clusters"]
		Array[String] features = ["doublet_scores", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"]

		String container_registry
	}

	scatter (project in projects) {
		scatter (sample in project.samples) {
			call Preprocess.preprocess {
				input:
					sample = sample,
					cellranger_reference_data = cellranger_reference_data,
					soup_rate = soup_rate,
					container_registry = container_registry
			}
		}

		if (project.run_cohort_analysis) {
			call CohortAnalysis.cohort_analysis as project_cohort_analysis {
				input:
					cohort_id = project.project_id,
					preprocessed_seurat_objects = preprocess.preprocessed_seurat_object,
					clustering_algorithm = clustering_algorithm,
					clustering_resolution = clustering_resolution,
					cell_type_markers_list = cell_type_markers_list,
					groups = groups,
					features = features,
					container_registry = container_registry
			}
		}
	}

	if (run_cohort_analysis) {
		call CohortAnalysis.cohort_analysis as cross_team_cohort_analysis {
			input:
				cohort_id = cohort_id,
				preprocessed_seurat_objects = flatten(preprocess.preprocessed_seurat_object),
				clustering_algorithm = clustering_algorithm,
				clustering_resolution = clustering_resolution,
				cell_type_markers_list = cell_type_markers_list,
				groups = groups,
				features = features,
				container_registry = container_registry
		}
	}

	output {
		# Sample-level outputs
		## Cellranger
		Array[Array[File]] raw_counts = preprocess.raw_counts
		Array[Array[File]] filtered_counts = preprocess.filtered_counts
		Array[Array[File]] molecule_info = preprocess.molecule_info
		Array[Array[File]] cellranger_metrics_csv = preprocess.cellranger_metrics_csv


		# Project cohort analysis outputs
		## QC plots
		Array[File?] project_qc_violin_plots = project_cohort_analysis.qc_violin_plots
		Array[File?] project_qc_umis_genes_plot = project_cohort_analysis.qc_umis_genes_plot

		## Clustering and sctyping output
		Array[File?] project_cluster_seurat_object = project_cohort_analysis.cluster_seurat_object
		Array[File?] project_metadata = project_cohort_analysis.metadata

		## Group and feature plots for final metadata
		Array[Array[File]?] project_group_umap_plots = project_cohort_analysis.group_umap_plots
		Array[Array[File]?] project_feature_umap_plots = project_cohort_analysis.feature_umap_plots


		# Cross-team cohort analysis outputs
		## QC plots
		File? cohort_qc_violin_plots = cross_team_cohort_analysis.qc_violin_plots
		File? cohort_qc_umis_genes_plot = cross_team_cohort_analysis.qc_umis_genes_plot

		## Clustering and sctyping output
		File? cohort_cluster_seurat_object = cross_team_cohort_analysis.cluster_seurat_object
		File? cohort_metadata = cross_team_cohort_analysis.metadata

		## Group and feature plots for final metadata
		Array[File]? cohort_group_umap_plots = cross_team_cohort_analysis.group_umap_plots
		Array[File]? cohort_feature_umap_plots = cross_team_cohort_analysis.feature_umap_plots
	}

	meta {
		description: "Harmonized postmortem-derived brain sequencing (PMDBS) workflow"
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files"}
		projects: {help: "The project ID, set of samples and their associated reads and metadata, as well as output bucket locations"}
		cellranger_reference_data: {help: "Cellranger transcriptome reference data; see https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest."}
		soup_rate: {help: "Dataset contamination rate fraction; used to remove mRNA contamination from the RNAseq data [0.2]"}
		run_cohort_analysis: {help: "Whether to run downstream harmonization steps. If set to false, only preprocessing steps (cellranger and generating the initial seurat object(s)) will run for samples. [false]"}
		clustering_algorithm: {help: "Clustering algorithm to use. [3]"}
		clustering_resolution: {help: "Clustering resolution to use during clustering. [0.3]"}
		cell_type_markers_list: {help: "Seurat object RDS file containing a list of major cell type markers; used to annotate clusters."}
		groups: {help: "Groups to produce umap plots for. ['sample', 'batch', 'seurat_clusters']"}
		features: {help: "Features to produce umap plots for. ['doublet_scores', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb']"}
		container_registry: {help: "Container registry where Docker images are hosted"}
	}
}
