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

		Boolean regenerate_preprocessed_seurat_objects = false

		Boolean run_cross_team_cohort_analysis = false
		String cohort_raw_data_bucket
		Array[String] cohort_staging_data_buckets

		Int clustering_algorithm = 3
		Float clustering_resolution = 0.3
		File cell_type_markers_list

		Array[String] groups = ["sample", "batch", "seurat_clusters"]
		Array[String] features = ["doublet_scores", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"]

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	Int multiome_container_revision = 14

	String workflow_execution_path = "workflow_execution"

	call get_workflow_metadata {
		input:
			zones = zones
	}

	scatter (project in projects) {
		String project_raw_data_path_prefix = "~{project.raw_data_bucket}/~{workflow_execution_path}"

		call Preprocess.preprocess {
			input:
				project_id = project.project_id,
				samples = project.samples,
				cellranger_reference_data = cellranger_reference_data,
				soup_rate = soup_rate,
				regenerate_preprocessed_seurat_objects = regenerate_preprocessed_seurat_objects,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = project_raw_data_path_prefix,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				multiome_container_revision = multiome_container_revision,
				zones = zones
		}

		if (project.run_project_cohort_analysis) {
			call CohortAnalysis.cohort_analysis as project_cohort_analysis {
				input:
					cohort_id = project.project_id,
					project_sample_ids = preprocess.project_sample_ids,
					preprocessed_seurat_objects = preprocess.seurat_object, # !FileCoercion
					preprocessing_output_file_paths = preprocess.preprocessing_output_file_paths,
					group_by_vars = ["batch"],
					clustering_algorithm = clustering_algorithm,
					clustering_resolution = clustering_resolution,
					cell_type_markers_list = cell_type_markers_list,
					groups = groups,
					features = features,
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = project_raw_data_path_prefix,
					staging_data_buckets = project.staging_data_buckets,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
					multiome_container_revision = multiome_container_revision,
					zones = zones
			}
		}
	}

	if (run_cross_team_cohort_analysis) {
		String cohort_raw_data_path_prefix = "~{cohort_raw_data_bucket}/~{workflow_execution_path}"

		call CohortAnalysis.cohort_analysis as cross_team_cohort_analysis {
			input:
				cohort_id = cohort_id,
				project_sample_ids = flatten(preprocess.project_sample_ids),
				preprocessed_seurat_objects = flatten(preprocess.seurat_object), # !FileCoercion
				group_by_vars = ["batch_id"],
				clustering_algorithm = clustering_algorithm,
				clustering_resolution = clustering_resolution,
				cell_type_markers_list = cell_type_markers_list,
				groups = groups,
				features = features,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = cohort_raw_data_path_prefix,
				staging_data_buckets = cohort_staging_data_buckets,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				multiome_container_revision = multiome_container_revision,
				zones = zones
		}
	}

	output {
		# Sample-level outputs
		## Cellranger
		Array[Array[File]] raw_counts = preprocess.raw_counts
		Array[Array[File]] filtered_counts = preprocess.filtered_counts
		Array[Array[File]] molecule_info = preprocess.molecule_info
		Array[Array[File]] cellranger_metrics_csvs = preprocess.metrics_csv

		# Project cohort analysis outputs
		## List of samples included in the cohort
		Array[File?] project_cohort_sample_list = project_cohort_analysis.cohort_sample_list

		## QC plots
		Array[Array[File]?] project_qc_plots_pdf = project_cohort_analysis.qc_plots_pdf
		Array[Array[File]?] project_qc_plots_png = project_cohort_analysis.qc_plots_png

		## Clustering and sctyping output
		Array[File?] project_integrated_seurat_object = project_cohort_analysis.integrated_seurat_object
		Array[File?] project_neighbors_seurat_object = project_cohort_analysis.neighbors_seurat_object
		Array[File?] project_umap_seurat_object = project_cohort_analysis.umap_seurat_object
		Array[File?] project_cluster_seurat_object = project_cohort_analysis.cluster_seurat_object
		Array[File?] project_major_cell_type_plot_pdf = project_cohort_analysis.major_cell_type_plot_pdf
		Array[File?] project_major_cell_type_plot_png = project_cohort_analysis.major_cell_type_plot_png
		Array[File?] project_metadata = project_cohort_analysis.metadata

		## Group and feature plots for final metadata
		Array[Array[File]?] project_group_umap_plots_pdf = project_cohort_analysis.group_umap_plots_pdf
		Array[Array[File]?] project_group_umap_plots_png = project_cohort_analysis.group_umap_plots_png
		Array[Array[File]?] project_feature_umap_plots_pdf = project_cohort_analysis.feature_umap_plots_pdf
		Array[Array[File]?] project_feature_umap_plots_png = project_cohort_analysis.feature_umap_plots_png

		Array[Array[File]?] preprocess_manifests = project_cohort_analysis.preprocess_manifest_tsvs
		Array[Array[File]?] project_manifests = project_cohort_analysis.cohort_analysis_manifest_tsvs

		# Cross-team cohort analysis outputs
		## List of samples included in the cohort
		File? cohort_sample_list = cross_team_cohort_analysis.cohort_sample_list

		## QC plots
		Array[File]? cohort_qc_plots_pdf = cross_team_cohort_analysis.qc_plots_pdf
		Array[File]? cohort_qc_plots_png = cross_team_cohort_analysis.qc_plots_png

		## Clustering and sctyping output
		File? cohort_integrated_seurat_object = cross_team_cohort_analysis.integrated_seurat_object
		File? cohort_neighbors_seurat_object = cross_team_cohort_analysis.neighbors_seurat_object
		File? cohort_umap_seurat_object = cross_team_cohort_analysis.umap_seurat_object
		File? cohort_cluster_seurat_object = cross_team_cohort_analysis.cluster_seurat_object
		File? cohort_major_cell_type_plot_pdf = cross_team_cohort_analysis.major_cell_type_plot_pdf
		File? cohort_major_cell_type_plot_png = cross_team_cohort_analysis.major_cell_type_plot_png
		File? cohort_metadata = cross_team_cohort_analysis.metadata

		## Group and feature plots for final metadata
		Array[File]? cohort_group_umap_plots_pdf = cross_team_cohort_analysis.group_umap_plots_pdf
		Array[File]? cohort_group_umap_plots_png = cross_team_cohort_analysis.group_umap_plots_png
		Array[File]? cohort_feature_umap_plots_pdf = cross_team_cohort_analysis.feature_umap_plots_pdf
		Array[File]? cohort_feature_umap_plots_png = cross_team_cohort_analysis.feature_umap_plots_png

		Array[File]? cohort_manifests = cross_team_cohort_analysis.cohort_analysis_manifest_tsvs
	}

	meta {
		description: "Harmonized postmortem-derived brain sequencing (PMDBS) workflow"
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files during cross-team cohort analysis."}
		projects: {help: "The project ID, set of samples and their associated reads and metadata, output bucket locations, and whether or not to run project-level cohort analysis."}
		cellranger_reference_data: {help: "Cellranger transcriptome reference data; see https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest."}
		soup_rate: {help: "Dataset contamination rate fraction; used to remove mRNA contamination from the RNAseq data. [0.2]"}
		regenerate_preprocessed_seurat_objects: {help: "Regenerate the preprocessed Seurat objects, even if these files already exist. [false]"}
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only preprocessing steps (cellranger and generating the initial seurat object(s)) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team cohort intermediate files to."}
		cohort_staging_data_bucket: {help: "Bucket to stage cross-team cohort analysis outputs in."}
		clustering_algorithm: {help: "Clustering algorithm to use. [3]"}
		clustering_resolution: {help: "Clustering resolution to use during clustering. [0.3]"}
		cell_type_markers_list: {help: "RDS file containing a list of major cell type markers; used to annotate clusters."}
		groups: {help: "Groups to produce umap plots for. ['sample', 'batch', 'seurat_clusters']"}
		features: {help: "Features to produce umap plots for. ['doublet_scores', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.rb']"}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones to spin up compute in."}
	}
}

task get_workflow_metadata {
	input {
		String zones
	}

	command <<<
		set -euo pipefail

		# UTC timestamp for the running workflow
		date -u +"%FT%H-%M-%SZ" > timestamp.txt

		# Billing project to use for file requests (matches the billing project used for compute)
		curl "http://metadata.google.internal/computeMetadata/v1/project/project-id" \
				-H "Metadata-Flavor: Google" \
		> billing_project.txt
	>>>

	output {
		String timestamp = read_string("timestamp.txt")
		String billing_project = read_string("billing_project.txt")
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 10 HDD"
		preemptible: 3
		zones: zones
	}
}
