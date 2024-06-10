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

		# Preprocess
		Float cellbender_fpr = 0.0

		# Cohort analysis
		Boolean run_cross_team_cohort_analysis = false
		String cohort_raw_data_bucket
		Array[String] cohort_staging_data_buckets

		Int n_top_genes = 3000
		String scvi_latent_key = "X_scvi"
		String batch_key = "batch_id"
		File cell_type_markers_list

		Array[String] groups = ["sample", "batch", "cell_type"]
		Array[String] features = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_rb", "doublet_score", "S_score", "G2M_score"]

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

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
				cellbender_fpr = cellbender_fpr,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = project_raw_data_path_prefix,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}

		Array[String] preprocessing_output_file_paths = flatten([
			preprocess.raw_counts,
			preprocess.filtered_counts,
			preprocess.molecule_info,
			preprocess.metrics_summary_csv,
			preprocess.report_html,
			preprocess.removed_background_counts,
			preprocess.filtered_removed_background_counts,
			preprocess.cell_barcodes_csv,
			preprocess.graph_pdf,
			preprocess.log,
			preprocess.metrics_csv,
			preprocess.posterior_probability,
			preprocess.adata_object
		]) #!StringCoercion

		if (project.run_project_cohort_analysis) {
			call CohortAnalysis.cohort_analysis as project_cohort_analysis {
				input:
					cohort_id = project.project_id,
					project_sample_ids = preprocess.project_sample_ids,
					preprocessed_adata_objects = preprocess.adata_object,
					preprocessing_output_file_paths = preprocessing_output_file_paths,
					n_top_genes = n_top_genes,
					scvi_latent_key =scvi_latent_key,
					batch_key = batch_key,
					cell_type_markers_list = cell_type_markers_list,
					groups = groups,
					features = features,
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = project_raw_data_path_prefix,
					staging_data_buckets = project.staging_data_buckets,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
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
				preprocessed_adata_objects = flatten(preprocess.adata_object),
				preprocessing_output_file_paths = flatten(preprocessing_output_file_paths),
				n_top_genes = n_top_genes,
				scvi_latent_key =scvi_latent_key,
				batch_key = batch_key,
				cell_type_markers_list = cell_type_markers_list,
				groups = groups,
				features = features,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = cohort_raw_data_path_prefix,
				staging_data_buckets = cohort_staging_data_buckets,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}
	}

	output {
		# Sample-level outputs
		# Sample list
		Array[Array[Array[String]]] project_sample_ids = preprocess.project_sample_ids

		# Cellranger
		Array[Array[File]] cellranger_raw_counts = preprocess.raw_counts
		Array[Array[File]] cellranger_filtered_counts = preprocess.filtered_counts
		Array[Array[File]] cellranger_molecule_info = preprocess.molecule_info
		Array[Array[File]] cellranger_metrics_summary_csv = preprocess.metrics_summary_csv

		# Preprocess
		Array[Array[File]] cellbender_report_html = preprocess.report_html
		Array[Array[File]] cellbender_removed_background_counts = preprocess.removed_background_counts
		Array[Array[File]] cellbender_filtered_removed_background_counts = preprocess.filtered_removed_background_counts
		Array[Array[File]] cellbender_cell_barcodes_csv = preprocess.cell_barcodes_csv
		Array[Array[File]] cellbender_graph_pdf = preprocess.graph_pdf
		Array[Array[File]] cellbender_log = preprocess.log
		Array[Array[File]] cellbender_metrics_csv = preprocess.metrics_csv
		Array[Array[File]] cellbender_posterior_probability = preprocess.posterior_probability
		Array[Array[File]] adata_object = preprocess.adata_object

		# Project cohort analysis outputs
		## List of samples included in the cohort
		Array[File?] project_cohort_sample_list = project_cohort_analysis.cohort_sample_list

		# Merged adata objects, filtered and normalized adata objects, QC plots
		Array[File?] project_merged_adata_object = project_cohort_analysis.merged_adata_object
		Array[Array[File]?] project_qc_plots_png = project_cohort_analysis.qc_plots_png
		Array[File?] project_filtered_adata_object = project_cohort_analysis.filtered_adata_object
		Array[File?] project_normalized_adata_object = project_cohort_analysis.normalized_adata_object
		Array[File?] project_final_validation_metrics = project_cohort_analysis.final_validation_metrics

		# Clustering outputs
		Array[File?] project_integrated_adata_object = project_cohort_analysis.integrated_adata_object
		Array[File?] project_scvi_model_tar_gz = project_cohort_analysis.scvi_model_tar_gz
		Array[File?] project_umap_cluster_adata_object = project_cohort_analysis.umap_cluster_adata_object
		Array[File?] project_cell_annotated_adata_object = project_cohort_analysis.cell_annotated_adata_object
		Array[File?] project_cell_types_csv = project_cohort_analysis.cell_types_csv
		Array[File?] project_cell_annotated_metadata = project_cohort_analysis.cell_annotated_metadata

		# PCA and Harmony integrated adata objects and artifact metrics
		Array[File?] project_harmony_integrated_adata_object = project_cohort_analysis.harmony_integrated_adata_object
		Array[File?] project_scib_report_results_csv = project_cohort_analysis.scib_report_results_csv

		# Groups and features plots
		Array[File?] project_groups_umap_plot_png = project_cohort_analysis.groups_umap_plot_png
		Array[File?] project_features_umap_plot_png = project_cohort_analysis.features_umap_plot_png

		Array[Array[File]?] preprocess_manifests = project_cohort_analysis.preprocess_manifest_tsvs
		Array[Array[File]?] project_manifests = project_cohort_analysis.cohort_analysis_manifest_tsvs

		# Cross-team cohort analysis outputs
		## List of samples included in the cohort
		File? cohort_sample_list = cross_team_cohort_analysis.cohort_sample_list

		# Merged adata objects, filtered and normalized adata objects, QC plots
		File? cohort_merged_adata_object = cross_team_cohort_analysis.merged_adata_object
		Array[File]? cohort_qc_plots_png = cross_team_cohort_analysis.qc_plots_png
		File? cohort_filtered_adata_object = cross_team_cohort_analysis.filtered_adata_object
		File? cohort_normalized_adata_object = cross_team_cohort_analysis.normalized_adata_object
		File? cohort_final_validation_metrics = cross_team_cohort_analysis.final_validation_metrics

		# Clustering outputs
		File? cohort_integrated_adata_object = cross_team_cohort_analysis.integrated_adata_object
		File? cohort_scvi_model_tar_gz = cross_team_cohort_analysis.scvi_model_tar_gz
		File? cohort_umap_cluster_adata_object = cross_team_cohort_analysis.umap_cluster_adata_object
		File? cohort_cell_annotated_adata_object = cross_team_cohort_analysis.cell_annotated_adata_object
		File? cohort_cell_types_csv = cross_team_cohort_analysis.cell_types_csv
		File? cohort_cell_annotated_metadata = cross_team_cohort_analysis.cell_annotated_metadata

		# PCA and Harmony integrated adata objects and artifact metrics
		File? cohort_harmony_integrated_adata_object = cross_team_cohort_analysis.harmony_integrated_adata_object
		File? cohort_scib_report_results_csv = cross_team_cohort_analysis.scib_report_results_csv

		# Groups and features plots
		File? cohort_groups_umap_plot_png = cross_team_cohort_analysis.groups_umap_plot_png
		File? cohort_features_umap_plot_png = cross_team_cohort_analysis.features_umap_plot_png

		Array[File]? cohort_manifests = cross_team_cohort_analysis.cohort_analysis_manifest_tsvs
	}

	meta {
		description: "Harmonized postmortem-derived brain sequencing (PMDBS) workflow"
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files during cross-team cohort analysis."}
		projects: {help: "The project ID, set of samples and their associated reads and metadata, output bucket locations, and whether or not to run project-level cohort analysis."}
		cellranger_reference_data: {help: "Cellranger transcriptome reference data; see https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest."}
		cellbender_fpr :{help: "Cellbender false positive rate. [0.0]"}
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only preprocessing steps (cellranger and generating the initial adata object(s)) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team cohort intermediate files to."}
		cohort_staging_data_buckets: {help: "Set of buckets to stage cross-team cohort analysis outputs in."}
		n_top_genes: {help: "Number of HVG genes to keep. [8000]"}
		scvi_latent_key: {help: "Latent key to save the scVI latent to. ['X_scvi']"}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		cell_type_markers_list: {help: "CSV file containing a list of major cell type markers; used to annotate clusters."}
		groups: {help: "Groups to produce umap plots for. ['sample', 'batch', 'cell_type']"}
		features: {help: "Features to produce umap plots for. ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rb', 'doublet_score', 'S_score', 'G2M_score']"}
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
