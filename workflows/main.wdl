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

		Int n_top_genes = 8000
		String scvi_latent_key = "X_scvi"
		String clustering_method = "umap"
		# TODO - double check these defaults once clustering_mde.py is complete
		Int clustering_algorithm = 3
		Float clustering_resolution = 0.3
		File cell_type_markers_list

		Array[String] groups = ["sample", "batch", "cell_type"]
		Array[String] features = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_rb", "doublet_score", "S_score", "G2M_score"]

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	# Task and subworkflow versions
	String cellranger_task_version = "v1.1.0"

	String workflow_execution_path = "workflow_execution"

	call get_workflow_metadata {
		input:
			zones = zones
	}

	scatter (project in projects) {
		String project_raw_data_path_prefix = "~{project.raw_data_bucket}/~{workflow_execution_path}"

		scatter (sample_object in project.samples) {
			String cellranger_count_output = "~{project_raw_data_path_prefix}/cellranger/~{cellranger_task_version}/~{sample_object.sample_id}.raw_feature_bc_matrix.h5"
		}

		# For each sample, outputs an array of true/false
		call check_output_files_exist {
		input:
			cellranger_count_output_files = cellranger_count_output,
			billing_project = get_workflow_metadata.billing_project,
			zones = zones
		}

		scatter (index in range(length(project.samples))) {
			Sample sample = project.samples[index]
			String cellranger_count_complete = check_output_files_exist.sample_cellranger_complete[index][0]

			String cellranger_raw_counts = "~{project_raw_data_path_prefix}/cellranger/~{cellranger_task_version}/~{sample.sample_id}.raw_feature_bc_matrix.h5"
			String cellranger_filtered_counts = "~{project_raw_data_path_prefix}/cellranger/~{cellranger_task_version}/~{sample.sample_id}.filtered_feature_bc_matrix.h5"
			String cellranger_molecule_info = "~{project_raw_data_path_prefix}/cellranger/~{cellranger_task_version}/~{sample.sample_id}.molecule_info.h5"
			String cellranger_metrics_csv = "~{project_raw_data_path_prefix}/cellranger/~{cellranger_task_version}/~{sample.sample_id}.metrics_summary.csv"

			if (cellranger_count_complete == "false") {
				call cellranger_count {
					input:
						sample_id = sample.sample_id,
						fastq_R1s = sample.fastq_R1s,
						fastq_R2s = sample.fastq_R2s,
						fastq_I1s = sample.fastq_I1s,
						fastq_I2s = sample.fastq_I2s,
						cellranger_reference_data = cellranger_reference_data,
						raw_data_path = "~{project_raw_data_path_prefix}/cellranger/~{cellranger_task_version}",
						workflow_info = [[get_workflow_metadata.timestamp, "cellranger", cellranger_task_version]],
						billing_project = get_workflow_metadata.billing_project,
						container_registry = container_registry,
						zones = zones
				}
			}

			File raw_counts_output = select_first([cellranger_count.raw_counts, cellranger_raw_counts]) #!FileCoercion
			File filtered_counts_output = select_first([cellranger_count.filtered_counts, cellranger_filtered_counts]) #!FileCoercion
			File molecule_info_output = select_first([cellranger_count.molecule_info, cellranger_molecule_info]) #!FileCoercion
			File metrics_csv_output = select_first([cellranger_count.metrics_csv, cellranger_metrics_csv]) #!FileCoercion
		}

		call Preprocess.preprocess {
			input:
				project_id = project.project_id,
				samples = project.samples,
				raw_counts = raw_counts_output,
				cellbender_fpr = cellbender_fpr,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = project_raw_data_path_prefix,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}

		Array[String] preprocessing_output_file_paths = flatten([
			raw_counts_output,
			filtered_counts_output,
			molecule_info_output,
			metrics_csv_output,
			preprocess.report_html,
			preprocess.remove_background_counts,
			preprocess.filtered_remove_background_counts,
			preprocess.cell_barcodes_csv,
			preprocess.graph_pdf,
			preprocess.log,
			preprocess.metrics_csv,
			preprocess.checkpoint_tar_gz,
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
					clustering_algorithm = clustering_algorithm,
					clustering_resolution = clustering_resolution,
					clustering_method = clustering_method,
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
				clustering_algorithm = clustering_algorithm,
				clustering_resolution = clustering_resolution,
				clustering_method = clustering_method,
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
		Array[Array[File]] raw_counts = raw_counts_output
		Array[Array[File]] filtered_counts = filtered_counts_output
		Array[Array[File]] molecule_info = molecule_info_output
		Array[Array[File]] metrics_csvs = metrics_csv_output

		# Preprocess
		Array[Array[File]] report_html = preprocess.report_html
		Array[Array[File]] remove_background_counts = preprocess.remove_background_counts
		Array[Array[File]] filtered_remove_background_counts = preprocess.filtered_remove_background_counts
		Array[Array[File]] cell_barcodes_csv = preprocess.cell_barcodes_csv
		Array[Array[File]] graph_pdf = preprocess.graph_pdf
		Array[Array[File]] log = preprocess.log
		Array[Array[File]] metrics_csv = preprocess.metrics_csv
		Array[Array[File]] checkpoint_tar_gz = preprocess.checkpoint_tar_gz
		Array[Array[File]] posterior_probability = preprocess.posterior_probability
		Array[Array[File]] adata_object = preprocess.adata_object

		# Project cohort analysis outputs
		## List of samples included in the cohort
		Array[File?] project_cohort_sample_list = project_cohort_analysis.cohort_sample_list

		## Merged adata objects and QC plots
		Array[File?] project_merged_adata_object = project_cohort_analysis.merged_adata_object
		Array[Array[File]?] project_qc_plots_png = project_cohort_analysis.qc_plots_png

		# Clustering outputs
		Array[File?] project_integrated_adata_object = project_cohort_analysis.integrated_adata_object
		Array[File?] project_scvi_model = project_cohort_analysis.scvi_model
		Array[File?] project_umap_cluster_adata_object = project_cohort_analysis.umap_cluster_adata_object
		Array[File?] project_mde_cluster_adata_object = project_cohort_analysis.mde_cluster_adata_object
		Array[File?] project_major_cell_type_plot_pdf = project_cohort_analysis.major_cell_type_plot_pdf
		Array[File?] project_major_cell_type_plot_png = project_cohort_analysis.major_cell_type_plot_png
		Array[File?] project_cellassign_model = project_cohort_analysis.cellassign_model
		Array[File?] project_cell_types_csv = project_cohort_analysis.cell_types_csv

		# Groups and features plots
		Array[File?] project_groups_umap_plot_png = project_cohort_analysis.groups_umap_plot_png
		Array[File?] project_features_umap_plot_png = project_cohort_analysis.features_umap_plot_png

		Array[Array[File]?] preprocess_manifests = project_cohort_analysis.preprocess_manifest_tsvs
		Array[Array[File]?] project_manifests = project_cohort_analysis.cohort_analysis_manifest_tsvs

		# Cross-team cohort analysis outputs
		## List of samples included in the cohort
		File? cohort_sample_list = cross_team_cohort_analysis.cohort_sample_list

		## QC plots
		File? cohort_merged_adata_object = cross_team_cohort_analysis.merged_adata_object
		Array[File]? cohort_qc_plots_png = cross_team_cohort_analysis.qc_plots_png

		# Clustering outputs
		File? cohort_integrated_adata_object = cross_team_cohort_analysis.integrated_adata_object
		File? cohort_scvi_model = cross_team_cohort_analysis.scvi_model
		File? cohort_umap_cluster_adata_object = cross_team_cohort_analysis.umap_cluster_adata_object
		File? cohort_mde_cluster_adata_object = cross_team_cohort_analysis.mde_cluster_adata_object
		File? cohort_major_cell_type_plot_pdf = cross_team_cohort_analysis.major_cell_type_plot_pdf
		File? cohort_major_cell_type_plot_png = cross_team_cohort_analysis.major_cell_type_plot_png
		File? cohort_cellassign_model = cross_team_cohort_analysis.cellassign_model
		File? cohort_cell_types_csv = cross_team_cohort_analysis.cell_types_csv

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
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only preprocessing steps (cellranger and generating the initial seurat object(s)) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team cohort intermediate files to."}
		cohort_staging_data_buckets: {help: "Set of buckets to stage cross-team cohort analysis outputs in."}
		n_top_genes: {help: "Number of HVG genes to keep. [8000]"}
		scvi_latent_key: {help: "Latent key to save the scVI latent to. ['X_scvi']"}
		clustering_method: {help: "Clustering method; options are 'umap' or 'mde'. ['umap']"}
		clustering_algorithm: {help: "Clustering algorithm to use. [3]"}
		clustering_resolution: {help: "Clustering resolution to use during clustering. [0.3]"}
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

task check_output_files_exist {
	input {
		Array[String] cellranger_count_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do 
			if gsutil -u ~{billing_project} ls "${output_files}"; then
				echo "true" >> sample_cellranger_complete.tsv
			else
				echo "false" >> sample_cellranger_complete.tsv
			fi
		done < ~{write_lines(cellranger_count_output_files)}
	>>>

	output {
		Array[Array[String]] sample_cellranger_complete = read_tsv("sample_cellranger_complete.tsv")
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
		zones: zones
	}
}

task cellranger_count {
	input {
		String sample_id

		Array[File] fastq_R1s
		Array[File] fastq_R2s
		Array[File] fastq_I1s
		Array[File] fastq_I2s

		File cellranger_reference_data

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 16
	Int mem_gb = 24
	Int disk_size = ceil((size(fastq_R1s, "GB") + size(fastq_R2s, "GB") + size(fastq_I1s, "GB") + size(fastq_I2s, "GB") + size(cellranger_reference_data, "GB")) * 4 + 50)

	command <<<
		set -euo pipefail

		# Unpack refdata
		mkdir cellranger_refdata
		tar \
			-zxvf ~{cellranger_reference_data} \
			-C cellranger_refdata \
			--strip-components 1

		# Ensure fastqs are in the same directory
		mkdir fastqs
		while read -r fastq || [[ -n "${fastq}" ]]; do
			if [[ -n "${fastq}" ]]; then
				validated_fastq_name=$(fix_fastq_names --fastq "${fastq}" --sample-id "~{sample_id}")
				if [[ -e "fastqs/${validated_fastq_name}" ]]; then
					echo "[ERROR] Something's gone wrong with fastq renaming; trying to create fastq [${validated_fastq_name}] but it already exists. Exiting."
					exit 1
				else
					ln -s "${fastq}" "fastqs/${validated_fastq_name}"
				fi
			fi
		done < <(cat \
			~{write_lines(fastq_R1s)} \
			~{write_lines(fastq_R2s)} \
			~{write_lines(fastq_I1s)} \
			~{write_lines(fastq_I2s)})

		cellranger --version

		/usr/bin/time \
		cellranger count \
			--id=~{sample_id} \
			--transcriptome="$(pwd)/cellranger_refdata" \
			--fastqs="$(pwd)/fastqs" \
			--localcores ~{threads} \
			--localmem ~{mem_gb - 4}

		# Rename outputs to include sample ID
		mv ~{sample_id}/outs/raw_feature_bc_matrix.h5 ~{sample_id}.raw_feature_bc_matrix.h5
		mv ~{sample_id}/outs/filtered_feature_bc_matrix.h5 ~{sample_id}.filtered_feature_bc_matrix.h5
		mv ~{sample_id}/outs/molecule_info.h5 ~{sample_id}.molecule_info.h5
		mv ~{sample_id}/outs/metrics_summary.csv ~{sample_id}.metrics_summary.csv

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.raw_feature_bc_matrix.h5" \
			-o "~{sample_id}.filtered_feature_bc_matrix.h5" \
			-o "~{sample_id}.molecule_info.h5" \
			-o "~{sample_id}.metrics_summary.csv"
	>>>

	output {
		String raw_counts = "~{raw_data_path}/~{sample_id}.raw_feature_bc_matrix.h5"
		String filtered_counts = "~{raw_data_path}/~{sample_id}.filtered_feature_bc_matrix.h5"
		String molecule_info = "~{raw_data_path}/~{sample_id}.molecule_info.h5"
		String metrics_csv = "~{raw_data_path}/~{sample_id}.metrics_summary.csv"
	}

	runtime {
		docker: "~{container_registry}/cellranger:7.1.0"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}
