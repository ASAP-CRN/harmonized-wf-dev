version 1.0

# Identify doublets and generate QC plots

workflow run_quality_control {
	input {
		String cohort_id
		Array[File] preprocessed_seurat_objects
		Int n_samples

		String raw_data_path
		String billing_project
		String container_registry
		Int multiome_container_revision
	}

	call identify_doublets {
		input:
			cohort_id = cohort_id,
			preprocessed_seurat_objects = preprocessed_seurat_objects,
			n_samples = n_samples,
			raw_data_path = raw_data_path,
			billing_project = billing_project,
			container_registry = container_registry,
			multiome_container_revision = multiome_container_revision
	}

	call plot_qc_metrics {
		input:
			cohort_id = cohort_id,
			unfiltered_metadata = identify_doublets.unfiltered_metadata, #!FileCoercion
			n_samples = n_samples,
			raw_data_path = raw_data_path,
			billing_project = billing_project,
			container_registry = container_registry,
			multiome_container_revision = multiome_container_revision
	}

	output {
		# QC plots
		Array[File] qc_plots_pdf = plot_qc_metrics.qc_plots_pdf #!FileCoercion
		Array[File] qc_plots_png = plot_qc_metrics.qc_plots_png #!FileCoercion

		File unfiltered_metadata = identify_doublets.unfiltered_metadata #!FileCoercion
	}
}

task identify_doublets {
	input {
		String cohort_id
		Array[File] preprocessed_seurat_objects
		Int n_samples

		String raw_data_path
		String billing_project
		String container_registry
		Int multiome_container_revision
	}

	Int threads = 2
	Int disk_size = ceil(size(preprocessed_seurat_objects[0], "GB") * length(preprocessed_seurat_objects) * 2 + 30)
	Int mem_gb = ceil(0.2 * n_samples + threads * 4)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/gmm_doublet_calling.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-objects-fofn ~{write_lines(preprocessed_seurat_objects)} \
			--project-name ~{cohort_id} \
			--output-metadata-file ~{cohort_id}.unfiltered_metadata.csv

		# Upload outputs
		gsutil -u ~{billing_project} -m cp \
			~{cohort_id}.unfiltered_metadata.csv \
			~{raw_data_path}/
	>>>

	output {
		String unfiltered_metadata = "~{raw_data_path}/~{cohort_id}.unfiltered_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_~{multiome_container_revision}"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task plot_qc_metrics {
	input {
		String cohort_id
		File unfiltered_metadata
		Int n_samples

		String raw_data_path
		String billing_project
		String container_registry
		Int multiome_container_revision
	}

	Int threads = 2
	Int disk_size = ceil(size(unfiltered_metadata, "GB") * 2 + 20)
	Int mem_gb = ceil(0.02 * n_samples + threads * 2)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/plot_qc_metrics.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--metadata ~{unfiltered_metadata} \
			--project-name ~{cohort_id} \
			--output-violin-plot-prefix ~{cohort_id}.qc.violin_plots \
			--output-umis-genes-plot-prefix ~{cohort_id}.qc.umis_genes_plot

		# Upload outputs
		gsutil -u ~{billing_project} -m cp \
			~{cohort_id}.qc.violin_plots.pdf \
			~{cohort_id}.qc.violin_plots.png \
			~{cohort_id}.qc.umis_genes_plot.pdf \
			~{cohort_id}.qc.umis_genes_plot.png \
			~{raw_data_path}/
	>>>

	output {
		Array[String] qc_plots_pdf = [
			"~{raw_data_path}/~{cohort_id}.qc.violin_plots.pdf",
			"~{raw_data_path}/~{cohort_id}.qc.umis_genes_plot.pdf"
		]
		Array[String] qc_plots_png = [
			"~{raw_data_path}/~{cohort_id}.qc.violin_plots.png",
			"~{raw_data_path}/~{cohort_id}.qc.umis_genes_plot.pdf"
		]
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_~{multiome_container_revision}"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}
