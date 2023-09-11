version 1.0

# Identify doublets and generate QC plots

workflow run_quality_control {
	input {
		String cohort_id
		Array[File] preprocessed_seurat_objects

		String raw_data_path
		String curated_data_path
		String container_registry
	}

	call identify_doublets {
		input:
			cohort_id = cohort_id,
			preprocessed_seurat_objects = preprocessed_seurat_objects,
			raw_data_path = raw_data_path,
			container_registry = container_registry
	}

	call plot_qc_metrics {
		input:
			cohort_id = cohort_id,
			unfiltered_metadata = identify_doublets.unfiltered_metadata, #!FileCoercion
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	output {
		# QC plots
		File qc_violin_plots = plot_qc_metrics.qc_violin_plots #!FileCoercion
		File qc_umis_genes_plot = plot_qc_metrics.qc_umis_genes_plot #!FileCoercion

		File unfiltered_metadata = identify_doublets.unfiltered_metadata #!FileCoercion
	}
}

task identify_doublets {
	input {
		String cohort_id
		Array[File] preprocessed_seurat_objects

		String raw_data_path
		String container_registry
	}

	Int threads = 2
	Int disk_size = ceil(size(preprocessed_seurat_objects[0], "GB") * length(preprocessed_seurat_objects) * 2 + 30)
	Int mem_gb = ceil(threads * 4 + length(preprocessed_seurat_objects) * 0.2)

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
		gsutil -m cp \
			~{cohort_id}.unfiltered_metadata.csv \
			~{raw_data_path}/
	>>>

	output {
		String unfiltered_metadata = "~{raw_data_path}/~{cohort_id}.unfiltered_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_1"
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

		String curated_data_path
		String container_registry
	}

	Int threads = 2
	Int disk_size = ceil(size(unfiltered_metadata, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		/usr/bin/time \
		Rscript /opt/scripts/main/plot_qc_metrics.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--metadata ~{unfiltered_metadata} \
			--project-name ~{cohort_id} \
			--output-violin-plots ~{cohort_id}.qc.violin_plots.pdf \
			--output-umis-genes-plot ~{cohort_id}.qc.umis_genes_plot.pdf

		# Upload outputs
		gsutil -m cp \
			~{cohort_id}.qc.violin_plots.pdf \
			~{cohort_id}.qc.umis_genes_plot.pdf \
			~{curated_data_path}/
	>>>

	output {
		String qc_violin_plots = "~{curated_data_path}/~{cohort_id}.qc.violin_plots.pdf"
		String qc_umis_genes_plot = "~{curated_data_path}/~{cohort_id}.qc.umis_genes_plot.pdf"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84_1"
		cpu: threads
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}
