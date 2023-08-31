version 1.0

# Call doublets and generate QC plots

workflow quality_control {
	input {
		String cohort_id
		Array[File] preprocessed_seurat_objects

		String container_registry
	}

	call doublets {
		input:
			cohort_id = cohort_id,
			preprocessed_seurat_objects = preprocessed_seurat_objects,
			container_registry = container_registry
	}

	call plot_qc {
		input:
			cohort_id = cohort_id,
			unfiltered_metadata = doublets.unfiltered_metadata,
			container_registry = container_registry
	}

	output {
		# QC plots
		File qc_violin_plots = plot_qc.qc_violin_plots
		File qc_umis_genes_plot = plot_qc.qc_umis_genes_plot

		File unfiltered_metadata = doublets.unfiltered_metadata
	}
}

task doublets {
	input {
		String cohort_id
		Array[File] preprocessed_seurat_objects

		String container_registry
	}

	Int threads = 2
	Int disk_size = ceil(size(preprocessed_seurat_objects[0], "GB") * length(preprocessed_seurat_objects) * 2 + 30)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/gmm_doublet_calling.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--seurat-objects-fofn ~{write_lines(preprocessed_seurat_objects)} \
			--project-name ~{cohort_id} \
			--output-metadata-file ~{cohort_id}.unfiltered_metadata.csv
	>>>

	output {
		File unfiltered_metadata = "~{cohort_id}.unfiltered_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: threads
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task plot_qc {
	input {
		String cohort_id
		File unfiltered_metadata

		String container_registry
	}

	Int threads = 2
	Int disk_size = ceil(size(unfiltered_metadata, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/main/plot_qc_metrics.R \
			--working-dir "$(pwd)" \
			--script-dir /opt/scripts \
			--threads ~{threads} \
			--metadata ~{unfiltered_metadata} \
			--project-name ~{cohort_id} \
			--output-violin-plots ~{cohort_id}.qc.violin_plots.pdf \
			--output-umis-genes-plot ~{cohort_id}.qc.umis_genes_plot.pdf
	>>>

	output {
		File qc_violin_plots = "~{cohort_id}.qc.violin_plots.pdf"
		File qc_umis_genes_plot = "~{cohort_id}.qc.umis_genes_plot.pdf"
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: threads
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}
