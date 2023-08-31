version 1.0

# Call doublets and generate QC plots

workflow quality_control {
	input {
		String project_name
		Array[File] preprocessed_seurat_objects

		String container_registry
	}

	call doublets {
		input:
			project_name = project_name,
			preprocessed_seurat_objects = preprocessed_seurat_objects,
			container_registry = container_registry
	}

	call plot_qc {
		input:
			project_name = project_name,
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
		String project_name
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
			--project-name ~{project_name} \
			--output-metadata-file ~{project_name}.unfiltered_metadata.csv
	>>>

	output {
		File unfiltered_metadata = "~{project_name}.unfiltered_metadata.csv"
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
		String project_name
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
			--project-name ~{project_name} \
			--output-violin-plots ~{project_name}.qc.violin_plots.pdf \
			--output-umis-genes-plot ~{project_name}.qc.umis_genes_plot.pdf
	>>>

	output {
		File qc_violin_plots = "~{project_name}.qc.violin_plots.pdf"
		File qc_umis_genes_plot = "~{project_name}.qc.umis_genes_plot.pdf"
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
