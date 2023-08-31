version 1.0

# Plot groups and features from final metadata

workflow plot {
	input {
		String project_name

		File metadata

		Array[String] groups
		Array[String] features

		String container_registry
	}

	call plot_groups {
		input:
			project_name = project_name,
			metadata = metadata,
			groups = groups,
			container_registry = container_registry
	}

	call plot_features {
		input:
			project_name = project_name,
			metadata = metadata,
			features = features,
			container_registry = container_registry
	}

	output {
		# Group and feature plots
		Array[File] group_umap_plots = plot_groups.group_umap_plots
		Array[File] feature_umap_plots = plot_features.feature_umap_plots
	}
}

task plot_groups {
	input {
		String project_name
		File metadata

		Array[String] groups

		String container_registry
	}

	Int disk_size = ceil(size(metadata, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		while read -r group || [[ -n "${group}" ]]; do
			Rscript /opt/scripts/main/plot_groups.R \
				--working-dir "$(pwd)" \
				--metadata ~{metadata} \
				--group "${group}" \
				--output-group-umap-plot "~{project_name}.${group}_group_umap.pdf"
			done < ~{write_lines(groups)}
	>>>

	output {
		Array[File] group_umap_plots = glob("~{project_name}.*_group_umap.pdf")
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}

task plot_features {
	input {
		String project_name
		File metadata

		Array[String] features

		String container_registry
	}

	Int disk_size = ceil(size(metadata, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		while read -r feature || [[ -n "${feature}" ]]; do
			Rscript /opt/scripts/main/plot_features.R \
				--working-dir "$(pwd)" \
				--metadata ~{metadata} \
				--feature "${feature}" \
				--output-feature-umap-plot "~{project_name}.${feature}_feature_umap.pdf"
		done < ~{write_lines(features)}
	>>>

	output {
		Array[File] feature_umap_plots = glob("~{project_name}.*_feature_umap.pdf")
	}

	runtime {
		docker: "~{container_registry}/multiome:4a7fd84"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 20
	}
}
