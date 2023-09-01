version 1.0

# Plot groups and features from final metadata

workflow plot {
	input {
		String cohort_id

		File metadata

		Array[String] groups
		Array[String] features

		String curated_data_path
		String container_registry
	}

	call plot_groups {
		input:
			cohort_id = cohort_id,
			metadata = metadata,
			groups = groups,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	call plot_features {
		input:
			cohort_id = cohort_id,
			metadata = metadata,
			features = features,
			curated_data_path = curated_data_path,
			container_registry = container_registry
	}

	output {
		# Group and feature plots
		Array[File] group_umap_plots = plot_groups.group_umap_plots #!FileCoercion
		Array[File] feature_umap_plots = plot_features.feature_umap_plots #!FileCoercion
	}
}

task plot_groups {
	input {
		String cohort_id
		File metadata

		Array[String] groups

		String curated_data_path
		String container_registry
	}

	Int disk_size = ceil(size(metadata, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		declare -a group_plots
		while read -r group || [[ -n "${group}" ]]; do
			Rscript /opt/scripts/main/plot_groups.R \
				--working-dir "$(pwd)" \
				--metadata ~{metadata} \
				--group "${group}" \
				--output-group-umap-plot "~{cohort_id}.${group}_group_umap.pdf"

				group_plots+=("~{cohort_id}.${group}_group_umap.pdf")
			done < ~{write_lines(groups)}

		# Upload outputs
		gsutil -m cp \
			"${group_plots[@]}" \
			~{curated_data_path}/

		echo "${group_plots[@]/#/~{curated_data_path}/}" \
		| tr ' ' '\n' \
		> group_plot_locs.txt
	>>>

	output {
		Array[String] group_umap_plots = read_lines("group_plot_locs.txt")
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
		String cohort_id
		File metadata

		Array[String] features

		String curated_data_path
		String container_registry
	}

	Int disk_size = ceil(size(metadata, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		declare -a feature_plots
		while read -r feature || [[ -n "${feature}" ]]; do
			Rscript /opt/scripts/main/plot_features.R \
				--working-dir "$(pwd)" \
				--metadata ~{metadata} \
				--feature "${feature}" \
				--output-feature-umap-plot "~{cohort_id}.${feature}_feature_umap.pdf"

				feature_plots+=("~{cohort_id}.${feature}_feature_umap.pdf")
		done < ~{write_lines(features)}

		# Upload outputs
		gsutil -m cp \
			"${feature_plots[@]}" \
			~{curated_data_path}/

		echo "${feature_plots[@]/#/~{curated_data_path}/}" \
		| tr ' ' '\n' \
		> feature_plot_locs.txt
	>>>

	output {
		Array[String] feature_umap_plots = read_lines("feature_plot_locs.txt")
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
