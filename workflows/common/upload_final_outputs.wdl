version 1.0

task upload_final_outputs {
	input {
		String manifest_path
		Array[String] output_file_paths

		String workflow_name
		String workflow_version
		String run_timestamp

		String curated_data_path
		String billing_project
		String container_registry
	}

	command <<<
		set -euo pipefail

		echo -e "filename\tworkflow\tworkflow_version\ttimestamp" > new_files.manifest.tsv
		while read -r output_file || [[ -n "${output_file}" ]]; do
			gsutil -u ~{billing_project} -m cp "${output_file}" ~{curated_data_path}

			echo -e "$(basename "${output_file}")\t~{workflow_name}\t~{workflow_version}\t~{run_timestamp}" >> new_files.manifest.tsv
		done < ~{write_lines(output_file_paths)}

		if [[ ~{workflow_name} == "preprocess" ]] && gsutil ls ~{manifest_path}; then
			# If a manifest already exists, merge the previous and updated manifests
			#   and replace the existing manifest with the updated one
			gsutil -u ~{billing_project} -m cp ~{manifest_path} previous_manifest.tsv

			merge_manifests.py \
				--previous-manifest previous_manifest.tsv \
				--new-files-manifest new_files.manifest.tsv \
				--updated-manifest updated_manifest.tsv

			gsutil -u ~{billing_project} -m cp updated_manifest.tsv ~{manifest_path}
		else
			# If a manifest does not exist or if the workflow is not preprocess, create one (containing info about the files just uploaded)
			gsutil -u ~{billing_project} -m cp new_files_manifest.tsv ~{manifest_path}
		fi
	>>>

	output {
		String updated_manifest = manifest_path
	}

	runtime {
		docker: "~{container_registry}/util:0.0.1"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
	}
}
