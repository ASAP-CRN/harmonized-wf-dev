version 1.0

task upload_final_outputs {
	input {
		String manifest_path
		Array[String] output_file_paths

		String workflow_name
		String workflow_version
		String run_timestamp

		String staging_data_path
		String billing_project
		String container_registry
	}

	command <<<
		set -euo pipefail

		NEW_FILES_MANIFEST=new_files_manifest.tsv

		echo -e "filename\tworkflow\tworkflow_version\ttimestamp" > "${NEW_FILES_MANIFEST}"
		while read -r output_file || [[ -n "${output_file}" ]]; do
			gsutil -u ~{billing_project} -m cp "${output_file}" ~{staging_data_path}/

			echo -e "$(basename "${output_file}")\t~{workflow_name}\t~{workflow_version}\t~{run_timestamp}" >> "${NEW_FILES_MANIFEST}"
		done < ~{write_lines(output_file_paths)}

		if gsutil ls ~{manifest_path}; then
			# If a manifest already exists, merge the previous and updated manifests
			#   and replace the existing manifest with the updated one
			gsutil -u ~{billing_project} -m cp ~{manifest_path} previous_manifest.tsv

			update_manifest.py \
				--previous-manifest previous_manifest.tsv \
				--new-files-manifest "${NEW_FILES_MANIFEST}" \
				--updated-manifest updated_manifest.tsv

			gsutil -u ~{billing_project} -m cp updated_manifest.tsv ~{manifest_path}
		else
			# If a manifest does not exist, create one (containing info about the files just uploaded)
			gsutil -u ~{billing_project} -m cp "${NEW_FILES_MANIFEST}" ~{manifest_path}
		fi
	>>>

	output {
		String updated_manifest = manifest_path
	}

	runtime {
		docker: "~{container_registry}/util:1.0.0"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
	}
}

task sync_buckets {
	input {
		Array[String] source_buckets
		String target_bucket

		String billing_project
		String container_registry
	}

	command <<<
		set -euo pipefail

		while read -r source_bucket || [[ -n "${source_bucket}" ]]; do
			gsutil -u ~{billing_project} \
				-m rsync -r \
				"${source_bucket}/" \
				"~{target_bucket}/"
		done < ~{write_lines(source_buckets)}
	>>>

	output {
	}

	runtime {
		docker: "~{container_registry}/util:1.0.0"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
	}
}
