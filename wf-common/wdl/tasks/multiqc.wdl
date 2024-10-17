version 1.0

task multiqc {
	input {
		String project_id
		Array[File] output_files

		String output_name

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 2
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size(output_files, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		while read -r file || [[ -n "${file}" ]]; do
			if [[ "${file}" == *.tar.gz ]]; then
				tar -xzvf "${file}" --strip-components=1
			fi
		done < ~{write_lines(output_files)}

		multiqc . \
			--filename ~{project_id}.~{output_name} \
			--flat \
			--zip-data-dir

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{project_id}.~{output_name}.html" \
			-o "~{project_id}.~{output_name}_data.zip"
	>>>

	output {
		String multiqc_report_html =  "~{raw_data_path}/~{project_id}.~{output_name}.html"
		String multiqc_data_zip =  "~{raw_data_path}/~{project_id}.~{output_name}_data.zip"
	}

	runtime {
		docker: "~{container_registry}/multiqc:1.24.1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}
