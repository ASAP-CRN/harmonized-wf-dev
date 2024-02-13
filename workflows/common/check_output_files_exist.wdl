version 1.0

# For each sample, outputs an array of true/false

task check_output_files_exist {
	input {
		Array[String] output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			if gsutil -u ~{billing_project} ls "${output_files}"; then
				echo "true" >> sample_complete.tsv
			else
				echo "false" >> sample_complete.tsv
			fi
		done < ~{write_lines(output_files)}
	>>>

	output {
		Array[Array[String]] sample_complete = read_tsv("sample_complete.tsv")
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
