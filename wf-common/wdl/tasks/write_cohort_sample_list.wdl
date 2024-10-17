version 1.0

# Upload the list of samples used for this cohort analysis to the output bucket
task write_cohort_sample_list {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	command <<<
		set -euo pipefail

		echo -e "project_id\tsample_id" > ~{cohort_id}.sample_list.tsv
		cat ~{write_tsv(project_sample_ids)} >> ~{cohort_id}.sample_list.tsv

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.sample_list.tsv"
	>>>

	output {
		String cohort_sample_list = "~{raw_data_path}/~{cohort_id}.sample_list.tsv"
	}

	runtime {
		docker: "~{container_registry}/util:1.1.0"
		cpu: 1
		memory: "1 GB"
		disks: "local-disk 10 HDD"
		preemptible: 3
		zones: zones
	}
}
