version 1.0

task get_workflow_metadata {
	input {
		String zones
	}

	command <<<
		set -euo pipefail

		# UTC timestamp for the running workflow
		date -u +"%FT%H-%M-%SZ" > timestamp.txt

		# Billing project to use for file requests (matches the billing project used for compute)
		curl "http://metadata.google.internal/computeMetadata/v1/project/project-id" \
				-H "Metadata-Flavor: Google" \
		> billing_project.txt
	>>>

	output {
		String timestamp = read_string("timestamp.txt")
		String billing_project = read_string("billing_project.txt")
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 10 HDD"
		preemptible: 3
		zones: zones
	}
}
