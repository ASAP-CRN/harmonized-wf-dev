version 1.0

struct Sample {
	String sample_id
	String? batch

	File fastq_R1
	File fastq_R2
	File? fastq_I1
	File? fastq_I2
}

struct Project {
	String project_id
	Array[Sample] samples

	Boolean run_project_cohort_analysis

	String raw_data_bucket
	String staging_data_bucket
}
