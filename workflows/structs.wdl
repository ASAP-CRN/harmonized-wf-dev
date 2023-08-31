version 1.0

struct Sample {
	String sample_id
	String batch

	File fastq_R1
	File fastq_R2
}

struct Project {
	String project_id
	Array[Sample] samples

	String raw_data_bucket
	String curated_data_bucket
}
