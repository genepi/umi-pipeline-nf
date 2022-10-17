merged_fastq="merged.fastq"
process MERGE_FASTQ {
    input:
        path sample
    output:
        tuple val( "${sample.baseName}" ), val( "target" ), path( "${merged_fastq}"), emit: merged_fastq
    """
        catfishq ${sample} > ${merged_fastq}
    """
}