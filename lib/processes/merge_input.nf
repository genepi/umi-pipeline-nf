merged_fastq="merged.fastq"
process MERGE_FASTQ {
    input:
        path sample
    output:
        tuple val( "${sample.baseName}" ), val( "target" ), path( "${merged_fastq}"), emit: merged_fastq
    """
        catfishq --min-length ${params.min_read_length} --min-qscore ${params.min_qscore} ${sample} > ${merged_fastq}
    """
}