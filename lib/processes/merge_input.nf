merged_fastq="merged.fastq"
process MERGE_FASTQ {
    
    input:
        tuple val ( sample ), path ( fastq_path )
    
    output:
        tuple val( "${sample}" ), val( "target" ), path( "${merged_fastq}"), emit: merged_fastq
    
    script:
    """
        catfishq \
            --min-length ${params.min_read_length} \
            --min-qscore ${params.min_qscore} \
            ${fastq_path} > ${merged_fastq}
    """
}