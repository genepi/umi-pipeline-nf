process MERGE_FASTQ {
    input:
        tuple val ( sample ), path ( fastq_path )
    
    output:
        tuple val( "${sample}" ), val( "target" ), path( "*fastq" ), emit: merged_fastq
    
    script:
    def fastq_name = "filtered_$fastq_path.Name"
    """
        catfishq \
            --min-length ${params.min_read_length} \
            --min-qscore ${params.min_qscore} \
            ${fastq_path} > ${fastq_name}
    """
}