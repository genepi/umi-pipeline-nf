process MERGE_FASTQ {
    input:
        tuple val ( sample ), path ( fastq_path )
    
    output:
        tuple val( "${sample}" ), path( "*fastq" ), emit: merged_fastq
    
    script:
    def fastq_name = params.live ? "filtered_${fastq_path.baseName.replaceAll(/(\.fastq|\.fq)(\.gz)?$/, "")}.fastq" : "${sample}_merged.fastq"
    """
        catfishq \
            --min-length ${params.min_read_length} \
            --min-qscore ${params.min_qscore} \
            ${fastq_path} > ${fastq_name}
    """
}