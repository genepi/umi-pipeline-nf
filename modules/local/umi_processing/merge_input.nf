process MERGE_FASTQ {
    input:
        tuple val ( sample ), path ( fastq_path )
    
    output:
        tuple val( "${sample}" ), path( "*fastq" ), emit: merged_fastq
    
    script:
    def fastq_name = params.subsampling ? "${sample}_merged.fastq" : "filtered_${fastq_path.baseName.replaceAll(/(\.fastq|\.fq)(\.gz)?$/, "")}.fastq"
    """
        catfishq \
            --min-length ${params.min_read_length} \
            --min-qscore ${params.min_qscore} \
            ${fastq_path} > ${fastq_name}
    """
}