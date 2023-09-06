merged_fastq="merged_consensus.fastq"
process MERGE_CONSENSUS_FASTQ {
    
    input:
        tuple val( sample ), val( target ), path( consensus_fastqs )
    
    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${merged_fastq}"), emit: merged_consensus_fastq
    
    script:
    """
        cat ${consensus_fastqs} > ${merged_fastq}
    """
}