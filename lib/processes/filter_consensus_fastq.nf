filtered_fastq="masked_consensus.fastq"
process FILTER_CONSENSUS_FASTQ {
    publishDir "${params.output}/${sample}/fastq/${type}", mode: 'copy'
    
    input:
        tuple val( sample ), val( target ), path( merged_consensus_fastq )
        val ( type )
    
    when:
        params.output_format == "fastq"
    
    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${filtered_fastq}" ), emit: filtered_consensus_fastq
    
    script:
    """
        vsearch \
        --fastx_mask ${merged_consensus_fastq} \
        --fastq_qmax 90 \
        --fastq_qmin ${params.min_consensus_quality} \
        --fastqout ${filtered_fastq}
    """
    // fastq_qmax hardcoded to 90 (maximal possible Q-score value) to have no upper limit of quality
}