process REFORMAT_CONSENSUS_CLUSTER {
    publishDir "${params.output}/${sample}/${target}/fastq/${type}", mode: 'copy', enabled: "${params.verbose}"

    input:
      tuple val( sample ), val( target ), path( cluster_consensus_fasta )
      val (type)
      path umi_reformat_consensus

    output:
      tuple val( "${sample}" ), val( "${target}" ), path( "*.fastq" ), emit:consensus_fastq

    script:
    """
        python ${umi_reformat_consensus} \
          --consensus_fasta ${cluster_consensus_fasta} \
          -o .
        
    """
}