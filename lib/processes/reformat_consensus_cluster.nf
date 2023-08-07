process REFORMAT_CONSENSUS_CLUSTER {
    publishDir "${params.output}/${sample}/fasta/${type}", mode: 'copy'

    input:
      tuple val( sample ), val( target ), path( cluster_consensus_fasta )
      val (type)
      path umi_reformat_consensus

    output:
      tuple val( "${sample}" ), val( "${target}" ), path( "${type}.fasta" ), emit:consensus_fasta

    script:
    """
        cat ${cluster_consensus_fasta} | python "${umi_reformat_consensus}" > ${type}.fasta
    """
}