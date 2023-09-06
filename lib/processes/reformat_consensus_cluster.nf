process REFORMAT_CONSENSUS_CLUSTER {
    publishDir "${params.output}/${sample}/fastq/${type}", mode: 'copy'

    input:
      tuple val( sample ), val( target ), path( cluster_consensus_fastq )
      val (type)
      path umi_reformat_consensus

    output:
      tuple val( "${sample}" ), val( "${target}" ), path( "${type}.fastq" ), emit:consensus_fastq

    script:
    """
        cat ${cluster_consensus_fastq} | python "${umi_reformat_consensus}" > ${type}.fastq
    """
}