process MAP_CONSENSUS {

    publishDir "${params.output}/${sample}/align", mode: 'copy'

  input:
    tuple val( sample ), val( target ), path( consensus_fasta )
    path reference
  output:
    tuple path ( "*.bam" ), val( "${target}" ), path ( "*.bam.bai" ), emit: bam_consensus

  script:
  """
    minimap2 ${params.minimap2_param} -t ${params.threads} ${reference} ${consensus_fasta} | 
    samtools sort -@ 5 -o ${consensus_fasta.baseName}.bam - && samtools index -@ ${params.threads} ${consensus_fasta.baseName}.bam
  """
}