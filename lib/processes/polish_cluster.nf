process POLISH_CLUSTER {

    publishDir "${params.output}/${sample}/fasta/${type}", pattern: "*.bam.bai", saveAs: { "${target}_consensus.bam.bai" }, mode: 'copy'
    publishDir "${params.output}/${sample}/fasta/${type}", pattern: "*.bam", saveAs: {"${target}_consensus.bam"}, mode: 'copy'
    publishDir "${params.output}/${sample}/fasta/${type}", pattern: "*.fasta", saveAs: {"${target}_consensus.fasta"}, mode: 'copy'
    
  input:
    tuple val( sample ), val( target ), path( smolecule_clusters_fasta )
    val ( type )
  output:
    tuple path( "*.bam" ), path( "*bam.bai")
    tuple val( "${sample}" ), val( "${target}" ), path( "consensus.fasta" ), emit: consensus_fasta

  script:
  """
    medaka smolecule --threads ${params.threads} --length 50 --depth 2 --model ${params.medaka_model} --method spoa . ${smolecule_clusters_fasta} 2> ${target}_smolecule.log
  """
}