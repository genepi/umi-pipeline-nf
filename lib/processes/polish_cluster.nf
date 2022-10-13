process POLISH_CLUSTER {

    publishDir "${params.output}/${sample}/fasta", pattern: "*.bam.bai", saveAs: { "${target}_consensus.bam.bai" }, mode: 'copy'
    publishDir "${params.output}/${sample}/fasta", pattern: "*.bam", saveAs: {"${target}_consensus.bam"}, mode: 'copy'
    publishDir "${params.output}/${sample}/fasta", pattern: "*.fasta", saveAs: {"${target}_consensus.fasta"}, mode: 'copy'
    
  input:
    tuple val( sample ), val( target ), path( smolecule_clusters_fasta )
  output:
    tuple val( "${sample}" ), val( "${target}" ), path( "*.bam" ), path( "*bam.bai")
    path "consensus.fasta"
    
  script:
  """
    medaka smolecule --threads ${params.threads} --length 50 --depth 2 --model ${params.medaka_model} --method spoa ${smolecule_clusters_fasta} . 2> ${target}_smolecule.log
  """
}