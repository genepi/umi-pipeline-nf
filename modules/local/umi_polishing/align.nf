process ALIGN_CLUSTER {
    tag "${sample}"
    publishDir "${params.output}/${sample}/${target}/polishing/${type}", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path( smolecule_clusters_fastq )
        val ( type )

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${smolecule_clusters_fastq.baseName}_consensus.fastq" ), emit: consensus_fastq
        path( "${smolecule_clusters_fastq.baseName}_smolecule.log" )
    
    script:
    """
    
    """
}