process ALIGN_CLUSTER {
    tag "${sample}"
    publishDir "${params.output}/${sample}/${target}/polishing/${type}", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path( smolecule_clusters_fastq )
        val ( type )
        path ( reference )

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${smolecule_clusters_fastq.baseName}.bam" ), path( "${smolecule_clusters_fastq.baseName}.bam.bai" ), emit: smolecule_clusters_bam
    
    script:
    """
        mini_align \
            -d map-ont \
            -m \
            -r ${reference}
            -i ${smolecule_clusters_fastq}
            -p ${smolecule_clusters_fastq.baseName}.bam
    """
}