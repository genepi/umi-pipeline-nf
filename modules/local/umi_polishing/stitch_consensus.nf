process STITCH_CONSENSUS {
    tag "${sample}"
    publishDir "${params.output}/${sample}/${target}/polishing/${type}", mode: 'copy'

    input:
        tuple val( sample ), val( target ), val( smolecule_cluster_name ), path( smolecule_consensus ), path ( parsed_reference )
        val ( type )

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${smolecule_cluster_name}_consensus.fastq" ), emit: consensus_fastq

    script:
    """
        medaka stitch \
            $qualities \
            ${smolecule_consensus} \
            ${parsed_reference} \
            ${smolecule_cluster_name}_consensus.fastq
    """
}
