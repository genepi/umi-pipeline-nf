process STITCH_CONSENSUS {
    tag "${sample}"
    publishDir "${params.output}/${sample}/${target}/polishing/${type}", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path( smolecule_consensus )
        val ( type )
        path ( parsed_reference )

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${smolecule_consensus.baseName}_consensus.fastq" ), emit: consensus_fastq

    script:
        def medaka_subtool = workflow.manifest.version.matches( '<= 0.2.1') ? 'stitch' : 'sequence'
        def qualities = params.output_format == "fastq" ? "--qualities" : ""
    """
        medaka $medaka_subtool \
            $qualities \
            ${smolecule_consensus} \
            ${parsed_reference} \
            ${smolecule_consensus.baseName}_consensus.fastq
    """
}