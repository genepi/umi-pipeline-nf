process STITCH_CONSENSUS {
    tag "${sample}"
    publishDir "${params.output}/${sample}/${target}/polishing/${type}", mode: 'copy', enabled: "${params.verbose}"

    input:
        tuple val( sample ), val( target ), val( smolecule_cluster_name ), path( smolecule_consensus ), path ( parsed_reference )
        val ( type )

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${smolecule_cluster_name}_consensus.fastq" ), emit: consensus_fastq

    script:
        def medaka_subtool = workflow.manifest.version.matches( '<= 0.2.1') ? 'stitch' : 'sequence'
        def qualities = params.output_format == "fastq" ? "--qualities" : ""
    """
        medaka $medaka_subtool \
            $qualities \
            ${smolecule_consensus} \
            ${parsed_reference} \
            ${smolecule_cluster_name}_consensus.fastq
    """}
