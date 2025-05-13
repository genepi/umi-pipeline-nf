process CREATE_CONSENSUS {
    input:
        tuple val( sample ), val( target ), path( smolecule_clusters_parsed_bam ), path( smolecule_clusters_parsed_bam_bai )
        val( type )

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${smolecule_clusters_parsed_bam.baseName}.hdf" ), emit: smolecule_consensus

    script:
        def medaka_subtool = workflow.manifest.version.matches( '<= 0.2.1') ? 'consensus' : 'inference'
    """
        medaka $medaka_subtool \
            --check_output \
            --save_features \
            --batch_size ${params.chunk_size} \
            --model ${params.medaka_model} \
            ${smolecule_clusters_parsed_bam} \
            ${smolecule_clusters_parsed_bam.baseName}.hdf
    """
}