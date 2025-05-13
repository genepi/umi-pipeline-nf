process CREATE_CONSENSUS {
    tag "${sample}"
    publishDir "${params.output}/${sample}/${target}/polishing/${type}", mode: 'copy'
    
    input:
        tuple val( sample ), val( target ), path( smolecule_clusters_parsed_bam ), path( smolecule_clusters_parsed_bam_bai )
        val( type )

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${smolecule_clusters_parsed_bam.baseName}.hdf" ), emit: smolecule_consensus
        path( "${smolecule_clusters_parsed_bam.baseName}.log" )


    script:
        def medaka_subtool = workflow.manifest.version.matches( '<= 0.2.1') ? 'consensus' : 'inference'
        def precision = workflow.manifest.version.matches( '<= 0.2.1') ? '' : '--full_precision'
    """
        medaka $medaka_subtool \
            --check_output \
            --save_features \
            $precision \
            --batch_size ${params.chunk_size} \
            --model ${params.medaka_model} \
            ${smolecule_clusters_parsed_bam} \
            ${smolecule_clusters_parsed_bam.baseName}.hdf 2> "${smolecule_clusters_parsed_bam.baseName}.log"
    """
}