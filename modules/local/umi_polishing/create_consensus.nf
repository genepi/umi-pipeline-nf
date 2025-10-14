process CREATE_CONSENSUS {
    tag "${sample}"
    publishDir "${params.output}/${sample}/${target}/polishing/${type}", mode: 'copy', enabled: "${params.verbose}"

    input:
    tuple val(sample), val(target), val(smolecule_cluster_name), path(smolecule_clusters_parsed_bam), path(smolecule_clusters_parsed_bam_bai)
    val type

    output:
    tuple val("${sample}"), val("${target}"), val("${smolecule_cluster_name}"), path("${smolecule_cluster_name}.hdf"), emit: smolecule_consensus
    path "${smolecule_cluster_name}.log"

    script:
    def medaka_subtool = workflow.manifest.version.matches('<= 0.2.1') ? 'consensus' : 'inference'
    def precision = workflow.manifest.version.matches('<= 0.2.1') ? '' : '--full_precision'
    """
        medaka ${medaka_subtool} \
            --check_output \
            --save_features \
            ${precision} \
            --batch_size ${params.chunk_size} \
            --model ${params.medaka_model} \
            ${smolecule_clusters_parsed_bam} \
            ${smolecule_cluster_name}.hdf 2> "${smolecule_cluster_name}.log"
    """
}
