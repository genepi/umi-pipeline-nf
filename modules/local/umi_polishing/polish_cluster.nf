process POLISH_CLUSTER {
    tag "${sample}"
    publishDir "${params.output}/${sample}/${target}/polishing/${type}", mode: 'copy', enabled: "${params.verbose}"

    input:
    tuple val(sample), val(target), path(smolecule_clusters_fastq)
    val type

    output:
    tuple val("${sample}"), val("${target}"), path("${smolecule_clusters_fastq.baseName}_consensus.fastq"), emit: consensus_fastq
    path "${smolecule_clusters_fastq.baseName}_smolecule.log"

    script:
    """
        medaka smolecule \
        --threads ${task.cpus} \
        --batch_size ${params.clusters_per_polishing_file} \
        --length 50 \
        --depth 2 \
        --model ${params.medaka_model} \
        --method spoa . \
        --qualities \
        ${smolecule_clusters_fastq} 2> "${smolecule_clusters_fastq.baseName}_smolecule.log"

        mv consensus.fastq ${smolecule_clusters_fastq.baseName}_consensus.fastq
    """
}
