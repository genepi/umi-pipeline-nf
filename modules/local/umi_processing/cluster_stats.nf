process CLUSTER_STATS {
    tag "${sample}"
    publishDir params.live_mode ? "${params.output}/${sample}/${target}/stats/${type}/${new Date().format('yyyyMMdd_HHmm')}/" : "${params.output}/${sample}/${target}/stats/${type}/", mode: 'copy'

    input:
    tuple val(sample), val(target), path(smolecule_cluster_stats)
    val type
    path cluster_stats_python

    output:
    path "${sample}_cluster_stats.tsv", emit: cluster_stats
    path "${sample}_cluster_report.pdf"

    script:
    def min_reads_per_cluster_stats = params.min_reads_per_cluster - 2
    """
    python ${cluster_stats_python} \
        --input ${smolecule_cluster_stats} \
        --output-pdf ${sample}_cluster_report.pdf \
        --output-tsv ${sample}_cluster_stats.tsv \
        --min-reads ${params.min_reads_per_cluster} \
        --threshold ${min_reads_per_cluster_stats}
    """
}
