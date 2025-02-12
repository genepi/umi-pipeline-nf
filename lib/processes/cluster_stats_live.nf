process CLUSTER_STATS_LIVE {
    tag "$sample"
    publishDir "${params.output}/${sample}/stats/${type}", mode: 'copy'
    
    input:
    tuple val(sample), val(type), path(smolecule_cluster_stats)  // Input tuple from the pipeline
    path cluster_stats_python
    
    output:
    path "${sample}_cluster_stats.tsv"
    path "${sample}_cluster_report.pdf"
    
    script:
    def min_read_per_cluster_stats = params.min_reads_per_cluster - 2
    """
    python ${cluster_stats_python} \
        --input ${smolecule_cluster_stats} \
        --output-pdf ${sample}_cluster_report.pdf \
        --output-tsv ${sample}_cluster_stats.tsv \
        --min-reads ${params.min_reads_per_cluster} \
        --threshold $min_read_per_cluster_stats
    """
}