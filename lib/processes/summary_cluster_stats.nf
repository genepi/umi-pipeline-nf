process SUMMARY_CLUSTER_STATS {    
    publishDir "${params.output}/cluster_stats/", mode: 'copy'
    
    input:
    path cluster_stats_file
    path umi_cluster_stats_summary

    output:
    path "summary_cluster_stats.tsv"
    path "summary_cluster_report.pdf"
    
    script:
    def current_summary = "${launchDir}/${params.output}/cluster_stats/summary_cluster_stats.tsv"
    """
    echo ${task.index}

    python ${umi_cluster_stats_summary} \
        --cluster-stat ${cluster_stats_file} \
        --current-summary $current_summary \
        --output-tsv summary_cluster_stats.tsv \
        --output-pdf summary_cluster_report.pdf 
    """
}