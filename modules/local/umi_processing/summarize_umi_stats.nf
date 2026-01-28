process SUMMARIZE_UMI_STATS {
    publishDir params.live_mode ? "${params.output}/${sample}/${target}/stats/${type}/${new Date().format('yyyyMMdd_HHmm')}/" : "${params.output}/${sample}/${target}/stats/${type}/", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(sample), val(target), path(stats_files)
    val type
    path python_summarize_umi_stats

    output:
    tuple val(sample), val(target), path("${sample}_${target}_umi_summary.tsv"), emit: umi_summary

    script:
    """
        python ${python_summarize_umi_stats} \
            ${stats_files.join(' ')} \
            --sample ${sample} \
            --target ${target} \
            -o ${sample}_${target}_umi_summary.tsv
    """
}
