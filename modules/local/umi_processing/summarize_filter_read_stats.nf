process SUMMARIZE_FILTER_READ_STATS {
    publishDir params.live ? "${params.output}/${sample}/${target}/stats/${type}/${new Date().format('yyyyMMdd_HHmm')}/" : "${params.output}/${sample}/${target}/stats/${type}/", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(sample), val(target), path(stats_files)
    val type
    path python_summarize_filter_reads

    output:
    tuple val(sample), val(target), path("${sample}_${target}_filtered_reads.tsv"), emit: summary_stats

    script:
    """
        python ${python_summarize_filter_reads} \
        ${stats_files.join(' ')} \
        --sample ${sample} \
        --target ${target} \
        -o ${sample}_${target}_filtered_reads.tsv
    """
}
