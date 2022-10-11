stats_filename = "umi_filter_reads_stats.txt"
process SPLIT_READS {

    publishDir "${params.output}/stats", mode: 'copy'

    input:
        path bam_1d
        path bai_1d
        path bed
        path python_filter_reads
    output:
        path "${stats_filename}"

    """
        python ${python_filter_reads} --min_overlap ${params.min_overlap} -o . ${bed} ${bam_1d} 2>&1 | tee ${stats_filename}
    """
}