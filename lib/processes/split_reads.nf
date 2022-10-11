stats_filename = "umi_filter_reads_stats.txt"
process SPLIT_READS {

    publishDir "${params.output}/stats", mode: 'copy', pattern: "${stats_filename}"
    publishDir "${params.output}/fasta_filtered", mode: 'copy', pattern: '*.fastq'

    input:
        path bam_1d
        path bai_1d
        path bed
        path python_filter_reads
    output:
        path "${stats_filename}"
        path "*.fastq", emit: split_reads_fastq

    """
        python ${python_filter_reads} --min_overlap ${params.min_overlap} -o . ${bed} ${bam_1d} 2>&1 | tee ${stats_filename}
    """
}