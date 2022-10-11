process SPLIT_READS {
    input:
        path 1d_bam
        path bed
    output:

    """
        python umi_amplicon_tools/filter_reads_py ${params.min_overlap} -o . ${bed} ${1d_bam} 2>&1 | tee 
    """
}