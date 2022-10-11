process SPLIT_READS {
    input:
        path bam_1d
        path bed
    output:

    """
        python ${projectDir}/bin/filter_reads.py --min_overlap ${params.min_overlap} -o . ${bed} ${bam_1d} 2>&1 | tee 
    """
}