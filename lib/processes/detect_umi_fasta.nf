process DETECT_UMI_FASTA {

    input:
        path split_reads_fastq
        path umi_extract_python
    """
        python ${umi_extract_python} --fwd-context ${params.fwd_universal_primer} --rev-context ${params.rev_universal_primer} --fwd-umi ${params.fwd_umi} --rev-umi ${params.rev_umi} --max-error ${params.umi_errors} ${split_reads_fastq} -o . --tsv test.tsv
    """
}