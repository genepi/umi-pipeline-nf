fasta_filename = "_detected_umis.fasta"
process DETECT_UMI_FASTA {

    publishDir "${params.output}/fasta_umi", mode: 'copy'

    input:
        path split_reads_fastq
        path umi_extract_python
    output:
        path "*${fasta_filename}", emit: umi_extract_fasta
        path "*${fasta_filename}.tsv"
    """
        python ${umi_extract_python} --fwd-context ${params.fwd_universal_primer} --rev-context ${params.rev_universal_primer} --fwd-umi ${params.fwd_umi} --rev-umi ${params.rev_umi} --max-error ${params.umi_errors} ${split_reads_fastq} -o ${fasta_filename} --tsv ${fasta_filename}.tsv
    """
}