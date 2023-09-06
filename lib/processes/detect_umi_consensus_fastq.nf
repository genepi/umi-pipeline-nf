process DETECT_UMI_CONSENSUS_FASTQ {
    publishDir "${params.output}/${sample}/stats/${type}", pattern: "*.tsv", mode: 'copy'
    publishDir "${params.output}/${sample}/fastq_umi/${type}", pattern: "*fastq", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( fastq )
        val ( type )
        path umi_extract_python
    
    output:
        tuple val( "${sample}" ), val( "${fastq.baseName}" ), path ( "*fastq" ), emit: umi_extract_fastq
        path "*.tsv"

    script:
        def write_report = "${params.write_reports}" ? "--tsv" : ""
    """
        python ${umi_extract_python} \
        --fwd-umi ${params.fwd_umi} \
        --rev-umi ${params.rev_umi} \
        --max-error ${params.umi_errors} \
        $write_report \
        -o . \
        ${fastq}
    """
}