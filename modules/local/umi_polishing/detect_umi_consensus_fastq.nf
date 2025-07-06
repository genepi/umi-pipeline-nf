process DETECT_UMI_CONSENSUS_FASTQ {
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.tsv", mode: 'copy'
    publishDir "${params.output}/${sample}/${target}/${params.output_format}_umi/${type}", pattern: "*${params.output_format}", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( fastq )
        val ( type )
        path umi_extract_python
    
    output:
        tuple val( "${sample}" ), val( "${target}" ), path ( "*${params.output_format}" ), emit: umi_extract_fastq
        path "*.tsv"

    script:
        def write_report = params.write_reports ? "--tsv" : ""

    """
        python ${umi_extract_python} \
        --fwd_umi ${params.fwd_umi} \
        --rev_umi ${params.rev_umi} \
        --fwd_primer ${params.fwd_context} \
        --rev_primer ${params.rev_context} \
        --max_error ${params.umi_errors} \
        --adapter_length ${params.adapter_length} \
        --output_format ${params.output_format} \
        --output_filename ${fastq.baseName}_umis \
        $write_report \
        -o . ${fastq}
    """
}