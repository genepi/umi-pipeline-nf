process DETECT_UMI_CONSENSUS_FASTA {
    publishDir "${params.output}/${sample}/stats/${type}", pattern: "*.tsv", mode: 'copy'
    publishDir "${params.output}/${sample}/fasta_umi/${type}", pattern: "*fasta", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( fasta )
        val ( type )
        path umi_extract_python
    
    output:
        tuple val( "${sample}" ), val( "${fasta.baseName}" ), path ( "*fasta" ), emit: umi_extract_fasta
        path "*.tsv"

    script:
        def write_report = "${params.write_reports}" ? "--tsv" : ""
    """
        python ${umi_extract_python} \
        --fwd-umi ${params.fwd_umi} \
        --rev-umi ${params.rev_umi} \
        --max-error ${params.umi_errors} \
        $write_report
        -o . \
        ${fasta} \
    """
}