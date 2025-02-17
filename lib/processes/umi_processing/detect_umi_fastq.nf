process DETECT_UMI_FASTQ {
    publishDir "${params.output}/${sample}/stats/${type}", pattern: "*.tsv", mode: 'copy'
    
    input:
        tuple val(sample), val(target), path(fastq)  // New fastq input
        path (cache_dir)
        val(type)
        path(umi_extract_python)

    output:
        tuple val("${sample}"), val("${target}"), path("${cache_dir}/${sample}/*fastq"), emit: umi_extract_fastq
        path("${cache_dir}/${sample}/*tsv")

    script:
        def output_filename = "extracted_umis_${sample}_${task.index}"  // New state directory
        def output_dir = "${cache_dir}/${sample}" 
        def write_report = params.write_reports ? "--tsv" : ""
        
    """
    mkdir -p $output_dir
    
    python ${umi_extract_python} \
        --fwd-umi ${params.fwd_umi} \
        --rev-umi ${params.rev_umi} \
        --max-error ${params.umi_errors} \
        --adapter_length ${params.adapter_length} \
        --output_format ${params.output_format} \
        --output_filename $output_filename \
        $write_report \
        -o $output_dir ${fastq}

    """
}

/*process DETECT_UMI_FASTQ {
    publishDir "${params.output}/${sample}/stats/${type}", pattern: "*.tsv", mode: 'copy'
    publishDir "${output_dir_clusters}", pattern: "*${params.output_format}", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path ( fastq )
        val ( type )
        path umi_extract_python
    
    output:
        tuple val( "${sample}" ), val( "${fastq.baseName}" ), path ( "*${params.output_format}" ), path("${output_dir_clusters}"), emit: umi_extract_fastq
        path "*.tsv"

    script:
        def write_report = params.write_reports ? "--tsv" : ""
        def output_dir_clusters="${params.output}/${sample}/${params.output_format}_umi/${type}"

    """
        python ${umi_extract_python} \
        --fwd-umi ${params.fwd_umi} \
        --rev-umi ${params.rev_umi} \
        --max-error ${params.umi_errors} \
        --adapter_length ${params.adapter_length} \
        --output_format ${params.output_format} \
        --output_filename ${fastq.baseName}_umis \
        $write_report \
        -o . ${fastq}
    """
}
*/