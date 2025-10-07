process DETECT_UMI_FASTQ {
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.tsv", mode: 'copy', enabled: "${params.verbose}"
    publishDir "${params.output}/${sample}/${target}/${params.output_format}_umi/${type}", pattern: "*${params.output_format}", mode: 'copy', enabled: "${params.verbose}"

    input:
        tuple val(sample), val(target), path(fastq)
        path (cache_dir)
        val(type)
        path(umi_extract_python)

    output:
        tuple val("${sample}"), val("${target}"), path("*${params.output_format}"), emit: umi_extract_fastq
        tuple val( "${sample}" ), val( "${target}" ), path("*tsv"), emit: umi_extract_fastq_stats

    script:
        def output_filename = "extracted_umis_${sample}_${task.index}"
        def output_dir = "${cache_dir}/${sample}/${target}" 
        def write_report = params.write_reports ? "--tsv" : ""
        def use_context = params.use_context ? "--use_context" : ""
    """
    mkdir -p $output_dir
    
    python ${umi_extract_python} \
        --fwd_umi ${params.fwd_umi} \
        --rev_umi ${params.rev_umi} \
        --fwd_primer ${params.fwd_context} \
        --rev_primer ${params.rev_context} \
        --max_error ${params.umi_errors} \
        --adapter_length ${params.adapter_length} \
        --output_format ${params.output_format} \
        --output_filename $output_filename \
        $use_context \
        $write_report \
        -o ./ ${fastq}

    cp ${output_filename}.${params.output_format} $output_dir/
    cp ${output_filename}.tsv $output_dir/

    """
}