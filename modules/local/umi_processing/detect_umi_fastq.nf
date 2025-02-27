process DETECT_UMI_FASTQ {
    publishDir "${params.output}/${sample}/${target}/stats/${type}", pattern: "*.tsv", mode: 'copy'
    
    input:
        tuple val(sample), val(target), path(fastq)
        path (cache_dir)
        val(type)
        path(umi_extract_python)

    output:
        tuple val("${sample}"), val("${target}"), path("${cache_dir}/${sample}/${target}/*fastq"), emit: umi_extract_fastq
        path("${cache_dir}/${sample}/${target}/*tsv")

    script:
        def output_filename = "extracted_umis_${sample}_${task.index}"
        def output_dir = "${cache_dir}/${sample}/${target}" 
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