process SPLIT_READS {
    publishDir "${params.output}/${sample}/${target}/stats/${type}", mode: 'copy', pattern: "*.tsv", enabled: "${params.verbose}"
    publishDir "${params.output}/${sample}/${target}/${params.output_format}_filtered/${type}/${bam.baseName}", mode: 'copy', pattern: "*${params.output_format}", enabled: "${params.verbose}"

    input:
    tuple val(sample), path(bam), path(bam_bai), val(target), path(bed)
    val type
    path python_filter_reads

    output:
    path "*", optional: true
    tuple val("${sample}"), val("${target}"), path("*.tsv"), emit: split_reads_stats
    tuple val("${sample}"), val("${target}"), path("*filtered.${params.output_format}"), optional: true, emit: split_reads_fastx

    script:
    def include_secondary_reads = params.include_secondary_reads ? "--include_secondary_reads" : ""
    def write_report = params.write_reports ? "--tsv" : ""
    """
        python ${python_filter_reads} \
        --min_overlap ${params.min_overlap} \
        --output_format ${params.output_format} \
        --adapter_length ${params.adapter_length} \
        --output_filename ${sample}_${target}_${bam.baseName} \
        ${include_secondary_reads} \
        ${write_report} \
        -o ./ ${bed} \
        ${bam}
    """
}
