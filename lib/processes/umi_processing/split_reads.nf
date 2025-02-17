process SPLIT_READS {
    publishDir "${params.output}/${sample}/stats/${type}/${bam.baseName}", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output}/${sample}/${params.output_format}_filtered/${type}/${bam.baseName}", mode: 'copy', pattern: "*${params.output_format}"

    input:
        tuple val( sample ), val( target ), path ( bam ) , path ( bam_bai )
        path bed
        val( type )
        path python_filter_reads
    
    output:
        path "*"
        tuple val ( "${sample}" ), val( "target" ), path ( "*filtered.${params.output_format}" ), emit: split_reads_fastx
    
    script:
        def include_secondary_reads = params.include_secondary_reads ? "--include_secondary_reads" : ""
        def write_report = params.write_reports ? "--tsv" : ""
    """
        python ${python_filter_reads} \
          --min_overlap ${params.min_overlap} \
          --output_format ${params.output_format} \
          --adapter_length ${params.adapter_length} \
          --output_filename ${bam.baseName} \
          $include_secondary_reads \
          $write_report \
          -o ./ ${bed} \
          ${bam}
    """
}