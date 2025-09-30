process PARSE_BED {
    publishDir "${params.output}/bed", mode: 'copy', enabled: "${params.verbose}"
    
    input:
    tuple val( target ), val( line )

    output:
    tuple val( target ), path("*.bed"), emit: bed_channel

    script:
    """
    echo '$line' > ${target}.bed
    """
}