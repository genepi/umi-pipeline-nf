process COPY_BED {
    publishDir "${params.output}", mode: 'copy'
    
    input:
        file bed

    output:
        file bed

    shell:
    """
        mkdir -p ${params.output}
        cp ${bed} ${params.output}
    """

}