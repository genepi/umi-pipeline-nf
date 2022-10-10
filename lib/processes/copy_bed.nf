process COPY_BED {
    publishDir "${params.output}", mode: 'copy'
    
    input:
        path bed

    output:
        path bed

    """
        mkdir -p ${params.output}
        cp ${bed} ${params.output}
    """

}