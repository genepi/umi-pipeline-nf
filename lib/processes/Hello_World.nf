process HELLO_WORLD {
    publishDir "${params.output}", mode: 'symlink'

    input:
        file myFile

    output:
        path 'text.txt', emit : test
    
    """
        cat ${myFile} > text.txt
    """
 }