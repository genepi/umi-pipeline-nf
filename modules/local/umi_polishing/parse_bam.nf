process PARSE_BAM {
    input:
        tuple val( sample ), val( target ), val( smolecule_cluster_name ), path( smolecule_clusters_bam ), path( smolecule_clusters_bam_bai )
        val( type )
        path( reference )
        path umi_parse_bam

    output:
        tuple val( "${sample}" ), val( "${target}" ), val( "${smolecule_cluster_name}" ), path( "${smolecule_cluster_name}_parsed.bam" ), path( "${smolecule_cluster_name}_parsed.bam.bai" ), emit: smolecule_clusters_bam_parsed
        tuple val( "${sample}" ), val( "${target}" ), val( "${smolecule_cluster_name}" ), path ( "reference_parsed.fasta" ), emit: parsed_reference 

    script:
    """
        python3 ${umi_parse_bam} \
            --bam ${smolecule_clusters_bam} \
            --reference ${reference} \
            --output_reference reference_parsed.fasta \
            --output_bam ${smolecule_cluster_name}_parsed.bam
    """

}