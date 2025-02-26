def consensus_fasta="consensus.fasta"
def vsearch_dir="vsearch_clusters"

process CLUSTER {
    //publishDir "${params.output}/${sample}/${target}/clustering/${type}", pattern: "${consensus_fasta}", mode: 'copy'
    //publishDir "${params.output}/${sample}/${target}/clustering/${type}", pattern: "cluster*", mode: 'copy'
    
    input:
        tuple val( sample ), val( target ), path( detected_umis_fastqs )
        val ( type )
    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${consensus_fasta}" ), optional: true, emit:consensus_fasta
        tuple val( "${sample}" ), val( "${target}" ), path( "cluster*" ), optional: true, emit:cluster_fastas
        
    script:
        def id = "${type}" == "raw" ? params.vsearch_sequence_identity : 0.99
    """
        cat ${detected_umis_fastqs} |
        vsearch \
        --clusterout_id \
        --clusters cluster \
        --consout ${consensus_fasta} \
        --minseqlength ${params.min_length} \
        --maxseqlength ${params.max_length} \
        --threads ${params.threads} \
        --cluster_fast - \
        --clusterout_sort \
        --gapopen 0E/5I \
        --gapext 0E/2I \
        --mismatch -8 \
        --match 6 \
        --iddef 0 \
        --minwordmatches 0 \
        --qmask none \
        --id $id
    """ 
}
