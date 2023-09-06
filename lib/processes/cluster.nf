centroid_fasta="centroid.fasta"
consensus_fasta="consensus.fasta"
vsearch_dir="vsearch_clusters"

process CLUSTER {
    publishDir "${params.output}/${sample}/clustering/${type}", pattern: "${consensus_fasta}" mode: 'copy'

    input:
        tuple val( sample ), val( target ), path( detected_umis_fasta )
        val ( type )
    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${consensus_fasta}" ), emit:consensus_fasta
        tuple val( "${sample}" ), val( "${target}" ), path( "cluster*" ), emit:cluster_fastas
        
    """
        vsearch \
        --clusterout_id \
        --clusters cluster \
        --centroids ${centroid_fasta} \
        --consout ${consensus_fasta} \
        --minseqlength ${params.min_length} \
        --maxseqlength ${params.max_length} \
        --threads ${params.threads} \
        --cluster_fast ${detected_umis_fasta} \
        --clusterout_sort \
        --gapopen 0E/5I \
        --gapext 0E/2I \
        --mismatch -8 \
        --match 6 \
        --iddef 0 \
        --minwordmatches 0 \
        --qmask none \
        --id 0.85
    """
}