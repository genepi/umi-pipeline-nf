centroid_fasta="clusters_centroid.fasta"
consensus_fasta="clusters_consensus.fasta"
vsearch_dir="vsearch_clusters"

process CLUSTER {
    publishDir "${params.output}/${sample}/clustering/${type}", mode: 'copy'

    input:
        tuple val( sample ), val( target ), path( detected_umis_fasta )
        val ( type )
    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${vsearch_dir}/cluster*" ), emit:cluster_fastas
        path( "${centroid_fasta}" ) 
        
    """
        mkdir -p ${vsearch_dir} && \
        vsearch \
        --clusterout_id \
        --clusters ${vsearch_dir}/cluster \
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