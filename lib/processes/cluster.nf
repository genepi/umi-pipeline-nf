centroid_fasta="clusters_centroid.fasta"
consensus_fasta="clusters_consensus.fasta"
process CLUSTER {
    publishDir "${params.output}/cluster", mode: 'copy'

    input:
        path detected_umis_fasta
    output:
        path "${detected_umis_fasta.simpleName}/", emit: cluster

    """
        outdir_target=${detected_umis_fasta.simpleName}
        outdir_vsearch=\$outdir_target/vsearch_clusters
        
        mkdir -p \$outdir_vsearch && \
        vsearch --clusterout_id --clusters \$outdir_vsearch/cluster \
        --centroids \$outdir_target/${centroid_fasta} --consout \$outdir_target/${consensus_fasta} --minseqlength ${params.min_length} \
        --maxseqlength ${params.max_length} --threads ${params.threads} --cluster_fast ${detected_umis_fasta} --clusterout_sort \
        --gapopen 0E/5I --gapext 0E/2I --mismatch -8 --match 6 --iddef 0 --minwordmatches 0 --qmask none -id 0.85
    """
}