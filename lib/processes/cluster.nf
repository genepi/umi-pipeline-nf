centroid_fasta="clusters_centroid.fasta"
consensus_fasta="clusters_consensus.fasta"
process CLUSTER {

    input:
        path detected_umis_fasta
    output:

    """
        outdir_target=clustering/${detected_umis_fasta.baseName}
        outdir_vsearch=\$outdir_target/vsearch_clusters
        
        mkdir -p \$outdir_vsearch && \
        vsearch --clusterout_id --clusters \$outdir_vsearch \
        --centroids \$outdir_target/${centroid_fasta} --consout \$outdir_target/${consensus_fasta} --minseqlength ${params.min_length} \
        --maxseqlength ${params.max_length} --threads ${params.threads} --cluster_fast ${detected_umis_fasta} --clusterout_sort \
        --gapopen 0E/5I --gapext 0E/2I --mismatch -8 --match 6 --iddef 0 --minwordmatches 0 --qmask none -id 0.85
    """
}