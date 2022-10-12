process REFORMAT_FILTER_CLUSTERS {
    
    input:
        path vsearch_clusters
        path consensus_fasta
        path parse_clusters
    """
        python ${parse_clusters} --smolecule_out {output.out_file} {params.balance_strands_param} \
        --min_reads_per_clusters {params.min_reads_per_cluster} --max_reads_per_clusters {params.max_reads_per_cluster} \
        --stats_out {output.stats} -o {output.out_dir} {input}
    """

}