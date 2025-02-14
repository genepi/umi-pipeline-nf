// ----------------------------------------------------------------------------
// Subworkflow 2: auto_finish
// This subworkflow receives the output of the live_feedback phase and then
// automatically glues clusters, polishes consensus, maps consensus, and calls
// variants.
// ----------------------------------------------------------------------------
workflow UMI_POLISHING {
    take:
        smolecule_cluster_fastqs_list

    main:
         GLUE_CLUSTERS(smolecule_cluster_fastqs_list)
            .map{ sample, type, clusters -> tuple(sample, type, clusters instanceof List ? clusters : [clusters]) }
            .set{ glued_clusters }
        glued_clusters
        .map{ sample, _type, clusters -> n_parsed_cluster.put("$sample", clusters.size)}
        
        glued_clusters
            .transpose(by: 2)
            .set { glued_clusters_transposed }

        POLISH_CLUSTER( glued_clusters_transposed, consensus )
        
        POLISH_CLUSTER.out.consensus_fastq
        .map{ sample, type, fastq -> tuple( groupKey(sample, n_parsed_cluster.get("$sample")), type, fastq) }
        .groupTuple( )
        .set{ merge_consensus }

        
        if ( params.output_format == "fastq"){
            MERGE_CONSENSUS_FASTQ(merge_consensus, consensus)
            FILTER_CONSENSUS_FASTQ(MERGE_CONSENSUS_FASTQ.out.merged_consensus_fastq, consensus)
            FILTER_CONSENSUS_FASTQ.out.filtered_consensus_fastq
            .set{ consensus_fastq }
        } else {
            MERGE_CONSENSUS_FASTQ(merge_consensus, consensus)
            .set{ consensus_fastq }
        }

        MAP_CONSENSUS( consensus_fastq, consensus, reference )
        DETECT_UMI_CONSENSUS_FASTQ( consensus_fastq, consensus, umi_extract )
        CLUSTER_CONSENSUS( DETECT_UMI_CONSENSUS_FASTQ.out.umi_extract_fastq , consensus )
        REFORMAT_CONSENSUS_CLUSTER( CLUSTER_CONSENSUS.out.consensus_fasta, final_consensus, umi_reformat_consensus )
        MAP_FINAL_CONSENSUS( REFORMAT_CONSENSUS_CLUSTER.out.consensus_fastq, final_consensus, reference )
        
        if( params.call_variants ){
            if( params.variant_caller == "lofreq" ){
                LOFREQ_CONSENSUS( MAP_CONSENSUS.out.bam_consensus, consensus, reference, reference_fai )
                LOFREQ_FINAL_CONSENSUS( MAP_FINAL_CONSENSUS.out.bam_consensus, final_consensus, reference, reference_fai )
            }else if( params.variant_caller == "mutserve"){
                MUTSERVE_CONSENSUS( MAP_CONSENSUS.out.bam_consensus, consensus, COPY_BED.out.bed, reference, reference_fai )
                MUTSERVE_FINAL_CONSENSUS( MAP_FINAL_CONSENSUS.out.bam_consensus, final_consensus, COPY_BED.out.bed, reference, reference_fai )
            }else if( params.variant_caller == "freebayes"){
                FREEBAYES_CONSENSUS( MAP_CONSENSUS.out.bam_consensus, consensus, reference, reference_fai )
                FREEBAYES_FINAL_CONSENSUS( MAP_FINAL_CONSENSUS.out.bam_consensus, final_consensus, reference, reference_fai )
            }else{
                exit 1, "${params.variant_caller} is not a valid option. \nPossible variant caller are <lofreq/mutserve/freebayes>"
            
            }
        }  

    emit: 
        MAP_CONSENSUS.out.bam_consensus
        MAP_FINAL_CONSENSUS.out.bam_consensus

}