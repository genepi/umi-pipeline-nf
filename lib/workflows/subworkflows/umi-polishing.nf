include {GLUE_CLUSTERS} from '../../processes/glue_clusters.nf'
include {POLISH_CLUSTER} from '../../processes/polish_cluster.nf'
include {MERGE_CONSENSUS_FASTQ} from '../../processes/merge_consensus_fastq.nf'
include {FILTER_CONSENSUS_FASTQ} from '../../processes/filter_consensus_fastq.nf'
include {REFORMAT_CONSENSUS_CLUSTER} from '../../processes/reformat_consensus_cluster.nf'
include {MAP_READS as MAP_CONSENSUS; MAP_READS as MAP_FINAL_CONSENSUS} from '../../processes/map_reads.nf'
include {SPLIT_READS} from  '../../processes/split_reads.nf'
include {DETECT_UMI_CONSENSUS_FASTQ} from '../../processes/detect_umi_consensus_fastq.nf'
include {CLUSTER as CLUSTER_CONSENSUS} from '../../processes/cluster.nf'


// ----------------------------------------------------------------------------
// Subworkflow 2: auto_finish
// This subworkflow receives the output of the live_feedback phase and then
// automatically glues clusters, polishes consensus, maps consensus, and calls
// variants.
// ----------------------------------------------------------------------------
workflow UMI_POLISHING {
    take:
        processed_umis
        n_parsed_cluster
        consensus
        final_consensus
        reference
        umi_extract
        umi_reformat_consensus

    main:
         GLUE_CLUSTERS(processed_umis)
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
        


    emit: 
        consensus_bam = MAP_CONSENSUS.out.bam_consensus
        final_consensus_bam = MAP_FINAL_CONSENSUS.out.bam_consensus

}