include { REFORMAT_CONSENSUS_CLUSTER           } from '../modules/local/umi_polishing/reformat_consensus_cluster.nf'
include { MAP_CONSENSUS as MAP_FINAL_CONSENSUS } from '../modules/local/umi_polishing/map_consensus.nf'
include { DETECT_UMI_CONSENSUS_FASTQ           } from '../modules/local/umi_polishing/detect_umi_consensus_fastq.nf'
include { CLUSTER_CONSENSUS                    } from '../modules/local/umi_polishing/cluster_consensus.nf'
include { SUMMARIZE_UMI_STATS                  } from '../modules/local/umi_processing/summarize_umi_stats.nf'

workflow CONSENSUS_POLISHING {
    take:
    consensus_fastq
    consensus
    final_consensus
    reference
    umi_extract
    umi_reformat_consensus
    umi_summarize_umi_stats

    main:
    DETECT_UMI_CONSENSUS_FASTQ(consensus_fastq, consensus, umi_extract)

    DETECT_UMI_CONSENSUS_FASTQ.out.umi_extract_fastq_stats
        .groupTuple(by: [0, 1])
        .set { grouped_umi_stats }

    SUMMARIZE_UMI_STATS(grouped_umi_stats, consensus, umi_summarize_umi_stats)

    CLUSTER_CONSENSUS(DETECT_UMI_CONSENSUS_FASTQ.out.umi_extract_fastq, consensus)
    REFORMAT_CONSENSUS_CLUSTER(CLUSTER_CONSENSUS.out.consensus_fasta, final_consensus, umi_reformat_consensus)
    MAP_FINAL_CONSENSUS(REFORMAT_CONSENSUS_CLUSTER.out.consensus_fastq, final_consensus, reference)

    emit:
    final_consensus_bam = MAP_FINAL_CONSENSUS.out.bam_consensus
}
