include { GLUE_CLUSTERS                        } from '../modules/local/umi_polishing/glue_clusters.nf'
include { POLISH_CLUSTER                       } from '../modules/local/umi_polishing/polish_cluster.nf'
include { ALIGN_CLUSTER                        } from '../modules/local/umi_polishing/align.nf'
include { PARSE_BAM                            } from '../modules/local/umi_polishing/parse_bam.nf'
include { CREATE_CONSENSUS                     } from '../modules/local/umi_polishing/create_consensus.nf'
include { STITCH_CONSENSUS                     } from '../modules/local/umi_polishing/stitch_consensus.nf'
include { MERGE_CONSENSUS_FASTQ                } from '../modules/local/umi_polishing/merge_consensus_fastq.nf'
include { FILTER_CONSENSUS_FASTQ               } from '../modules/local/umi_polishing/filter_consensus_fastq.nf'
include { REFORMAT_CONSENSUS_CLUSTER           } from '../modules/local/umi_polishing/reformat_consensus_cluster.nf'
include { MAP_CONSENSUS ; MAP_CONSENSUS as MAP_FINAL_CONSENSUS } from '../modules/local/umi_polishing/map_consensus.nf'
include { DETECT_UMI_CONSENSUS_FASTQ           } from '../modules/local/umi_polishing/detect_umi_consensus_fastq.nf'
include { CLUSTER_CONSENSUS                    } from '../modules/local/umi_polishing/cluster_consensus.nf'

workflow UMI_POLISHING {
    take:
    processed_umis
    n_parsed_cluster
    consensus
    reference
    umi_parse_bam

    main:
    GLUE_CLUSTERS(processed_umis, consensus)
        .map { sample, target, clusters -> tuple(sample, target, clusters instanceof List ? clusters : [clusters]) }
        .set { glued_clusters }

    glued_clusters.map { sample, _target, clusters -> n_parsed_cluster.put("${sample}", clusters.size) }

    glued_clusters
        .transpose(by: 2)
        .set { glued_clusters_transposed }

    if (params.reference_based_polishing) {
        ALIGN_CLUSTER(glued_clusters_transposed, consensus, reference)
        PARSE_BAM(ALIGN_CLUSTER.out.smolecule_clusters_bam, consensus, reference, umi_parse_bam)
        CREATE_CONSENSUS(PARSE_BAM.out.smolecule_clusters_bam_parsed, consensus)

        CREATE_CONSENSUS.out.smolecule_consensus
            .join(PARSE_BAM.out.parsed_reference, by: [0, 1, 2])
            .set { smolecule_consensus_parsed_reference }

        STITCH_CONSENSUS(smolecule_consensus_parsed_reference, consensus)
        STITCH_CONSENSUS.out.consensus_fastq.set { consensus_fastq }
    }
    else {
        POLISH_CLUSTER(glued_clusters_transposed, consensus)
        POLISH_CLUSTER.out.consensus_fastq.set { consensus_fastq }
    }

    consensus_fastq
        .map { sample, target, fastq -> tuple(groupKey(sample, n_parsed_cluster.get("${sample}")), target, fastq) }
        .groupTuple(by: [0, 1])
        .set { merge_consensus }


    if (params.output_format == "fastq") {
        MERGE_CONSENSUS_FASTQ(merge_consensus, consensus)
        FILTER_CONSENSUS_FASTQ(MERGE_CONSENSUS_FASTQ.out.merged_consensus_fastq, consensus)
        FILTER_CONSENSUS_FASTQ.out.filtered_consensus_fastq.set { consensus_fastq }
    }
    else {
        MERGE_CONSENSUS_FASTQ(merge_consensus, consensus).set { consensus_fastq }
    }

    MAP_CONSENSUS(consensus_fastq, consensus, reference)

    emit:
    consensus_bam   = MAP_CONSENSUS.out.bam_consensus
    consensus_fastq
}
