nextflow.enable.dsl = 2

include { LIVE_UMI_PROCESSING                          } from './live-umi-processing.nf'
include { OFFLINE_UMI_PROCESSING                       } from './offline-umi-processing.nf'
include { UMI_POLISHING                                } from './umi-polishing.nf'
include { CONSENSUS_POLISHING                          } from './consensus_polishing.nf'
include { VARIANT_CALLING ; VARIANT_CALLING as CONSENSUS_VARIANT_CALLING } from './variant_calling.nf'
include { PARSE_BED                                    } from '../modules/local/parse_bed.nf'

workflow UMI_PIPELINE {

    //validate input parameters
    WorkflowMain.validate(params)

    // file paths
    bed = file("${params.bed}", checkIfExists: true)
    reference = file("${params.reference}", checkIfExists: true)
    reference_fai = file("${params.reference_fai}", checkIfExists: true)

    // python scripts
    umi_filter_reads = file("${projectDir}/bin/filter_reads.py", checkIfExists: true)
    umi_extract = file("${projectDir}/bin/extract_umis.py", checkIfExists: true)
    umi_parse_clusters = file("${projectDir}/bin/parse_clusters.py", checkIfExists: true)
    umi_reformat_consensus = file("${projectDir}/bin/reformat_consensus.py", checkIfExists: true)
    umi_cluster_report = file("${projectDir}/bin/cluster_report.py", checkIfExists: true)
    umi_cluster_stats_summary = file("${projectDir}/bin/summary_cluster_report.py", checkIfExists: true)
    umi_parse_bam = file("${projectDir}/bin/parse_cluster_alignment.py", checkIfExists: true)
    umi_summarize_filter_reads = file("${projectDir}/bin/summarize_split_read_stats.py", checkIfExists: true)
    umi_summarize_umi_stats = file("${projectDir}/bin/summarize_umi_stats.py", checkIfExists: true)
    // subdirectory and file prefixes
    raw = "raw"
    consensus = "consensus"
    final_consensus = "final"
    n_parsed_cluster = [:]
    // cluster_summary_output_path = "${params.output}/cluster_stats/summary_cluster_stats.tsv"
    def extracted_fastq_cache_dir = new File(".nextflow/cache/${workflow.sessionId}/extracted_fastq_cache_dir")
    extracted_fastq_cache_dir.mkdir()
    extracted_fastq_cache_dir_nf = file(extracted_fastq_cache_dir)
    def cluster_summary_cache_dir = new File(".nextflow/cache/${workflow.sessionId}/cluster_summary_cache_dir")
    cluster_summary_cache_dir.mkdir()
    cluster_summary_cache_dir_nf = file(cluster_summary_cache_dir)
    bed_input_ch = Channel.fromPath(bed)

    // make sure output directory exists
    def outputDir = file(params.output)
    outputDir.mkdirs()

    PARSE_BED(bed_input_ch)

    PARSE_BED.out.bed_files
        .flatten()
        .map { bed_file ->
            def target = bed_file.baseName
            tuple(target, bed_file)
        }
        .set { bed_ch }

    if (params.live) {
        LIVE_UMI_PROCESSING(
            raw,
            reference,
            umi_filter_reads,
            extracted_fastq_cache_dir_nf,
            umi_extract,
            umi_parse_clusters,
            umi_cluster_report,
            umi_cluster_stats_summary,
            umi_summarize_filter_reads,
            umi_summarize_umi_stats,
            cluster_summary_cache_dir_nf,
            bed_ch,
        )

        LIVE_UMI_PROCESSING.out.processed_umis.set { processed_umis }
    }
    else {
        OFFLINE_UMI_PROCESSING(
            raw,
            reference,
            umi_filter_reads,
            umi_extract,
            umi_parse_clusters,
            umi_cluster_report,
            umi_cluster_stats_summary,
            umi_summarize_filter_reads,
            umi_summarize_umi_stats,
            cluster_summary_cache_dir_nf,
            bed_ch,
        )

        OFFLINE_UMI_PROCESSING.out.processed_umis.set { processed_umis }
    }

    UMI_POLISHING(
        processed_umis,
        n_parsed_cluster,
        consensus,
        reference,
        umi_parse_bam,
    )


    if (params.call_variants) {
        VARIANT_CALLING(
            UMI_POLISHING.out.consensus_bam,
            consensus,
            reference,
            reference_fai,
            bed_ch,
        )
    }

    if (!params.reference_based_polishing) {
        CONSENSUS_POLISHING(
            UMI_POLISHING.out.consensus_fastq,
            consensus,
            final_consensus,
            reference,
            umi_extract,
            umi_reformat_consensus,
            umi_summarize_umi_stats,
        )

        if (params.call_variants) {
            CONSENSUS_VARIANT_CALLING(
                CONSENSUS_POLISHING.out.final_consensus_bam,
                final_consensus,
                reference,
                reference_fai,
                bed_ch,
            )
        }
    }
}
