nextflow.enable.dsl = 2

include {COPY_BED} from '../../processes/copy_bed.nf'
include {CONTINUE_PIPELINE} from '../../processes/continue_pipeline.nf'
include {MERGE_FASTQ} from '../../processes/merge_input.nf'
include {MERGE_CONSENSUS_FASTQ} from '../../processes/merge_consensus_fastq.nf'
include {MAP_READS} from '../../processes/map_reads.nf'
include {SPLIT_READS} from  '../../processes/split_reads.nf'
include {DETECT_UMI_FASTQ} from '../../processes/detect_umi_fastq.nf'
include {CLUSTER as CLUSTER_CONSENSUS} from '../../processes/cluster.nf'
include {CLUSTER_LIVE} from '../../processes/cluster_live.nf'
include {REFORMAT_FILTER_CLUSTER} from '../../processes/reformat_filter_cluster.nf'
include {CLUSTER_STATS_LIVE} from '../../processes/cluster_stats_live.nf'
include {SUMMARY_CLUSTER_STATS} from '../../processes/summary_cluster_stats.nf'

// ----------------------------------------------------------------------------
// Subworkflow 1: live_feedback
// This subworkflow processes inputs, provides live user feedback via the
// clustering and summary_cluster_stats processes, and then emits outputs for
// further processing.
// ----------------------------------------------------------------------------
workflow LIVE_UMI_PROCESSING {
    /*
        Steps:
        1. Copy BED file & wait for a "continue" signal.
        2. Merge input FASTQs, map reads, split reads, and extract UMIs.
        3. Run live clustering and reformat/filter clusters.
        4. Launch cluster stats and summary processes to provide live feedback.
    */
    take:
        inputDir = params.input

    main:
        // Step 1: Copy BED and wait for continue signal
        COPY_BED(bed)
        Channel
            .watchPath("${params.output}/CONTINUE")
            .take(1)
            .set { continue_ch }
        continue_pipeline(continue_ch)

        // Step 2: Merge FASTQs (existing & newly added files)
        Channel.fromPath("${inputDir}/barcode*/*.fastq").set { existing_files_ch }
        Channel.watchPath("${inputDir}/barcode*/*.fastq", events: 'create,modify')
                .until { it.getFileName().toString().toLowerCase().contains("continue") }
                .set { watched_files_ch }
        existing_files_ch
            .concat( watched_files_ch )
            .map { fastq -> tuple( fastq.parent.name, fastq ) }
            .splitFastq( by: params.chunk_size, file: true )
            .set { ch_input_files }
        merge_fastq(ch_input_files).set { merged_fastq }

        // Step 3: Map reads and split them
        map_reads( merged_fastq, raw, reference )
            .set { map_out }
        split_reads( map_out.out.bam_consensus, bed, raw, umi_scripts.filter_reads )

        split_reads.out.split_reads_fastx
            .filter { sample, target, fastq -> fastq.countFastq() > 0 }
            .set { split_reads_filtered }

        detect_umi( split_reads_filtered, extractedCacheDir, raw, umi_scripts.extract_umis )
            .set { umi_out }

        // Step 4: Live clustering and reformatting
        cluster_live( umi_out.out.umi_extract_fastq, raw )
            .set { cluster_live_out }
        cluster_live_out.out.cluster_fastas
            .map { barcode, target, clusters ->
                def filtered = clusters.findAll { fasta -> fasta.countFasta() > params.min_reads_per_cluster }
                filtered ? tuple( barcode, target, filtered ) : null
            }
            .filter { it != null }
            .set { cluster_fastas }
        reformat_cluster( cluster_fastas, raw, umi_scripts.parse_clusters )
            .set { reformatted }

        // Step 5: Live feedback: report cluster stats & show summary
        cluster_stats( reformatted.out.smolecule_cluster_stats, umi_scripts.cluster_report )
        summary_stats( reformatted.out.smolecule_cluster_stats, umi_scripts.summary_cluster )

                // Recombine with continue signal and group cluster FASTQs
        reformattedOutput
            .combine( continueSignal )
            .map { sample, type, fastqs, task_index, _continue_file ->
                tuple( sample, type, fastqs, [task_index] )
            }
            .filter { _sample, _type, fastqs, _task_index -> fastq instanceof List || fastqs instanceof List }
            .groupTuple(by: [0,1], sort: { it[3] })
            .map { sample, type, fastqs, _task_index -> tuple( sample, type, fastqs[0] ) }
            .set { smolecule_cluster_fastqs_list }

    emit:
        // Emit the reformatted cluster FASTQs and the continue channel
        reformattedOutput = reformatted.out.smolecule_cluster_fastqs
        continueSignal    = continue_ch
}
