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

    main:
       println("We are live!")

        COPY_BED( bed )

        Channel
            .watchPath("${params.output}/CONTINUE")
            .take(1)
            .set{ continue_ch }

        CONTINUE_PIPELINE( continue_ch )
        

        Channel
        .fromPath("${params.input}/barcode*/*.fastq")
        .set{ existing_files_ch }

        Channel
        .watchPath("${params.input}/barcode*/*.fastq", 'create, modify')
        .until { it.getFileName().toString().toLowerCase().contains("continue") } 
        .set { watched_files_ch }

        existing_files_ch
        .concat( watched_files_ch )
        .map{ 
            fastq -> 
            def barcode = fastq.parent.name
            tuple(barcode, fastq)
            }
        .splitFastq( by: params.chunk_size , file: true)
        .set{ ch_input_files }
        
        MERGE_FASTQ( ch_input_files )
        .set { merged_fastq }

        MAP_READS( merged_fastq, raw, reference )
        SPLIT_READS( MAP_READS.out.bam_consensus, COPY_BED.out.bed, raw, umi_filter_reads )

        SPLIT_READS.out.split_reads_fastx
        .filter{ _sample, _target, fastq -> fastq.countFastq() > 0 }
        .set{ split_reads_filtered }

        // DETECT_UMI_FASTQ( split_reads_filtered, raw, umi_extract )
        DETECT_UMI_FASTQ( split_reads_filtered, extracted_fastq_cache_dir_nf, raw, umi_extract )
    

        CLUSTER_LIVE( DETECT_UMI_FASTQ.out.umi_extract_fastq, raw )

        CLUSTER_LIVE.out.cluster_fastas
            .map { barcode, target, clusters -> 
                def filtered_clusters = clusters.findAll { fasta -> fasta.countFasta() > params.min_reads_per_cluster }
                filtered_clusters ? [barcode, target, filtered_clusters] : null
            }
            .filter { it != null }
            .set{ cluster_fastas }
        
        REFORMAT_FILTER_CLUSTER( cluster_fastas, raw, umi_parse_clusters )

        // Launch the reporting process for each sample
        CLUSTER_STATS_LIVE( REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_stats, umi_cluster_report )
        SUMMARY_CLUSTER_STATS( REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_stats, umi_cluster_stats_summary)

        REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastqs
            .combine( continue_ch )
            .map { sample, type, fastqs, task_index, _continue_file ->
                tuple(sample, type, fastqs, [task_index])
            }
            .filter { _sample, _type, fastqs, _task_index -> fastqs instanceof List }
            .groupTuple(by: [0, 1], sort: { it[3] })
            .map { sample, type, fastqs, _task_index ->
                tuple(sample, type, fastqs[0])
            }
            .set { smolecule_cluster_fastqs_list }

    emit:
    // Emit the reformatted cluster FASTQs and the continue channel
    smolecule_cluster_fastqs_list

}