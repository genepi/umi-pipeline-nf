include {MERGE_FASTQ} from '../modules/local/umi_processing/merge_input.nf'
include {SUBSAMPLING} from '../modules/local/umi_processing/subsampling.nf'
include {MAP_READS} from '../modules/local/map_reads.nf'
include {SPLIT_READS} from  '../modules/local/umi_processing/split_reads.nf'
include {DETECT_UMI_CONSENSUS_FASTQ as DETECT_UMI_FASTQ} from '../modules/local/umi_polishing/detect_umi_consensus_fastq.nf'
include {CLUSTER} from '../modules/local/umi_processing/cluster.nf'
include {REFORMAT_FILTER_CLUSTER} from '../modules/local/umi_processing/reformat_filter_cluster.nf'

workflow OFFLINE_UMI_PROCESSING {

    take:
        raw
        reference
        umi_filter_reads
        umi_extract
        umi_parse_clusters
        bed

    main:       
        Channel
            .fromPath("${params.input}/barcode*/*.fastq")
            .map{ 
                fastqs -> 
                def barcode = fastqs.parent.name
                tuple(barcode, fastqs)
                }
            .groupTuple( by: 0 ) 
            .set{ existing_fastqs }

        if( params.subsampling ){
            MERGE_FASTQ( existing_fastqs )
            SUBSAMPLING( MERGE_FASTQ.out.merged_fastq )
            SUBSAMPLING.out.subsampled_fastq
            .set { input_fastqs }
        } else {
            MERGE_FASTQ( existing_fastqs )
            .set { input_fastqs }
        }

        input_fastqs
            .splitFastq( by: params.chunk_size , file: true)
            .set{ chunked_input_fastqs }

        MAP_READS( chunked_input_fastqs, raw, reference )
        SPLIT_READS( MAP_READS.out.bam_consensus, bed, raw, umi_filter_reads )

        SPLIT_READS.out.split_reads_fastx
        .filter{ _sample, _target, fastq -> fastq.countFastq() > 0 }
        .set{ split_reads_filtered }

        DETECT_UMI_FASTQ( split_reads_filtered, raw, umi_extract )
        
        DETECT_UMI_FASTQ.out.umi_extract_fastq
        .groupTuple( by: [0, 1])
        .set{ extracted_umis }

        CLUSTER( extracted_umis, raw )

        // Filters the clusters to only keep cluser with more or equal than min_reads_per_cluster, but keeps the grouping per sample
        CLUSTER.out.cluster_fastas
            .map { barcode, target, clusters -> 
                def filtered_clusters = clusters.findAll { fasta -> fasta.countFasta() >= params.min_reads_per_cluster }
                filtered_clusters ? [barcode, target, filtered_clusters] : null
            }
            .filter { it != null }
            .set{ cluster_fastas }

        REFORMAT_FILTER_CLUSTER( cluster_fastas, raw, umi_parse_clusters )

        REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastqs
        .filter{ _sample, _type, fastqs, _task_index -> fastqs instanceof List}            
        .map{ sample, type, fastqs, _task_index ->
                tuple(sample, type, fastqs)
        }
        .set{ processed_umis }

        emit:
            processed_umis

}