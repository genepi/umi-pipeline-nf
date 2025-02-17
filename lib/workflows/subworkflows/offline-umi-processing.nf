include {MERGE_FASTQ} from '../../processes/merge_input.nf'
include {SUBSAMPLING} from '../../processes/subsampling.nf'
include {MAP_READS} from '../../processes/map_reads.nf'
include {SPLIT_READS} from  '../../processes/split_reads.nf'
include {DETECT_UMI_CONSENSUS_FASTQ as DETECT_UMI_FASTQ} from '../../processes/detect_umi_consensus_fastq.nf'
include {CLUSTER_LIVE as CLUSTER} from '../../processes/cluster_live.nf'
//include {CLUSTER} from '../../processes/cluster.nf'
include {REFORMAT_FILTER_CLUSTER} from '../../processes/reformat_filter_cluster.nf'

workflow OFFLINE_UMI_PROCESSING {

    take:
        existing_fastqs
        raw
        reference
        umi_filter_reads
        umi_extract
        umi_parse_clusters
        bed

    main:       
        existing_fastqs
        .map{ 
            fastq -> 
            def barcode = fastq.parent.name
            tuple(barcode, fastq)
            }
        .set{ existing_fastqs_annotated }

        if( params.subsampling ){
            MERGE_FASTQ( existing_fastqs_annotated )
            SUBSAMPLING( MERGE_FASTQ.out.merged_fastq )
            .set { merged_fastq }
        } else {
            MERGE_FASTQ( existing_fastqs_annotated )
            .set { merged_fastq }
        }

        merged_fastq
            .splitFastq( by: params.chunk_size , file: true)
            .set{ input_fastqs }

        MAP_READS( input_fastqs, raw, reference )
        SPLIT_READS( MAP_READS.out.bam_consensus, bed, raw, umi_filter_reads )

        SPLIT_READS.out.split_reads_fastx
        .filter{ _sample, _target, fastq -> fastq.countFastq() > 0 }
        .set{ split_reads_filtered }

        DETECT_UMI_FASTQ( split_reads_filtered, raw, umi_extract )
        
        DETECT_UMI_FASTQ.out.umi_extract_fastq
        .groupTuple( by: [0, 1])
        .set{ extracted_umis }

        CLUSTER( extracted_umis, raw )

        CLUSTER.out.cluster_fastas
            .map { barcode, target, clusters -> 
                def filtered_clusters = clusters.findAll { fasta -> fasta.countFasta() > params.min_reads_per_cluster }
                filtered_clusters ? [barcode, target, filtered_clusters] : null
            }
            .filter { it != null }
            .set{ cluster_fastas }

        REFORMAT_FILTER_CLUSTER( cluster_fastas, raw, umi_parse_clusters )

        REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastqs.view() 

        REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastqs
        .filter{ _sample, _type, fastqs, _task_index -> fastqs instanceof List}            
        .map{ sample, type, fastqs, _task_index ->
                tuple(sample, type, fastqs)
        }
        .set{ processed_umis }

        emit:
            processed_umis

}