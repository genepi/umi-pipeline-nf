#!/usr/bin/env nextflow

// ENABLE DSL2
nextflow.enable.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\


         Usage: 
              nextflow run AmstlerStephan/umi-pipeline-nf [OPTIONS]...

         Options: BASIC PARAMETERS
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    

         Options: GENERAL - Required	
              --input [path/to/input/dir]     [REQUIRED] Input directory containing (zipped) FASTQ files
              --output [STR]                  [REQUIRED] A string that can be given to name the output directory
              --reference [path/to/ref.fa]    [REQUIRED] Path to the reference genome in fasta format
              --reference_fai [path/to/fai]   [REQUIRED] Path to the reference index file
              --bed [path/to/data.bed]        [REQUIRED] Path to the bed file
              --threads                       Number of maximum threads to use [default: availableProcessors -1]
          
         Options: READ FILTERING
              --min_read_length               flag to enable subsampling [default: 0]
              --min_qscore                    Seed to produce pseudorandom numbers [default: 0]

         Options: SUBSAMPLING
              --subsampling                   flag to enable subsampling [default: false]
              --subsampling_seed              Seed to produce pseudorandom numbers [default: 11]
              --subsampling_readnumber        Number of reads after subsampling [default: 100000]

         Options: VARIANT CALLING
              --call_variants                 flag to enable variant calling [default: false]
              --variant_caller [STR]          [REQUIRED if call_variants is set] Variant caller [lofreq | mutserve | freebayes ]

           Options: ADVANCED
              --min_reads_per_barcode         Minimal number of fastq reads for each barcode [default: 1000]
              --umi_errors                    Max differences between extracted UMIs of the read and UMI pattern [default: 3]
              --min_reads_per_cluster         Min number of raw reads required for a consensus read [default: 20]
              --max_reads_per_cluster         Max number of raw reads used for a consensus read [default: 60]
              --filter_strategy_clusters      Filtering strategy for clusters with more than max_reads_per_cluster reads [random | quality] [default: quality]
              --output_format                 Output format until the cluster filtering step [fasta | fastq] [default: fastq]
              --write_reports                 Write stats of cluster and cluster filtering [default: true]
              --min_overlap                   Min overlap with target region [default: 0.90]
              --balance_strands               Balance forward and reverse raw reads in clusters [default: true]
              --medaka_model                  Medaka model used to compute consensus reads [default: "r1041_e82_400bps_hac_g615"]
              --fwd_umi                       Forward UMI (Ftail...UMI...primer) [default: "TTTVVVVTTVVVVTTVVVVTTVVVVTTT"]
              --rev_umi                       Reverse UMI (Rtail...UMI...primer) [default: "AAABBBBAABBBBAABBBBAABBBBAAA"]
              --min_length                    Minimum combined UMI length [default: 40]
              --max_length                    Maximum combined UMI length [default: 60]
              --minimap2_param                 Set the parameters for minimap2 [default: "-ax map-ont -k 13"]
             --include_secondary_reads       Include secondary reads in the analysis [default: false]


         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    

         Example: 
              nextflow run AmstlerStephan/umi-pipeline-nf -r main -profile test,docker
              nextflow run AmstlerStephan/umi-pipeline-nf -r main -c <custom.config> -profile docker 

    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// PRINT VERSION AND EXIT
if(params.version){
    println """\

         ===============================================
          Nextflow UMI amplicon Pipeline
         ===============================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}


// PRINT STANDARD LOGGING INFO
log.info ""
log.info "         ==============================================="
log.info "          U M I - P I P E L I N E - N F"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         ===============================================" }
else {
log.info "         ===============================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         input dir    : ${params.input}"
log.info "         reference    : ${params.reference}"
log.info "         bed          : ${params.bed}"
log.info "         output dir   : ${params.output}"
log.info ""
log.info "         ==============================================="
log.info "         RUN NAME: ${workflow.runName}"
log.info ""


include { UMI_PIPELINE } from './lib/workflows/subworkflows/umi-pipeline.nf'
// include { UMI_PIPELINE } from './lib/workflows/_umi-pipeline.nf'
// include { UMI_PIPELINE_LIVE } from './lib/workflows/umi-pipeline-live.nf'

workflow {
     UMI_PIPELINE()

/*
     if ( params.live ){
          UMI_PIPELINE_LIVE()
     }else {
          UMI_PIPELINE()
     }
*/
}
