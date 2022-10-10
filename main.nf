#!/usr/bin/env nextflow

// ENABLE DSL2
nextflow.enable.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\


         Usage: 
              nextflow run umi-pipeline-nf [OPTIONS]...

         Options: GENERAL
              --input [path/to/input/dir]     [REQUIRED] Input directory containing (zipped) FASTQ files

              --reference [path/to/ref.fa]    [REQUIRED] Path to the reference genome in fasta format

              --bed [path/to/data.bed]        [REQUIRED] Path to the bed file

              --output [STR]                  A string that can be given to name the output directory [default: "umi-pipeline-nf_results"]

              --threads                       Number of maximum threads to use [default: availableProcessors -1]

         Options: ADVANCED
              --umi_errors                    Max differences between UMI in read and UMI pattern [default: 3]

              --min_reads_per_cluster         Min number of reads required for a consensus read [default: 20]

              --max_reads_per_cluster         Max number of 1D used for a consensus read [default: 60]

              --min_overlap                   Min overlap with target region [default: 0.90]

              --balance_strands               Balance forward and reverse 1D reads in clusters [default: true]

              --medaka_model                  Medaka model used to compute consensus reads [default: "r941_min_high_g360"]

              --fwd_universal_primer          Forward tail of primer (Ftail...UMI...primer) [default: "GTATCGTGTAGAGACTGCGTAGG"]

              --rev_universal_primer          Reverse tail of primer (Rtail...UMI...primer) [default: "AGTGATCGAGTCAGTGCGAGTG"]

              --fwd_umi                       Forward UMI (Ftail...UMI...primer) [default: "TTTVVVVTTVVVVTTVVVVTTVVVVTTT"]

              --rev_umi                       Reverse UMI (Rtail...UMI...primer) [default: "AAABBBBAABBBBAABBBBAABBBBAAA"]

              --min_length                    Minimum combined UMI length [default: 40]

              --max_length                    Maximum combined UMI length [default: 60]

              --minimap_param                 Set the parameters for minimap2 [default: "-ax map-ont -k 13"]

         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    

         Example: 
              nextflow run umi-pipeline-nf \
              --input /path/to/input/dir --reference /path/to/genome.fa \
              --bed [path/to/data.bed] --output [STR]

              nextflow run umi-pipeline-nf \
              --project test --input data/example_egfr_single_cluster.fastq \
              --reference data/example_egfr_reference.fasta --bed data/example_egfr_amplicon.bed

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
log.info "          E C S E Q - t e m p l a t e   P I P E L I N E"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         ===============================================" }
else {
log.info "         ===============================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         input dir    : ${workflow.profile.tokenize(",").contains("test") ? "-" : "${params.input}"}"
log.info "         reference    : ${params.reference}"
log.info "         output dir   : ${params.output}"
log.info ""
log.info "         ==============================================="
log.info "         RUN NAME: ${workflow.runName}"
log.info ""


include { UMI_PIPELINE } from './lib/workflows/umi-pipeline.nf'

workflow {
    UMI_PIPELINE()
}