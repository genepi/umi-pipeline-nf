nextflow.enable.dsl = 2

requiredParams = [
    'input', 'reference', 'bed'
]

for (param in requiredParams) {
    if (params[param] == null) {
      exit 1, "Parameter ${param} is required."
    }
}

// DEFINE PATHS
bed = file("${params.bed}", checkIfExists: true)
reference = file("${params.reference}", checkIfExists: true)

// python scripts
umi_filter_reads = file( "${projectDir}/bin/filter_reads.py", checkIfExists: true)
umi_extract = file( "${projectDir}/bin/extract_umis.py", checkIfExists: true)

// STAGE CHANNELS
fastq_files_ch = Channel.fromPath("${params.input}")

////////////////////
// BEGIN PIPELINE //
////////////////////

// INCLUDES # here you must give the relevant process files from the lib directory 
include {COPY_BED} from '../processes/copy_bed.nf'
include {MAP_1D} from '../processes/map_1d.nf'
include {SPLIT_READS} from  '../processes/split_reads.nf'
include {DETECT_UMI_FASTA} from '../processes/detect_umi_fasta.nf'

// SUB-WORKFLOWS
workflow UMI_PIPELINE {

    main:

        COPY_BED( bed )
        MAP_1D( fastq_files_ch, reference )
        SPLIT_READS( MAP_1D.out.bam_1d, MAP_1D.out.bai_1d, COPY_BED.out.bed, umi_filter_reads )
        DETECT_UMI_FASTA( SPLIT_READS.out.split_reads_fastq, umi_extract )

}


//////////////////
// END PIPELINE //
//////////////////

// WORKFLOW TRACING # what to display when the pipeline finishes
// eg. with errors
workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

// eg. in general
workflow.onComplete {

    log.info ""
    log.info "         Pipeline execution summary"
    log.info "         ---------------------------"
    log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
    log.info "         Profile      : ${workflow.profile}"
    log.info "         Launch dir   : ${workflow.launchDir}"    
    log.info "         Work dir     : ${workflow.workDir} ${!params.debug && workflow.success ? "(cleared)" : "" }"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    // run a small clean-up script to remove "work" directory after successful completion 
    if (!params.debug && workflow.success) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}
