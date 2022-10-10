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
test = file("${baseDir}/data/test.txt", checkIfExists : true)

// STAGE CHANNELS

/*
// STAGE BAM FILES FROM TEST PROFILE # this establishes the test data to use with -profile test
if ( workflow.profile.tokenize(",").contains("test") ){

        include { check_test_data } from './lib/functions.nf' params(readPaths: params.readPaths, singleEnd: params.SE)
        INPUT = check_test_data(params.readPaths, params.SE)

} else {

    // STAGE INPUT CHANNELS # this defines the normal input when test profile is not in use
    INPUT = Channel.fromPath(params.input)

}
*/

fastq_files_ch = Channel.fromPath(params.input)

////////////////////
// BEGIN PIPELINE //
////////////////////

/*
 *   Workflows are where you define how different processes link together. They
 *    may be modularised into "sub-workflows" which must be named eg. 'DNAseq'
 *    and there must always be one MAIN workflow to link them together, which
 *    is always unnamed.
 */

// INCLUDES # here you must give the relevant process files from the lib directory 
include {HELLO_WORLD} from '../processes/Hello_World.nf'
include {COPY_BED} from '../processes/copy_bed.nf'

// SUB-WORKFLOWS
workflow UMI_PIPELINE {

    // here we define the structure of our workflow i.e. how the different processes lead into each other
    // eg. process(input1, input2, etc.)
    // eg. process.out[0], process.out[1], etc.
    // index numbers [0],[1],etc. refer to different outputs defined for processes in process.nf
    // You can forgo the index number if there is only 1 output.
    // ALWAYS PAY ATTENTION TO CARDINALITY!!

    main:
        //HELLO_WORLD( test )
        COPY_BED( bed )
        
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
