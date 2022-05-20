#!/usr/bin/env nextflow
// This file for loading custom functions into the main.nf script (separated for portability)

// FUNCTION TO LOAD DATASETS IN TEST PROFILE
def check_test_data(readPaths, singleEnd) {

    // Set READS testdata
    if( singleEnd ) {
        READS = Channel.from(readPaths)
            .transpose()
            .map { row -> [ file(row[1]).getName().tokenize(".").get(0), file(row[1]) ] }
            .ifEmpty { exit 1, "test profile readPaths was empty - no input files supplied" }
    } else {
        READS = Channel.from(readPaths)
            .map { row -> [ row[0], [file(row[1][0]),file(row[1][1])] ] }
            .ifEmpty { exit 1, "test profile readPaths was empty - no input files supplied" }
    }
    // Return READS channel
    return READS
}
