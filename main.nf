#!/usr/bin/env nextflow
/*
=======================================================================
                        TAYLOR
=======================================================================
 Trimming Amplications You LOve Rapidly
 #### Homepage / Documentation
https://github.com/greninger-lab/covid_swift_pipeline
-----------------------------------------------------------------------
*/

// Using the Nextflow DSL-2 to account for the logic flow of this workflow

// Print help message
def helpMessage() {
    log.info"""
    This is a short pipeline made to redo consensuses with more IUPAC output.
    
    Parameters:
        --INPUT         Input folder where all bams to redo are located. [REQUIRED]
        --OUTDIR        Output directory. [REQUIRED]
        -with-docker    ubuntu:18.04   [REQUIRED]
        -resume         [RECOMMENDED]
        -profile        Specify which profile to run. For AWS, run with -profile cloud_big. For large memory-intensive runs on AWS, run with -profile cloud_bigger.
        
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Check Nextflow version for enabling DSL2
nextflow_dsl2_v = '20.07.1'
if ( nextflow.version.matches(">= $nextflow_dsl2_v") ) {
    nextflow.enable.dsl=2
} else {
    nextflow.preview.dsl=2
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Initializing flags
params.INPUT = false
params.OUTDIR= false

// Checking for argument validity
// Throw error if --INPUT not set
if (params.INPUT == false) {
    println( "Must provide an input directory with --INPUT") 
    exit(1)
}
// Make sure INPUT ends with trailing slash
if (!params.INPUT.endsWith("/")){
    println("Make sure your input directory ends with trailing slash.")
   exit(1)
}
// Throw error if --OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
   println("Make sure your output directory ends with trailing slash.")
   exit(1)
}

// Import processes 
include { IvarConsensus } from './modules.nf'

// Import bams from input folder
input_read_ch = Channel
    .fromPath("${params.INPUT}*.bam")
    .map { it -> file(it)}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    IvarConsensus (
        input_read_ch
    )
}
