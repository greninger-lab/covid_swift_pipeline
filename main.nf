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
nextflow.preview.dsl=2

// Print help message
def helpMessage() {
    log.info"""
    This is a short pipeline made to redo consensuses with more IUPAC output.
    
    Parameters:
        --INPUT         Input folder where all bams to redo are located. [REQUIRED]
        --OUTDIR        Output directory. [REQUIRED]
        -with-docker ubuntu:18.04   [REQUIRED]
        -resume [RECOMMENDED]
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

// Setting up files 
// REFERENCE_FASTA = file("${baseDir}/NC_045512.2.fasta")
// REFERENCE_FASTA_FAI = file("${baseDir}/NC_045512.2.fasta.fai")
// VCFUTILS=file("${baseDir}/vcfutils.pl")
// SPLITCHR=file("${baseDir}/splitchr.txt")
// FIX_COVERAGE = file("${baseDir}/fix_coverage.py")
// PROTEINS = file("${baseDir}/NC_045512_proteins.txt")
// AT_REFGENE = file("${baseDir}/annotation/AT_refGene.txt")
// AT_REFGENE_MRNA = file("${baseDir}/annotation/AT_refGeneMrna.fa")
// LAVA_GFF = file("${baseDir}/annotation/lava_ref.gff")
// MAT_PEPTIDES = file("${baseDir}/annotation/mat_peptides.txt")
// MAT_PEPTIDE_ADDITION = file("${baseDir}/annotation/mat_peptide_addition.py")
// RIBOSOMAL_START = file("${baseDir}/annotation/ribosomal_start.txt")
// RIBOSOMAL_SLIPPAGE = file("${baseDir}/annotation/ribosomal_slippage.py")
// PROTEINS = file("${baseDir}/annotation/proteins.csv")
// CORRECT_AF = file("${baseDir}/annotation/correct_AF.py")
// CORRECT_AF_BCFTOOLS = file("${baseDir}/annotation/correct_AF_bcftools.py")
// SGRNAS = file("${baseDir}/sgRNAs_60.fasta")
// FULL_SGRNAS=file("${baseDir}/sgRNAs.fasta")
// FIX_COMPLEX_MUTATIONS = file("${baseDir}/annotation/fix_complex_mutations.py")

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
