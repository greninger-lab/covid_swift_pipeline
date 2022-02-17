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
    Usage: 
    An example command for running the pipeline is as follows:
    nextflow run greninger-lab/covid_swift_pipeline -resume -with-docker ubuntu:18.04 --INPUT example/ --OUTDIR example/output/ --PRIMERS v2
    
    Parameters:
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
                        Fastqs should all be gzipped. This can be done with the command gzip *.fastq. [REQUIRED]
        --OUTDIR        Output directory. [REQUIRED]
        --PRIMERS       Primer masterfile to run. By default, this pipeline uses the Swift V2 primers.
                        This pipeline can also use QiaSeq and Artic v3 primers, specified by --PRIMERS qiaseq and --PRIMERS artic respectively.
        --SINGLE_END    Optional flag for single end reads. By default, this pipeline does 
                        paired-end reads.
        --NO_CLIPPING   Skip primerclip option.
        --SGRNA_COUNT   Add extra step to count sgRNAs. 
        --MIN_LEN       Set minimum length for Trimmomatic. Default is 75.
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
params.SINGLE_END = false
params.PRIMERS = false
params.SGRNA_COUNT = false
params.NO_CLIPPING = false
params.MIN_LEN = 75

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
// Use specified primer masterfile
if(params.NO_CLIPPING == false) {
    if (params.PRIMERS == false) {
        println("No primerset specified. Defaulting to Swift V2...")
    }
    else {
        if (params.PRIMERS.toUpperCase() == "QIASEQ") {
            MASTERFILE = file("${baseDir}/sarscov2_qiaseq_masterfile.txt")
            println("Using QiaSeq primerset...")
        }
        else if (params.PRIMERS.toUpperCase() == "ARTIC") {
            MASTERFILE = file("${baseDir}/sarscov2_artic_masterfile.txt")
            println("Using Artic v3 primerset...")
        }
        else {
            MASTERFILE = file("${baseDir}/sarscov2_swift_v2_masterfile.txt")
            println("Using Swift V2 primerset...")
        }
    }
}
else {
    MASTERFILE = file("${baseDir}/sarscov2_swift_v2_masterfile.txt")
    println("--NO_CLIPPING specified. Skipping primer clipping step...")
}

// Print when using SGRNA_COUNT
if (params.SGRNA_COUNT == false) {
    println("--SGRNA_COUNT not specified, skipping counting sgRNAs...")
} else {
    println("--SGRNA_COUNT specified. Will count sgRNAs...")
}

// Setting up files 
REFERENCE_FASTA = file("${baseDir}/NC_045512.2.fasta")
REFERENCE_FASTA_FAI = file("${baseDir}/NC_045512.2.fasta.fai")
TRIM_ENDS=file("${baseDir}/trim_ends.py")
VCFUTILS=file("${baseDir}/vcfutils.pl")
SPLITCHR=file("${baseDir}/splitchr.txt")
ADAPTERS = file("${baseDir}/All_adapters.fa")
FIX_COVERAGE = file("${baseDir}/fix_coverage.py")
PROTEINS = file("${baseDir}/NC_045512_proteins.txt")
AT_REFGENE = file("${baseDir}/annotation/AT_refGene.txt")
AT_REFGENE_MRNA = file("${baseDir}/annotation/AT_refGeneMrna.fa")
LAVA_GFF = file("${baseDir}/annotation/lava_ref.gff")
MAT_PEPTIDES = file("${baseDir}/annotation/mat_peptides.txt")
MAT_PEPTIDE_ADDITION = file("${baseDir}/annotation/mat_peptide_addition.py")
RIBOSOMAL_START = file("${baseDir}/annotation/ribosomal_start.txt")
RIBOSOMAL_SLIPPAGE = file("${baseDir}/annotation/ribosomal_slippage.py")
PROTEINS = file("${baseDir}/annotation/proteins.csv")
CORRECT_AF = file("${baseDir}/annotation/correct_AF.py")
CORRECT_AF_BCFTOOLS = file("${baseDir}/annotation/correct_AF_bcftools.py")
SGRNAS = file("${baseDir}/sgRNAs_60.fasta")
FULL_SGRNAS=file("${baseDir}/sgRNAs.fasta")
FIX_COMPLEX_MUTATIONS = file("${baseDir}/annotation/fix_complex_mutations.py")

// Import processes 
include { Trimming } from './modules.nf'
include { Fastqc } from './modules.nf'
include { Aligning } from './modules.nf'
include { Trimming_SE } from './modules.nf' 
include { Fastqc_SE } from './modules.nf'
include { CountSubgenomicRNAs } from './modules.nf'
include { MapSubgenomics } from './modules.nf'
include { NameSorting } from './modules.nf'
include { Clipping } from './modules.nf'
include { BamSorting } from './modules.nf'
include { GenerateVcf } from './modules.nf'
include { GenerateConsensus } from './modules.nf'
include { PostProcessing } from './modules.nf'
include { AnnotateVariants } from './modules.nf'

// Import reads depending on single end vs. paired end
if(params.SINGLE_END == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.INPUT}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // Looks for gzipped files, assumes all separate samples
    input_read_ch = Channel
        .fromPath("${params.INPUT}*_R1.fastq.gz")
        //.map { it -> [ file(it)]}
        .map { it -> file(it)}
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    // Paired end first few steps
    if(params.SINGLE_END == false) {
        Trimming (
            input_read_ch, 
            ADAPTERS,
            params.MIN_LEN
        )
        Fastqc (
            Trimming.out[1]
        )
        Aligning (
            Trimming.out[0],
            REFERENCE_FASTA
        )
        // Optional step for counting sgRNAs 
        if (params.SGRNA_COUNT != false) {
            CountSubgenomicRNAs (
                Trimming.out[2],
                SGRNAS
            )
            MapSubgenomics (
                CountSubgenomicRNAs.out[1],
                FULL_SGRNAS
            )
        }
    } else {
    // Single end first few steps
        Trimming_SE (
            input_read_ch,
            ADAPTERS,
            params.MIN_LEN
        )
        Fastqc_SE (
            Trimming_SE.out[1]
        )
        Aligning (
            Trimming_SE.out[0],
            REFERENCE_FASTA
        )
        // Optional step for counting sgRNAs 
        if (params.SGRNA_COUNT != false) {
            CountSubgenomicRNAs (
                Trimming_SE.out[2],
                SGRNAS
            )
            MapSubgenomics (
                CountSubgenomicRNAs.out[1],
                FULL_SGRNAS
            )
        }
    }

    // Primerclip options for Swift runs
    if(params.NO_CLIPPING == false) {
        NameSorting (
            Aligning.out[0]
        )
        Clipping (
            NameSorting.out[0],
            MASTERFILE
        )
        GenerateVcf (
            Clipping.out[0],
            REFERENCE_FASTA,
            VCFUTILS,
            REFERENCE_FASTA_FAI,
            SPLITCHR
        )
    } else {
    // Skip primerclip for non-Swift runs
        BamSorting (
            Aligning.out[0]
        )
        
        GenerateVcf (
            BamSorting.out[0],
            REFERENCE_FASTA,
            VCFUTILS,
            REFERENCE_FASTA_FAI,
            SPLITCHR
        )
    }
    GenerateConsensus (
        GenerateVcf.out[0],
        REFERENCE_FASTA,
        REFERENCE_FASTA_FAI,
    )
    PostProcessing (
        GenerateConsensus.out[0],
        REFERENCE_FASTA,
        TRIM_ENDS,
        FIX_COVERAGE,
        REFERENCE_FASTA_FAI
    )
    AnnotateVariants (
        GenerateVcf.out[1],
        MAT_PEPTIDES,
        MAT_PEPTIDE_ADDITION,
        RIBOSOMAL_SLIPPAGE,
        RIBOSOMAL_START,
        PROTEINS,
        AT_REFGENE,
        AT_REFGENE_MRNA,
        CORRECT_AF_BCFTOOLS,
        FIX_COMPLEX_MUTATIONS
    )
}
