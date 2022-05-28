# TAYLOR (Trimming Amplicons You LOve Rapidly)
This pipeline takes gzipped fastq files and outputs .bam files aligned to NC_045512.2 as well as consensus fastas. Notably this pipeline incorporates [primerclip](https://github.com/swiftbiosciences/primerclip/tree/deltest), and can handle both Swift and QiaSeq primersets. Without specifying any additional options, default input files are paired-end fastq files that cover the 116-255 bp amplicons produced from the Swift Amplicon SARS-CoV-2 Panel. `--SINGLE_END` can also be specified for single end reads. This pipeline can also be run without the primerclip option by specifying `--NO_CLIPPING` for consensus generation of non-Swift or non-QiaSeq SARS-CoV-2 samples. 

## Installation

1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation).
   - Make sure you move nextflow to a directory in your PATH variable.
2. Install [docker](https://docs.docker.com/get-docker/).

## Usage
| Command  | Description |
| ---      | ---         | 
| --INPUT  | __(Required)__ Input folder where gzipped fastqs are located. For current  directory, `./` can be used.
| --OUTDIR | __(Required)__ Output folder where .bams and consensus fastas will be piped into.
| --SINGLE_END | (Optional) Flag to indicate input reads are single end. By default this pipeline expects paired end reads.
| --NO_CLIPPING | (Optional) Skip primerclip option for shotgun/CovidSeq runs.
| --SGRNA_COUNT | (Optional) Add extra step to count sgRNAs.
| --PRIMERS | (Optional) Specify which primerset to use (e.g. `--PRIMERS qiaseq`). Default: Swift V2. Options: `qiaseq`, `artic_v3`, `artic_v4`, `artic_v4.1`.
| --MIN_LEN | (Optional) Set minimum length for trimming. Default: 75.
| -with-docker ubuntu:18.04 | __(Required)__ Runs command with Ubuntu docker.
| -resume  | __(Recommended)__ nextflow will pick up where it left off if the previous command was interrupted for some reason.
| -with-trace | __(Recommended)__ Outputs a trace.txt that shows which processes end up in which work/ directories. 
| -with-report | __(Recommended)__ Outputs a report.html that gives basic stats and work directories for each process.
| -latest | __(Recommended)__ Pulls the most recent github version.
| -profile | __(Recommended)__ Picks which profile in nextflow.config to run (e.g. `-profile cloud_big`). If running on AWS, recommended to run with `profile cloud_big` and for more memory-intensive runs, with `profile cloud_bigger`).

Example paired fastqs are provided in the example/ folder. These can be run with the command:
- Example command for example fastqs: ```nextflow run greninger-lab/covid_swift_pipeline -latest --INPUT example/ --OUTDIR output/ --PAIRED_END -resume -with-trace with-docker ubuntu:18.04 --PRIMERS v2```
