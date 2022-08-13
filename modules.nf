process IvarConsensus {
    container "quay.io/biocontainers/ivar:1.3.1--hecb563c_3"
    errorStrategy 'retry'
    maxRetries 2

    input:
        file(BAMFILE)

    output:
        tuple file("*_taylor.fasta"), file("*.bed")
    
    publishDir "${params.OUTDIR}ivar_consensus", mode: 'copy', pattern: '*.fasta'
    publishDir "${params.OUTDIR}coverage_bed", mode: 'copy', pattern: '*.bed'

    shell:
    '''
    #!/bin/bash
    
    base=$(echo !{BAMFILE} | cut -d. -f1)
    
    # generate coverage bed file
    samtools depth -a -H !{BAMFILE} -o ${base}.bed

    # call consensus genome 
    samtools mpileup -d 5000 -A -Q 0 !{BAMFILE} | ivar consensus -p ${base} -n 'N' -m 50 -t 0.25 -i ${base}
    cp ${base}.fa ${base}_taylor.fasta

    '''
}
