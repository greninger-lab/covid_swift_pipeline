process IvarConsensus {
    container "greningerlab/revica:latest"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 2

    input:
        file(BAMFILE)

    output:
        tuple file("*_taylor.fasta")
    
    publishDir params.OUTDIR, mode: 'copy', pattern: '*.fasta'

    shell:
    '''
    #!/bin/bash
    
    base=$(echo !{BAMFILE} | cut -d. -f1)
    echo "${base}...testtest" > test.txt
    
    samtools mpileup -d 5000 -A -Q 0 !{BAMFILE} | ivar consensus -p ${base} -n 'N' -m 50 -t 0.25 -i ${base}
    cp ${base}.fa ${base}_taylor.fasta

    '''
}
