process IvarConsensus {
    container "quay.io/biocontainers/ivar:1.3.1--hecb563c_3"
    errorStrategy 'retry'
    maxRetries 2

    input:
        file(BAMFILE)

    output:
        file("*_taylor.fasta")
    
    publishDir params.OUTDIR, mode: 'copy', pattern: '*.fasta'

    shell:
    '''
    #!/bin/bash
    
    base=$(echo !{BAMFILE} | cut -d. -f1)
    samtools mpileup -d 5000 -A -Q 0 !{BAMFILE} | ivar consensus -p ${base} -n 'N' -m 50 -t 0.25 -i ${base}
    cp ${base}.fa ${base}_taylor.fasta

    '''
}
