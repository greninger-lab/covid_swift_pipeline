// Generate pileups from input bams.
process GeneratePileup { 
    container "quay.io/greninger-lab/swift-pipeline:latest"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 2

    input:
        file(BAMFILE) // from input_read_ch
        file REFERENCE_FASTA
        file VCFUTILS
        file REFERENCE_FASTA_FAI
        file SPLITCHR
    output:
        tuple env(base),file(BAMFILE),file("*.mpileup")
    
    publishDir params.OUTDIR, mode: 'copy', pattern: '*.mpileup'

    shell:
    '''
    #!/bin/bash
    ls -latr

    R1=$(echo !{BAMFILE} | cut -d. -f1)
    base=$R1
    echo $R1

    # Index bam file
    /usr/local/miniconda/bin/samtools index !{BAMFILE} 

    # Parallelize pileup based on number of cores
    splitnum=$(($((29903/!{task.cpus}))+1))
    perl !{VCFUTILS} splitchr -l $splitnum !{REFERENCE_FASTA_FAI} | \\
    #cat !{SPLITCHR} | \\
        xargs -I {} -n 1 -P !{task.cpus} sh -c \\
            "/usr/local/miniconda/bin/samtools mpileup \\
                -f !{REFERENCE_FASTA} -r {} -aa -Q 0 \\
                --no-BAQ \\
                --count-orphans \\
                -d 5000 \\
            !{BAMFILE} > tmp.{}.pileup"
    
    mv tmp.NC_045512.2\\:1-* \${R1}_catted.pileup
    for file in tmp*.pileup; do grep -v "#" $file >> \${R1}_catted.pileup; done
    sort -k2 -n \${R1}_catted.pileup > \${R1}.sorted.mpileup
    '''
}

process IvarConsensus {
    container "quay.io/biocontainers/ivar:1.3.1--hecb563c_3"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 2

    input:
        tuple val(base),file(BAMFILE),file(PILEUP)

    output:
        tuple val(base),file("${base}_taylor.fasta")
    
    publishDir params.OUTDIR, mode: 'copy', pattern: '*.fasta'

    shell:
    '''
    #!/bin/bash
    ls -latr

    pileupsize=$(($(wc -c !{PILEUP} | awk '{print $1'})+0))
    if (( $pileupsize > 92 ))
    then
        cat !{PILEUP} | /usr/local/bin/ivar consensus -p !{base} -n N -m 50 -t 0.25 -i !{base}
        cp !{base}.fa !{base}_taylor.fasta
    else
       printf '>!{base}\n' > !{base}_taylor.fasta
       printf 'n%.0s' {1..29539} >> !{base}_taylor.fasta
    fi
    
    '''
}
