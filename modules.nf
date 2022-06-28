// Use Trimmomatic to trim files, above Q20, minlen of 75
// Initialize summary file and input trimming stats into summary file
process Trimming { 
    container "quay.io/biocontainers/trimmomatic:0.35--6"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(base), file(R1), file(R2) // from input_read_ch
        file ADAPTERS
        val MINLEN
    output: 
        tuple val(base), file("${base}.trimmed.fastq.gz"),file("${base}_summary.csv") //into Trim_out_ch
        tuple val(base), file(R1),file(R2),file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz"),file("${base}.R1.unpaired.fastq.gz"), file("${base}.R2.unpaired.fastq.gz") //into Trim_out_ch2
        tuple val(base), file("${base}.trimmed.fastq.gz") //into Trim_out_ch3

    publishDir "${params.OUTDIR}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'

    script:
    """
    #!/bin/bash

    trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:${MINLEN}

    num_r1_untrimmed=\$(gunzip -c ${R1} | wc -l)
    num_r2_untrimmed=\$(gunzip -c ${R2} | wc -l)
    num_untrimmed=\$((\$((num_r1_untrimmed + num_r2_untrimmed))/4))

    num_r1_paired=\$(gunzip -c ${base}.R1.paired.fastq.gz | wc -l)
    num_r2_paired=\$(gunzip -c ${base}.R2.paired.fastq.gz | wc -l)
    num_paired=\$((\$((num_r1_paired + num_r2_paired))/4))

    num_r1_unpaired=\$(gunzip -c ${base}.R1.unpaired.fastq.gz | wc -l)
    num_r2_unpaired=\$(gunzip -c ${base}.R2.unpaired.fastq.gz | wc -l)
    num_unpaired=\$((\$((num_r1_unpaired + num_r2_unpaired))/4))

    num_trimmed=\$((num_paired + num_unpaired))
    
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    
    echo Sample_Name,Raw_Reads,Trimmed_Paired_Reads,Trimmed_Unpaired_Reads,Total_Trimmed_Reads,Percent_Trimmed,Mapped_Reads,Clipped_Mapped_Reads,Mean_Coverage,Spike_Mean_Coverage,Spike_100X_Cov_Percentage,Spike_200X_Cov_Percentage,Lowest_Spike_Cov,Percent_N > ${base}_summary.csv
    printf "${base},\$num_untrimmed,\$num_paired,\$num_unpaired,\$num_trimmed,\$percent_trimmed" >> ${base}_summary.csv

    cat *paired.fastq.gz > ${base}.trimmed.fastq.gz
    
    """
}

// Fastqc our fastq files for quick sanity check
process Fastqc {
    container "quay.io/biocontainers/fastqc:0.11.9--0"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file(R1),file(R2),file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz"),file("${base}.R1.unpaired.fastq.gz"), file("${base}.R2.unpaired.fastq.gz") // from Trim_out_ch2
    output: 
    file("*fastqc*") //into Fastqc_ch 

    publishDir "${params.OUTDIR}fastqc", mode: 'copy'

    script:
    """
    #!/bin/bash

    /usr/local/bin/fastqc ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R2.paired.fastq.gz

    """
}


// Use Trimmomatic to trim files, above Q20, minlen of 75
// Initialize summary file and input trimming stats into summary file
process Trimming_SE { 
    container "quay.io/biocontainers/trimmomatic:0.35--6"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        file R1 //from input_read_ch
        file ADAPTERS
        val MINLEN
    output: 
        tuple env(base),file("*.trimmed.fastq.gz"),file("*summary.csv") //into Trim_out_ch_SE
        tuple env(base),file("*.trimmed.fastq.gz") //into Trim_out_ch2_SE
        tuple env(base),file("*.trimmed.fastq.gz") //into Trim_out_ch3_SE

    publishDir "${params.OUTDIR}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'

    script:
    """
    #!/bin/bash

    #base=`basename ${R1} ".fastq.gz"`
    base=\$(echo ${R1} | awk -F'_R1' '{print \$1}')

    echo \$base

    trimmomatic SE -threads ${task.cpus} ${R1} \$base.trimmed.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:${MINLEN}

    num_untrimmed=\$((\$(gunzip -c ${R1} | wc -l)/4))
    num_trimmed=\$((\$(gunzip -c \$base'.trimmed.fastq.gz' | wc -l)/4))
    
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    
    echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Mapped_Reads,Clipped_Mapped_Reads,Mean_Coverage,Spike_Mean_Coverage,Spike_100X_Cov_Percentage,Spike_200X_Cov_Percentage,Lowest_Spike_Cov,Percent_N > \$base'_summary.csv'
    printf "\$base,\$num_untrimmed,\$num_trimmed,\$percent_trimmed" >> \$base'_summary.csv'
    
    ls -latr
    """
}

// Fastqc our fastq files for quick sanity check
process Fastqc_SE {
    container "quay.io/biocontainers/fastqc:0.11.9--0"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(base),file("${base}.trimmed.fastq.gz") // from Trim_out_ch2_SE
    output: 
        file("*fastqc*") //into Fastqc_ch 

    publishDir "${params.OUTDIR}fastqc", mode: 'copy'

    script:
    """
    #!/bin/bash

    /usr/local/bin/fastqc ${base}.trimmed.fastq.gz
    """
}

// Align fastq files to Wuhan refseq using bbmap
process Aligning {
    container "quay.io/biocontainers/bbmap:38.86--h1296035_0"
    //container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7	"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz"),file("${base}_summary.csv")
        file REFERENCE_FASTA
        val BBMAP_INDEL_SKIP
    output:
        tuple val(base), file("${base}.bam"),file("${base}_summary2.csv") //into Aligned_bam_ch
        tuple val (base), file("*") //into Dump_ch

    script:
    """
    #!/bin/bash

    cat ${base}*.fastq.gz > ${base}_cat.fastq.gz
    /usr/local/bin/bbmap.sh in=${base}_cat.fastq.gz outm=${base}.bam ref=${REFERENCE_FASTA} local=true ${BBMAP_INDEL_SKIP} -Xmx6g > bbmap_out.txt 2>&1
    reads_mapped=\$(cat bbmap_out.txt | grep "mapped:" | cut -d\$'\\t' -f3)

    cp ${base}_summary.csv ${base}_summary2.csv
    printf ",\$reads_mapped" >> ${base}_summary2.csv

    """
}

// Optional step for counting sgRNAs
process CountSubgenomicRNAs {
    container "quay.io/biocontainers/bbmap:38.86--h1296035_0"

    // Retry on fail at most three times
    errorStrategy 'retry'
    maxRetries 3

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz") //from Trim_out_ch3
        file SGRNAS 
    output:
        file("*stats*")
        tuple val(base),file("*_sgrnas.fastq.gz")
    
    publishDir "${params.OUTDIR}sgRNAs", mode: 'copy'

    script:
    """
    #!/bin/bash
    /usr/local/bin/bbduk.sh in=${base}.trimmed.fastq.gz outm=${base}_sgrnas.fastq.gz ref=${SGRNAS} stats=${base}_sgrnas_stats.txt refstats=${base}_sgrnas_refstats.txt k=40 qhdist=1 -Xmx6g

    """
}

// Optional step for counting sgRNAs
process MapSubgenomics {
    container "quay.io/biocontainers/bbmap:38.86--h1296035_0"

    // Retry on fail at most three times
    errorStrategy 'retry'
    maxRetries 3

    input: 
        tuple val(base), file("${base}_sgrnas.fastq.gz") //from CountSubgenomics
        file FULL_SGRNAS 
    output:
        file("${base}_sgrnas_mapped.bam")
    
    publishDir "${params.OUTDIR}sgRNAs", mode: 'copy', pattern: '*.bam'

    script:
    """
    #!/bin/bash

    /usr/local/bin/bbmap.sh in=${base}_sgrnas.fastq.gz outm=${base}_sgrnas_mapped.bam ref=${FULL_SGRNAS} -Xmx6g 2>&1

    """
}


// Sort sam for input into primerclip
process NameSorting { 
    container "quay.io/biocontainers/samtools:1.3--h0592bc0_3"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val (base), file("${base}.bam"),file("${base}_summary2.csv") //from Aligned_bam_ch
    output:
        tuple val (base), file("${base}.sorted.sam"),file("${base}_summary2.csv") //into Sorted_sam_ch

    publishDir "${params.OUTDIR}inprogress_summary", mode: 'copy', pattern: '*summary.csv'

    script:
    """
    #!/bin/bash
    samtools sort -@ ${task.cpus} -n -O sam ${base}.bam > ${base}.sorted.sam

    """
}

// Use primerclip to trim Swift primers
process Clipping { 
    container "quay.io/greninger-lab/swift-pipeline:latest"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val (base), file("${base}.sorted.sam"),file("${base}_summary2.csv") //from Sorted_sam_ch
        file MASTERFILE
    output:
        tuple val (base), file("${base}.clipped.bam"), file("${base}.clipped.bam.bai"),file("${base}_summary3.csv"),env(bamsize) //into Clipped_bam_ch
        tuple val (base), file("${base}.clipped.bam"), file("${base}.clipped.bam.bai"),env(bamsize) //into Clipped_bam_ch2
        tuple val (base), file("${base}.clipped.bam"), file("${base}.clipped.bam.bai"),env(bamsize) //into Clipped_bam_ch3

    publishDir params.OUTDIR, mode: 'copy', pattern: '*.clipped.bam'
    publishDir "${params.OUTDIR}inprogress_summary", mode: 'copy', pattern: '*summary3.csv'

    script:
        """
        #!/bin/bash
        ls -latr
        /./root/.local/bin/primerclip -s ${MASTERFILE} ${base}.sorted.sam ${base}.clipped.sam
        #/usr/local/miniconda/bin/samtools sort -@ ${task.cpus} -n -O sam ${base}.clipped.sam > ${base}.clipped.sorted.sam
        #/usr/local/miniconda/bin/samtools view -@ ${task.cpus} -Sb ${base}.clipped.sorted.sam > ${base}.clipped.unsorted.bam
        #/usr/local/miniconda/bin/samtools sort -@ ${task.cpus} -o ${base}.clipped.unsorted.bam ${base}.clipped.bam
        /usr/local/miniconda/bin/samtools sort -@ ${task.cpus} ${base}.clipped.sam -o ${base}.clipped.bam
        /usr/local/miniconda/bin/samtools index ${base}.clipped.bam
        clipped_reads=\$(/usr/local/miniconda/bin/samtools flagstat ${base}.clipped.bam | grep "mapped (" | awk '{print \$1}')
        echo "clipped reads: \$clipped_reads"
        /usr/local/miniconda/bin/bedtools genomecov -d -ibam ${base}.clipped.bam > ${base}_coverage.txt
        meancoverage=\$(cat ${base}_coverage.txt | awk '{sum+=\$3} END { print sum/NR}')
        bamsize=\$((\$(wc -c ${base}.clipped.bam | awk '{print \$1'})+0))
        echo "bamsize: \$bamsize"
        if (( \$bamsize > 92 ))
        then
            # Spike protein coverage
            awk '\$2 ~ /21563/,\$2 ~ /25384/' ${base}_coverage.txt > ${base}_spike_coverage.txt
            avgcoverage=\$(cat ${base}_spike_coverage.txt | awk '{sum+=\$3} END { print sum/NR}')
            proteinlength=\$((25384-21563+1))
            cov100=\$((100*\$(cat ${base}_spike_coverage.txt | awk '\$3>=100' | wc -l)/3822))
            cov200=\$((100*\$(cat ${base}_spike_coverage.txt | awk '\$3>=200' | wc -l)/3822))
            mincov=\$(sort -nk 3 ${base}_spike_coverage.txt | head -n 1 | cut -f3)
        else
            avgcoverage=0
            cov100=0
            cov200=0
            mincov=0
        fi
        
        cp ${base}_summary2.csv ${base}_summary3.csv
        printf ",\$clipped_reads,\$meancoverage,\$avgcoverage,\$cov100,\$cov200,\$mincov" >> ${base}_summary3.csv
        """
    } 

// Sort bam
process BamSorting { 
    container "quay.io/greninger-lab/swift-pipeline:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
      tuple val (base), file("${base}.bam"),file("${base}_summary2.csv") //from Aligned_bam_ch
    output:
      tuple val (base), file("${base}.sorted.bam"),file("${base}.sorted.bam.bai"),file("${base}_summary3.csv"),env(bamsize) //into Clipped_bam_ch
    
    publishDir "${params.OUTDIR}inprogress_summary", mode: 'copy', pattern: '*summary.csv'
    publishDir params.OUTDIR, mode: 'copy', pattern: '*.sorted.bam'

    script:
    """
    #!/bin/bash
    /usr/local/miniconda/bin/samtools sort -@ ${task.cpus} ${base}.bam > ${base}.sorted.bam
    /usr/local/miniconda/bin/samtools index ${base}.sorted.bam

    clipped_reads=0
    /usr/local/miniconda/bin/bedtools genomecov -d -ibam ${base}.sorted.bam > ${base}_coverage.txt
    meancoverage=\$(cat ${base}_coverage.txt | awk '{sum+=\$3} END { print sum/NR}')
    bamsize=\$((\$(wc -c ${base}.sorted.bam | awk '{print \$1'})+0))
    echo "bamsize: \$bamsize"
    if (( \$bamsize > 92 ))
    then
        # Spike protein coverage
        awk '\$2 ~ /21563/,\$2 ~ /25384/' ${base}_coverage.txt > ${base}_spike_coverage.txt
        avgcoverage=\$(cat ${base}_spike_coverage.txt | awk '{sum+=\$3} END { print sum/NR}')
        proteinlength=\$((25384-21563+1))
        cov100=\$((100*\$(cat ${base}_spike_coverage.txt | awk '\$3>=100' | wc -l)/3822))
        cov200=\$((100*\$(cat ${base}_spike_coverage.txt | awk '\$3>=200' | wc -l)/3822))
        mincov=\$(sort -nk 3 ${base}_spike_coverage.txt | head -n 1 | cut -f3)
    else
        avgcoverage=0
        cov100=0
        cov200=0
        mincov=0
    fi

    cp ${base}_summary2.csv ${base}_summary3.csv
    printf ",\$clipped_reads,\$meancoverage,\$avgcoverage,\$cov100,\$cov200,\$mincov" >> ${base}_summary3.csv


    """
}

// Generate final consensus from pileup from bam.
process GenerateVcf {
    container "quay.io/greninger-lab/swift-pipeline:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val (base), file(BAMFILE),file(INDEX_FILE),file("${base}_summary3.csv"),val(bamsize) //from Clipped_bam_ch
        file REFERENCE_FASTA
        file VCFUTILS
        file REFERENCE_FASTA_FAI
        file SPLITCHR
    output:
        tuple val(base), val(bamsize), file(BAMFILE),file(INDEX_FILE),file("${base}_pre_bcftools.vcf"), file("${base}.vcf.gz"),file("${base}.vcf.gz.tbi"),file("${base}_summary3.csv")
        tuple val(base),val(bamsize),file("${base}_pre_bcftools.vcf")
    
    publishDir params.OUTDIR, mode: 'copy', pattern: '*_pre_bcftools.vcf'

    shell:
    '''
    #!/bin/bash
    ls -latr

    R1=!{base}

    echo "bamsize: !{bamsize}"

    #if [ -s !{BAMFILE} ]
    # More reliable way of checking bam size, because of aliases
    if (( !{bamsize} > 92 ))
    then
        # Parallelize pileup based on number of cores
        splitnum=$(($((29903/!{task.cpus}))+1))
        perl !{VCFUTILS} splitchr -l $splitnum !{REFERENCE_FASTA_FAI} | \\
        #cat !{SPLITCHR} | \\
            xargs -I {} -n 1 -P !{task.cpus} sh -c \\
                "/usr/local/miniconda/bin/bcftools mpileup \\
                    -f !{REFERENCE_FASTA} -r {} \\
                    --count-orphans \\
                    --no-BAQ \\
                    --max-depth 50000 \\
                    --max-idepth 500000 \\
                    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
                !{BAMFILE} | /usr/local/miniconda/bin/bcftools call -A -m -Oz - > tmp.{}.vcf.gz"
        
        # Concatenate parallelized vcfs back together
        gunzip tmp*vcf.gz
        mv tmp.NC_045512.2\\:1-* \${R1}_catted.vcf
        for file in tmp*.vcf; do grep -v "#" $file >> \${R1}_catted.vcf; done

        cat \${R1}_catted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | /usr/local/miniconda/bin/bcftools norm -m -any > \${R1}_pre_bcftools.vcf
        
        # Make sure variants are majority variants for consensus calling
        #/usr/local/miniconda/bin/bcftools filter -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)' --threads !{task.cpus} \${R1}_pre_bcftools.vcf -o \${R1}.vcf
        #/usr/local/miniconda/bin/bcftools filter -e 'IMF < 0.5' \${R1}_pre2.vcf -o \${R1}.vcf
	
        /usr/local/miniconda/bin/bcftools filter -i 'IMF > 0.5 || (DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)' --threads !{task.cpus} \${R1}_pre_bcftools.vcf -o \${R1}_pre2.vcf
        /usr/local/miniconda/bin/bcftools norm --check-ref s --fasta-ref !{REFERENCE_FASTA} -Ov \${R1}_pre2.vcf > \${R1}_pre3.vcf

        # pull out header
        grep "#" \${R1}_pre3.vcf > \${R1}.vcf
        # get rid of long sgRNAs that are called due to fake depths
        grep -v "#" \${R1}_pre3.vcf | awk -F'\t' 'length($4) <60 { print }' >> \${R1}.vcf

        # Index and generate consensus from vcf with majority variants
        /usr/local/miniconda/bin/bgzip \${R1}.vcf
        /usr/local/miniconda/bin/tabix \${R1}.vcf.gz
    else
        touch \${R1}_pre_bcftools.vcf
        touch \${R1}.vcf.gz
    fi

    '''
}

// Generate final consensus from vcf.
process GenerateConsensus {
    container "quay.io/biocontainers/bcftools:1.14--hde04aa1_1"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(base), val(bamsize), file(BAMFILE),file(INDEX_FILE),file("${base}_pre_bcftools.vcf"), file("${base}.vcf.gz"),file("${base}.vcf.gz.tbi"),file("${base}_summary3.csv")
        file REFERENCE_FASTA
        file REFERENCE_FASTA_FAI
    output:
        tuple val(base), val(bamsize), file("${base}_pre_bcftools.vcf"), file("${base}.vcf.gz"),file("${base}.vcf.gz.tbi"),file("${base}_summary3.csv"),file("${base}.consensus.fa"),file(BAMFILE),file(INDEX_FILE)

    script:
    """
        cat ${REFERENCE_FASTA} | /usr/local/bin/bcftools consensus ${base}.vcf.gz > ${base}.consensus.fa
    """
}

// Process consensus and annotate variants.
process PostProcessing {
    container "quay.io/greninger-lab/swift-pipeline:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(base), val(bamsize), file("${base}_pre_bcftools.vcf"), file("${base}.vcf.gz"),file("${base}.vcf.gz.tbi"),file("${base}_summary3.csv"),file("${base}.consensus.fa"),file(BAMFILE),file(INDEX_FILE)
        file REFERENCE_FASTA
        file TRIM_ENDS
        file FIX_COVERAGE
        file REFERENCE_FASTA_FAI
    output:
        file("${base}_taylor.fasta")
        file("${base}_bcftools.vcf")
        file("${base}_pre_bcftools.vcf")
        file(INDEX_FILE)
        file("${base}_summary.csv")

    publishDir params.OUTDIR, mode: 'copy'

    shell:
    '''
    R1=!{base}

    if (( !{bamsize} > 92 ))
    then
        # Create coverage file from bam for whole genome, then pipe anything that has less than 6 coverage to bed file,
        # to be masked later
        /usr/local/miniconda/bin/bedtools genomecov \\
            -bga \\
            -ibam !{BAMFILE} \\
            -g !{REFERENCE_FASTA} \\
            | awk '\$4 < 6' | /usr/local/miniconda/bin/bedtools merge > \${R1}.mask.bed
        # Get rid of anything outside of the genome we care about, to prevent some sgrnas from screwing with masking
        awk '{ if(\$3 > 200 && \$2 < 29742) {print}}' \${R1}.mask.bed > a.tmp && mv a.tmp \${R1}.mask.bed

        # Mask refseq fasta for low coverage areas based on bed file
        /usr/local/miniconda/bin/bedtools maskfasta \\
            -fi !{REFERENCE_FASTA} \\
            -bed \${R1}.mask.bed \\
            -fo ref.mask.fasta
        
        # If everything is below 6 coverage, just make an empty fasta
        if grep -q "(printf '\t')0$(printf '\t')29903" \${R1}.mask.bed; then
            printf '>!{base}\n' > \${R1}_taylor.fasta
            printf 'n%.0s' {1..29539} >> \${R1}_taylor.fasta
        else
            # Align to Wuhan refseq and unwrap fasta
            cat ref.mask.fasta \${R1}.consensus.fa > align_input.fasta
            /usr/local/miniconda/bin/mafft --auto --thread !{task.cpus} align_input.fasta > repositioned.fasta
            awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' repositioned.fasta > repositioned_unwrap.fasta
            
            # Trim ends and aligns masking of refseq to our consensus
            python3 !{TRIM_ENDS} \${R1}
        fi

        # Find percent ns, doesn't work, fix later in python script
        num_bases=$(grep -v ">" \${R1}_taylor.fasta | wc | awk '{print $3-$1}')
        num_ns=$(grep -v ">" \${R1}_taylor.fasta | awk -F"n" '{print NF-1}')
        percent_n=$(awk -v num_ns=$num_ns -v num_bases=$num_bases 'BEGIN { print ( num_ns * 100 / num_bases ) }')
        echo "num_bases=$num_bases"
        echo "num_ns=$num_ns"
        echo "percent_n=$percent_n"
        cp \${R1}.vcf.gz \${R1}_2.vcf.gz
        gunzip \${R1}_2.vcf.gz
        mv \${R1}_2.vcf \${R1}_bcftools.vcf
        #/usr/local/miniconda/bin/samtools view !{BAMFILE} -@ !{task.cpus} | awk -F: '$12 < 600' > \${R1}'.clipped.cleaned.bam'
    else
       echo "Empty bam detected. Generating empty consensus fasta file..."
       # Generate empty stats for empty bam
       printf '>!{base}\n' > \${R1}_taylor.fasta
       printf 'n%.0s' {1..29539} >> \${R1}_taylor.fasta
       percent_n=100
    fi
    
    cp \${R1}_summary3.csv \${R1}_summary.csv
    printf ",\$percent_n" >> \${R1}_summary.csv

    cat \${R1}_summary.csv | tr -d "[:blank:]" > a.tmp
    mv a.tmp \${R1}_summary.csv

    # Correctly calculates %ns and cleans up the summary file.
    if [[ !{bamsize} > 92 ]]
    then
        python3 !{FIX_COVERAGE} \${R1}
        mv \${R1}_summary_fixed.csv \${R1}_summary.csv
    fi

    [ -s \${R1}_taylor.fasta ] || echo "WARNING: \${R1} produced blank output. Manual review may be needed."

    '''
}

// Grab variants from vcf and output into standard format
process AnnotateVariants {
    errorStrategy 'retry'
    maxRetries 3

    container "quay.io/vpeddu/lava_image:latest"

    input:
        tuple val(base),val(bamsize),file("${base}_pre_bcftools.vcf") //from Vcf_ch
        file MAT_PEPTIDES
        file MAT_PEPTIDE_ADDITION
        file RIBOSOMAL_SLIPPAGE
        file RIBOSOMAL_START
        file PROTEINS
        file AT_REFGENE
        file AT_REFGENE_MRNA
        file CORRECT_AF_BCFTOOLS
        file FIX_COMPLEX_MUTATIONS
        
    output: 
        file("${base}_bcftools_variants.csv")
        file("*")
    
    publishDir params.OUTDIR, mode: 'copy', pattern:'*_bcftools_variants*'

    shell:
    '''
    #!/bin/bash
    ls -latr
    
    if (( !{bamsize} > 92))
    then
        # Fixes ploidy issues.
        #awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)gsub("1/1","0/1",$10)gsub("1/1","1/0",$11)}1\' !{base}_lofreq.vcf > !{base}_p.vcf
        awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)gsub("1/1","1/0",$10)gsub("1/1","1/0",$11)}1\' !{base}_pre_bcftools.vcf > !{base}_p.vcf

        # Converts VCF to .avinput for Annovar.
        file="!{base}""_p.vcf"
        #convert2annovar.pl -withfreq -format vcf4 -includeinfo !{base}_p.vcf > !{base}.avinput 
        convert2annovar.pl -withfreq -format vcf4 -includeinfo !{base}_p.vcf > !{base}.avinput 
        annotate_variation.pl -v -buildver AT -outfile !{base} !{base}.avinput .

        #awk -F":" '($26+0)>=1{print}' !{base}.exonic_variant_function > !{base}.txt
        cp !{base}.exonic_variant_function variants.txt
        #grep "SNV" !{base}.txt > a.tmp
        #grep "stop" !{base}.txt >> a.tmp
        #mv a.tmp variants.txt
    
        awk -v name=!{base} -F'[\t:,]' '{print name","$6" "substr($9,3)","$12","$44+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$41}' variants.txt > !{base}.csv

        grep -v "transcript" !{base}.csv > a.tmp && mv a.tmp !{base}.csv 
        grep -v "delins" !{base}.csv > final.csv
        # Sorts by beginning of mat peptide
        sort -k2 -t, -n mat_peptides.txt > a.tmp && mv a.tmp mat_peptides.txt
        # Adds mature peptide differences from protein start.
        python3 !{MAT_PEPTIDE_ADDITION}
        rm mat_peptides.txt
        python3 !{CORRECT_AF_BCFTOOLS} -name !{base}
        # Corrects for ribosomal slippage.
        python3 !{RIBOSOMAL_SLIPPAGE} filtered_variants.csv proteins.csv
        awk NF final.csv > a.tmp && mv a.tmp final.csv
        echo "SAMPLE,gene,AAPOS,AAREF,AASUB,TCOV,VCOV,AAFREQ,NTPOS,snpid,nsp,NSPPOS,NSPREF,NSPSUB" > !{base}_bcftools_variants.csv
        #sort -h -k2 -t, visualization.csv >> !{base}_bcftools_variants.csv
        cat visualization.csv >> !{base}_bcftools_variants.csv
        python3 !{FIX_COMPLEX_MUTATIONS} -name !{base} -file !{base}_bcftools_variants.csv -fasta !{AT_REFGENE_MRNA}

    else 
        echo "Bam is empty, skipping annotation."
        touch !{base}_bcftools_variants.csv
        touch !{base}_bcftools_variants_edited.csv
    fi

    '''
}

    // process annotateVariants_Lofreq {
    //     errorStrategy 'retry'
    //     maxRetries 3

    //     container "quay.io/vpeddu/lava_image:latest"

    //     input:
    //         tuple val(base),val(bamsize),file("${base}_lofreq.vcf") from Vcf_ch2
    //         file MAT_PEPTIDES
    //         file MAT_PEPTIDE_ADDITION
    //         file RIBOSOMAL_SLIPPAGE
    //         file RIBOSOMAL_START
    //         file PROTEINS
    //         file AT_REFGENE
    //         file AT_REFGENE_MRNA
    //         file CORRECT_AF
            
    //     output: 
    //         file("${base}_lofreq_variants.csv")
        
    //     publishDir params.OUTDIR, mode: 'copy'

    //     shell:
    //     '''
    //     #!/bin/bash
    //     ls -latr
        
    //     if (( !{bamsize} > 92))
    //     then
    //         # Fixes ploidy issues.
    //         #awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)gsub("1/1","0/1",$10)gsub("1/1","1/0",$11)}1\' !{base}_lofreq.vcf > !{base}_p.vcf
    //         #awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)gsub("1/1","1/0",$10)gsub("1/1","1/0",$11)}1\' !{base}_bcftools.vcf > !{base}_p.vcf
    //         cp !{base}_lofreq.vcf !{base}_p.vcf

    //         # Converts VCF to .avinput for Annovar.
    //         file="!{base}""_p.vcf"
    //         #convert2annovar.pl -withfreq -format vcf4 -includeinfo !{base}_p.vcf > !{base}.avinput 
    //         convert2annovar.pl -withfreq -format vcf4 -includeinfo !{base}_p.vcf > !{base}.avinput 
    //         annotate_variation.pl -v -buildver AT -outfile !{base} !{base}.avinput .

    //         #awk -F":" '($26+0)>=1{print}' !{base}.exonic_variant_function > !{base}.txt
    //         cp !{base}.exonic_variant_function !{base}.txt
    //         grep "SNV" !{base}.txt > a.tmp
    //         grep "stop" !{base}.txt >> a.tmp
    //         mv a.tmp variants.txt
        
    //         awk -v name=!{base} -F'[\t:,]' '{print name","$6" "substr($9,3)","$12","$44+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$41}' variants.txt > !{base}.csv

    //         grep -v "transcript" !{base}.csv > a.tmp && mv a.tmp !{base}.csv 
    //         grep -v "delins" !{base}.csv > final.csv
    //         # Sorts by beginning of mat peptide
    //         sort -k2 -t, -n mat_peptides.txt > a.tmp && mv a.tmp mat_peptides.txt
    //         # Adds mature peptide differences from protein start.
    //         python3 !{MAT_PEPTIDE_ADDITION}
    //         rm mat_peptides.txt
    //         # Corrects for ribosomal slippage.
    //         python3 !{RIBOSOMAL_SLIPPAGE} final.csv proteins.csv
    //         awk NF final.csv > a.tmp && mv a.tmp final.csv
    //         python3 !{CORRECT_AF}
    //         sort -h -k2 -t, fixed_variants.txt > !{base}_lofreq_variants.csv
    //     else 
    //         echo "Bam is empty, skipping annotation."
    //         touch !{base}_lofreq_variants.csv
    //     fi

    //     '''
    // }

//     process varscan2 { 
//         container "quay.io/vpeddu/lava_image:latest"

//         // Retry on fail at most three times 
//         errorStrategy 'retry'
//         maxRetries 3

//         input:
//             tuple val (base), file(BAMFILE), file(INDEX_FILE),val(bamsize) from Clipped_bam_ch3
//             file REFERENCE_FASTA
//             file REFERENCE_FASTA_FAI
//             file SPLITCHR
//         output:
//             tuple val(base),val(bamsize),file("${base}_varscan.vcf") into Varscan_ch

//         publishDir params.OUTDIR, mode: 'copy'

//         shell:
//         '''
//         #!/bin/bash
//         ls -latr
//         R1=`basename !{BAMFILE} .clipped.bam`
//         echo "bamsize: !{bamsize}"
//         #if [ -s !{BAMFILE} ]
//         # More reliable way of checking bam size, because of aliases
//         if (( !{bamsize} > 92 ))
//         then
//             # Parallelize pileup based on number of cores
//             splitnum=$(($((29903/!{task.cpus}))+1))
//             cat !{SPLITCHR} | \\
//                 xargs -I {} -n 1 -P !{task.cpus} sh -c \\
//                     "/usr/local/miniconda/bin/samtools mpileup \\
//                         -f !{REFERENCE_FASTA} -r {} \\
//                         -B \\
//                         --max-depth 50000 \\
//                         --max-idepth 500000 \\
//                     !{BAMFILE} | 
//                     java -jar /usr/local/bin/VarScan mpileup2cns --validation 1 --output-vcf 1 --min-coverage 2 --min-var-freq 0.001 --p-value 0.99 --min-reads2 1 > tmp.{}.vcf"
            
//             cat *.vcf > \${R1}_varscan.vcf
//         else
//         touch \${R1}_varscan.vcf
//         fi
//         '''
// }
// }
