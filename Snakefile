samples= 'NAB1_T0 NAB1_T1 NAB1_T2 NAB1_T3' #should be in order want output
JOBID = 'test'
RAW_SR = 'data/'
REFIN = 'data/yw_polished_anvio.fasta'

#run using snakemake --cluster "sbatch -t 02:00:00" -j 20

rule all:
    input:
        expand("{REFIN}.sa", REFIN=REFIN),
        expand('bwa_out/{samples}.bam', samples=samples.split(' ')),
        expand('inter/counts_{samples}.txt', samples=samples.split(' ')),
        expand('inter/{jobid}_read_counts.out', jobid= JOBID)

rule bwa_index:
    input:
        ref = REFIN
    output:
        '{REFIN}.sa'
    shell:
        'bwa index {input.ref}'

rule bwa_mem:
    input:
        fq1 = RAW_SR + '{samples}_R1.fastq.gz',
        fq2 = RAW_SR + '{samples}_R2.fastq.gz',
        ref = REFIN,
        ref_ind = expand("{reference}.sa", reference=REFIN)
    output:
        bam = 'bwa_out/{samples}.bam'
    threads: 16
    shell:
        """
        bwa mem -M -t {threads} {input.ref} {input.fq1} {input.fq2} | samtools view -bhS - | samtools sort -o {output.bam}
        """

rule samtools:
    input:
        bam = expand('bwa_out/{samples}.bam',samples=samples.split(' '))
    output:
        counts = expand('inter/counts_{samples}.txt', samples=samples.split(' ')),
    shell:
        """
        samtools index {input.bam}
        samtools idxstats {input.bam} > {output.counts}
        """

rule merge_filecounts:
    input:
        loc = "inter/",
        counts = expand('inter/counts_{samples}.txt', samples=samples.split(' ')),
    output:
        txt = 'inter/{JOBID}_read_counts.out'
    conda:
        "envs/py3"
    shell:
        """
        python scripts/merge_filecounts.py inter {JOBID}
        rm inter/*.txt
        """
