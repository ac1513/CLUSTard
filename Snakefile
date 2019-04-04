samples= 'NAB1_T0 NAB1_T1 NAB1_T2 NAB1_T3' #should be in order want output
JOBID = 'test'
RAW_SR = 'data/'
REFIN = 'data/yw_polished_anvio.fasta'
BWAOUT = 'bwa_out/'
INTER = 'inter_files/'

#run using snakemake --cluster "sbatch -t 02:00:00" -j 20

rule all:
    input:
        expand(BWAOUT + '{samples}.bam', samples=samples.split(' ')),
        expand(INTER + JOBID + '_read_counts.out')

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
        bam = BWAOUT + '{samples}.bam',
        counts = INTER + 'counts_{samples}.txt'
    threads: 16
    shell:
        """
        bwa mem -M -t {threads} {input.ref} {input.fq1} {input.fq2} | samtools view -bhS - | samtools sort -o {output.bam}
        sleep 2s
        samtools index {output.bam}
        samtools idxstats {output.bam} > {output.counts}
        """

rule merge_filecounts:
    input:
        loc = INTER
    output:
        txt = INTER + '{JOBID}_read_counts.out'
    conda:
        "envs/py3"
    shell:
        """
        python scripts/merge_filecounts.py {INTER} {JOBID}
        rm {INTER}.*txt
        """
