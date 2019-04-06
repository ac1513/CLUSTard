samples= 'NAB1_T0 NAB1_T1 NAB1_T2 NAB1_T3' #should be in order want output
JOBID = 'test'
RAW_SR = 'data/'
REFIN = 'data/yw_polished_anvio.fasta'
THRESH = '10000'

#run using snakemake --cluster "sbatch -t 02:00:00" -j 20

rule all:
    input:
        expand("{REFIN}.sa", REFIN=REFIN),
        expand('inter/counts_{samples}.txt', samples=samples.split(' ')),
        expand('inter/{jobid}_read_counts.out', jobid= JOBID),
        expand('inter/{jobid}_read_counts_derived.csv', jobid= JOBID),
        dynamic('inter/{JOBID}_read_counts_derived{PART}.csv'),
        expand('inter/{jobid}_values.csv', jobid = JOBID),
        expand('inter/{jobid}_diffs.csv', jobid = JOBID),
        "inter/test.txt"

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
        counts = 'inter/counts_{samples}.txt'
    params:
        bam = 'bwa_out/{samples}.bam'
    threads: 16
    shell:
        """
        mkdir -p bwa_out
        bwa mem -M -t {threads} {input.ref} {input.fq1} {input.fq2} | samtools view -bhS - | samtools sort -o {params.bam}
        samtools index {params.bam}
        samtools idxstats {params.bam} > {output.counts}
        """

rule merge_filecounts:
    input:
        'inter/'
    output:
        txt = 'inter/{JOBID}_read_counts.out'
    conda:
        "envs/py3.yaml"
    shell:
        """
        python scripts/merge_filecounts.py inter {JOBID}
        """

rule derive:
    input:
        expand('inter/{JOBID}_read_counts.out', JOBID=JOBID)
    output:
        csv = "inter/{JOBID}_read_counts_derived.csv"
    params:
        thresh = THRESH
    conda:
        "envs/py3.yaml"
    shell:
        """
        python scripts/derive.py inter {JOBID} {params.thresh}
        """

rule start_feeder:
    input:
        expand('inter/{JOBID}_read_counts_derived.csv', JOBID=JOBID)
    output:
        values = "inter/{JOBID}_values.csv",
        diffs = "inter/{JOBID}_diffs.csv"
    conda:
        "envs/py3.yaml"
    shell:
        """
        python scripts/start_feeder.py inter {JOBID}
        """

rule split_file:
    input:
        values = expand('inter/{JOBID}_values.csv', JOBID=JOBID),
        diffs = expand("inter/{JOBID}_diffs.csv", JOBID=JOBID)
    output:
        dynamic(expand('inter/{JOBID}_diffs{{PART}}.csv', JOBID = JOBID))
    params: out = expand("inter/{JOBID}_read_counts_derived", JOBID = JOBID)
    shell:
        """
        split -d -l 10000 --additional-suffix=.csv {input.values} {params.out}
        split -d -l 10000 --additional-suffix=.csv {input.diffs} {params.out}
        """
(job, part) = glob_wildcards('inter/{JOBID}_diffs{PART}.csv')

rule bin_feeder:
    input:
        diffs = dynamic(expand('inter/{JOBID}_diffs{{PART}}.csv', JOBID = JOBID),
        values = expand('inter/{jobid}_values{PART}.csv', jobid = JOBID, PART=part),
        all_diffs = expand("inter/{JOBID}_diffs.csv", JOBID = JOBID),
        all_values = expand("inter/{JOBID}_values.csv", JOBID = JOBID)
    output:
        "inter/test.txt"
    shell:
        """
        echo "{input.diffs} {input.values}" >> {output}
        """
