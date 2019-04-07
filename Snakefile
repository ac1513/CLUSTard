samples= 'NAB1_T0 NAB1_T1 NAB1_T2 NAB1_T3' #should be in order want output
JOBID = 'test'
RAW_SR = 'data/'
REFIN = 'data/yw_polished_anvio.fasta'
THRESH = '10000'
P_THRESH = '0.99'

#run using snakemake --cluster "sbatch -t 02:00:00" -j 20

rule all:
    input:
        expand("{REFIN}.sa", REFIN=REFIN),
        expand('inter/counts_{samples}.txt', samples=samples.split(' ')),
        expand('inter/{jobid}_read_counts.out', jobid= JOBID),
        expand('inter/{jobid}_read_counts_derived.csv', jobid= JOBID),
        expand('inter/{jobid}_values.csv', jobid = JOBID),
        expand('inter/{jobid}_diffs.csv', jobid = JOBID),
        dynamic(expand("bins/{JOBID}_diffs{{PART}}.csv", JOBID=JOBID)),
        expand("bins/{JOBID}_parallel_merged.out", JOBID = JOBID)

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
        diffs = expand("inter/{JOBID}_diffs.csv", JOBID=JOBID)
    output:
        dynamic(expand('inter/{JOBID}_diffs{{PART}}.csv', JOBID = JOBID))
    params:
        diffs = expand("inter/{JOBID}_diffs", JOBID = JOBID)
    shell:
        """
        split -d -l 10000 --additional-suffix=.csv {input.diffs} {params.diffs}
        """

(job, part) = glob_wildcards('inter/{JOBID}_diffs{PART}.csv')

rule bin_feeder:
    input:
        #diffs = expand('inter/{JOBID}_diffs{PART}.csv', JOBID = JOBID, PART=part),
        dyn_diffs = dynamic(expand("bins/{JOBID}_diffs{{PART}}.csv", JOBID=JOBID)),
        diffs = expand('inter/{JOBID}_diffs{PART}.csv', JOBID = JOBID, PART = part)
    output:
        all = expand("bins/{JOBID}_output_{PART}.csv", JOBID = JOBID, PART = part),
    params:
        thresh = P_THRESH, #add this in as a variable at the top later..
        all_diffs = expand("inter/{JOBID}_diffs.csv", JOBID = JOBID)
    conda:
        "envs/py3.yaml"
    shell:
        """
        python scripts/bin_feeder.py {input.diffs} {params.all_diffs} {params.thresh} {output.all}
        """

rule para_sets:
    input:
        bins = expand("bins/{JOBID}_output_{PART}.csv", JOBID = JOBID, PART = part)
    output:
        expand("bins/{JOBID}_parallel_sets_{PART}.csv", JOBID = JOBID, PART = part)
    params:
        thresh = P_THRESH
    conda:
        "envs/py3.yaml"
    shell:
        """
        python scripts/para_sets.py {input.bins} {output} {params.thresh}
        """

rule para_merge:
    input:
        expand("bins/{JOBID}_parallel_sets_{PART}.csv", JOBID=JOBID, PART = part)
    output:
        "bins/{JOBID}_parallel_merged.out"
    conda:
        "envs/py3.yaml"
    shell:
        """
        python scripts/parallel_merge_step2.py -i {input} -o {output}
        """
