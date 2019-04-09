samples= 'NAB1_T0 NAB1_T1 NAB1_T2 NAB1_T3 NAB1_T4 NAB1_T5 NAB1_T6 NAB1_T7 NAB1_T8 NAB1_T9 NAB1_T10 NAB1_T11 NAB1_T12 NAB1_T13 NAB1_T14 NAB1_T15 NAB1_T16 NAB1_T17 NAB2_T0 NAB2_T1 NAB2_T2 NAB2_T3 NAB2_T4 NAB2_T5 NAB2_T6 NAB2_T7 NAB2_T8 NAB2_T9 NAB2_T10 NAB2_T11 NAB2_T12 NAB2_T13 NAB2_T14 NAB2_T15 NAB2_T16 NAB2_T17 NAB3_T0 NAB3_T1 NAB3_T2 NAB3_T3 NAB3_T4 NAB3_T5 NAB3_T6 NAB3_T7 NAB3_T8 NAB3_T9 NAB3_T10 NAB3_T11 NAB3_T12 NAB3_T13 NAB3_T14 NAB3_T15 NAB3_T16 NAB3_T17 NAB4_T0 NAB4_T1 NAB4_T2 NAB4_T3 NAB4_T4 NAB4_T5 NAB4_T6 NAB4_T7 NAB4_T8 NAB4_T9 NAB4_T10 NAB4_T11 NAB4_T12 NAB4_T13 NAB4_T14 NAB4_T15 NAB4_T16 NAB4_T17 NAB_Feed_T9 NAB_Feed_T10 NAB_Feed_T11 NAB_Feed_T12 NAB_Feed_T13 NAB_Feed_T14 NAB_Feed_T15 NAB_Feed_T17' #should be in order want output
JOBID = 'NAB_all'
RAW_SR = 'data/'
REFIN = 'data/yw_polished_anvio.fasta'
THRESH = '1000'
P_THRESH = '0.99'
krakendb = "/mnt/lustre/groups/biol-chong-2019/databases/krakendb/kraken2_samstudio8/"

subworkflow bwa_split:
    snakefile:
        "bwa_Snakefile"

subworkflow para:
    snakefile:
        "para_Snakefile"

rule all:
    input:
        "test.txt",
        expand("results/{JOBID}_summary_stats.txt", JOBID = JOBID),
        expand("kraken/{JOBID}_top_kraken.out", JOBID = JOBID),
        expand("plots/1_{JOBID}_plot.pdf", JOBID = JOBID)

rule test:
    input:
        bwa_split(expand("inter/{JOBID}_output.txt", JOBID = JOBID))
    output:
        "test.txt"
    shell:
        """
        echo "Done BWA" > {output}
        """

rule file_parser:
    input:
        clusters = para(expand("bins/{JOBID}_non_red_list.out", JOBID = JOBID))
    output:
        results = "results/{JOBID}_summary_stats.txt"
    params:
        contigs = REFIN,
        csv = expand("inter/{JOBID}_read_counts_derived.csv", JOBID = JOBID),
        wd = "results/",
        header = samples
    shell:
        """
        python scripts/file_parser.py {params.contigs} {params.csv} {input.clusters} {params.wd} {output.results} -l {params.header}
        echo "Sort output at somepoint" >> {output.results}
        """

clusters = glob_wildcards("results/{clusters}.fasta")
rule kraken:
    input:
        expand("results/{clusters}.fasta", clusters=clusters)
    output:
        "kraken/{clusters}_kraken.out"
    threads:
        4
    params:
        db = krakendb
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        kraken2 -db {params.db} --threads 4 --report {output} --output ${NAME}_kraken_names.out --use-names $FILE
        """

rule kraken_merge:
    input: expand("kraken/{clusters}_kraken.out", clusters = clusters)
    output: expand("kraken/{JOBID}_top_kraken.out", JOBID = JOBID)
    params:
        level = 'F'
    shell:
    """
    find -name '*report.out' -type f -printf '\n%p\t' -exec sh -c 'echo {{}} | sort -k1nr {{}} | grep -P "\t{level}\t" | head -n1 ' \;
    """

rule plot:
    input:
         expand("results/{JOBID}_summary_stats.txt", JOBID = JOBID),
         expand("kraken/{JOBID}_top_kraken.out", JOBID = JOBID)
    output:
        "plots/1_{JOBID}_plot.pdf"
    params:
        files = "plot_in_files.txt"
    conda:
        "envs/plot2.yaml"
    shell:
        """
        ls -S results/Cluster*.fasta > {params.files}
        sed -i "s/.fasta/.csv/g" {params.files}
        python scripts/plot.py {params.files} {JOBID}
        rm {params.files}
        """
