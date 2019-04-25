configfile: "config.yaml"

import pandas as pd
df_samples = pd.read_csv(config["samples"], sep ='\t', index_col = 0)
samples = df_samples["sample"].to_list()

JOBID = config["jobid"]
RAW_SR = config["RAW_SR"]
REFIN = config["REFIN"]
CONTIG_T = config["CONTIG_T"]
P_THRESH = config["P_THRESH"]
krakendb = config["krakendb"]
kraken_level = config["kraken_level"]

subworkflow bwa_split:
    snakefile:
        "bwa_Snakefile"

subworkflow para:
    snakefile:
        "para_Snakefile"

subworkflow kraken2:
    snakefile:
        "kraken2_Snakefile"


rule all:
    input:
        "test.txt",
        expand("results/{JOBID}_summary_stats.txt", JOBID = JOBID),
        expand("plots/1_{JOBID}_plot.pdf", JOBID = JOBID)

rule test:
    input:
        bwa_split(expand("inter/{JOBID}_bwa_output.txt", JOBID = JOBID))
    output:
        "test.txt"
    shell:
        """
        echo "Done BWA" > {output}
        """

rule para_out:
    input:
        clusters = para(expand("results/{JOBID}_summary_stats.txt", JOBID = JOBID))
    output:
        "results/{JOBID}_summary_stats.txt"
    shell:
        """
        echo "Done" >> {output}
        """

rule plot:
    input:
         file_out = expand("results/{JOBID}_summary_stats.txt", JOBID = JOBID),
         kraken = kraken2(expand("kraken/{JOBID}_{kraken_level}_top_kraken.out", JOBID = JOBID, kraken_level = kraken_level))
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
        python scripts/plot.py {params.files} {JOBID} -k {input.kraken}
        rm {params.files}
        """
