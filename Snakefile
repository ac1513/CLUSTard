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

rule all:
    input:
        expand("plots/1_{JOBID}_plot.pdf", JOBID = JOBID)

include: "bwa_Snakefile"

include: "para_Snakefile"

include: "kraken2_Snakefile"

rule plot:
    input:
         file_out = expand("results/{JOBID}_summary_stats.txt", JOBID = JOBID),
         kraken = (expand("kraken/{JOBID}_{kraken_level}_top_kraken.out", JOBID = JOBID, kraken_level = kraken_level))
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
