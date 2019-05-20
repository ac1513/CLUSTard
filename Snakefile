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
        "scripts/bwa_Snakefile"

subworkflow para:
    snakefile:
        "scripts/para_Snakefile"


subworkflow kraken2:
    snakefile:
        "scripts/kraken2_Snakefile"


rule all:
    input:
        expand("logs/{JOBID}_all_bwa_output.txt", JOBID=JOBID),
        expand("logs/{JOBID}_para_out.txt", JOBID = JOBID),
        expand("output/plots/1_{JOBID}_{kraken_level}_plot.pdf", JOBID = JOBID, kraken_level = kraken_level)

localrules: plot

rule test:
    input:
        bwa_split(expand("output/clustering/{JOBID}_bwa_output.txt", JOBID = JOBID))
    output:
        "logs/{JOBID}_all_bwa_output.txt"
    shell:
        """
        more *.out > {output} 2> /dev/null
        rm *.out
        """

rule para_out:
    input:
        clusters = para(expand("logs/{JOBID}_para_done.txt", JOBID = JOBID))
    output:
        "logs/{JOBID}_para_out.txt"
    shell:
        """
        echo "Done" >> {output}
        """

rule plot:
    input:
         file_out = expand("logs/{JOBID}_para_out.txt", JOBID = JOBID),
         wait = kraken2("done.txt")
    output:
        "output/plots/1_{JOBID}_{kraken_level}_plot.pdf"
    params:
        files = "plot_in_files.txt",
        sample_file = config["samples"],
        kraken = expand("output/kraken/{JOBID}_{kraken_level}_top_kraken.out", JOBID = JOBID, kraken_level = kraken_level)
    conda:
        "envs/py3.yaml"
    shell:
        """
        ls -S output/results/Cluster*.fasta > {params.files}
        sed -i "s/.fasta/.csv/g" {params.files}
        python scripts/plot.py {params.files} {JOBID} -k {params.kraken} -k_l {kraken_level} {params.sample_file}
        rm {params.files}
        """
