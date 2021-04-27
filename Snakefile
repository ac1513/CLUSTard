configfile: "config.yaml"

import pandas as pd
df_samples = pd.read_csv(config["samples"], sep ='\t', index_col = 0)
samples = df_samples["sample"].to_list()

JOBID = config["jobid"]
RAW_SR = config["RAW_DIR"]
REFIN = config["REFIN"]
CONTIG_T = config["CONTIG_T"]
P_THRESH = config["P_THRESH"]
krakendb = config["krakendb"]
kraken_level = config["kraken_level"]
plot_order = config["plot_order"]
#for plotting
date_scale = config["date_scale"]
rel_or_abs = "r"
top20 = "y"

if 'y' in top20:
    out_abun = rel_or_abs + '_top20'
else:
    out_abun = rel_or_abs

subworkflow bwa_split:
    snakefile:
        "scripts/bwa_Snakefile"

subworkflow para:
    snakefile:
        "scripts/para_Snakefile"

subworkflow singleton:
    snakefile:
        "scripts/singleton_Snakefile"

subworkflow kraken2:
    snakefile:
        "scripts/kraken2_Snakefile"


rule all:
    input:
        expand("logs/{JOBID}_all_bwa_output.txt", JOBID=JOBID),
        expand("logs/{JOBID}_para_out.txt", JOBID = JOBID),
        expand("logs/{JOBID}_singleton_out.txt", JOBID=JOBID),
        expand("output/plots/1_{JOBID}_{kraken_level}_plot.png", JOBID = JOBID, kraken_level = kraken_level),
        expand("output/plots/{JOBID}_bin_contigs.png", JOBID = JOBID),
        expand("output/clustering/{JOBID}_read_counts_absolute.csv", JOBID = JOBID),
        expand("output/plots/{JOBID}_{out_abun}_abun_plot.png", JOBID = JOBID, out_abun = out_abun),
        expand("output/{JOBID}_cluster_summary_stats.tsv", JOBID=JOBID),
        expand("output/{JOBID}_qual_MAGs.txt", JOBID=JOBID)


localrules: test, para_out, singleton_out, plot, bin_plot, abs_derive, abun_plot, clus_stats, high_mags

rule test:
    input:
        bwa_split(expand("output/clustering/{JOBID}_bwa_output.txt", JOBID = JOBID))
    output:
        touch("logs/{JOBID}_all_bwa_output.txt")
    shell:
        """
        echo "Done with BWA"
        """

rule para_out:
    input:
        clusters = para(expand("logs/{JOBID}_singleton_step1.txt", JOBID=JOBID))
    output:
        touch("logs/{JOBID}_para_out.txt")
    shell:
        """
        echo "Clustering completed"
        """

rule singleton_out:
    input:
        singleton(expand("logs/{JOBID}_singleton_done.txt", JOBID=JOBID))
    output:
        touch("logs/{JOBID}_singleton_out.txt")
    shell:
        """
        echo "singleton analysis done"
        """

rule clus_stats:
    input:
        file_out = expand("logs/{JOBID}_singleton_out.txt", JOBID = JOBID),
        checkm = kraken2(expand("output/{JOBID}_checkm/{JOBID}_checkm.log", JOBID=JOBID)),
    output:
        csv = "output/{JOBID}_cluster_summary_stats.tsv"
    conda:
        "envs/py3.yaml" #change clustering (below) when add counts folder..
    params:
        checkm = expand("output/{JOBID}_checkm/{JOBID}_checkm.log", JOBID=JOBID),
        seqk = expand("output/results/{JOBID}_seqkit_stats.tsv", JOBID=JOBID),
        kraken = expand("output/kraken/{JOBID}_{kraken_level}_top_kraken.out", JOBID = JOBID, kraken_level = kraken_level),
        gtdb_ncbi = expand("output/kraken/{JOBID}_{kraken_level}_GTDB_lookup.json", JOBID = JOBID, kraken_level = kraken_level)
    shell:
        """
        ls output/results/C*.csv > stat_input.txt
        python scripts/clus_stats.py stat_input.txt {JOBID} -cm {params.checkm} -sk {params.seqk} -k {params.kraken} -g {params.gtdb_ncbi}
        rm stat_input.txt
        """

rule plot:
    input:
        stats_in = expand("output/{JOBID}_cluster_summary_stats.tsv", JOBID = JOBID),
    output:
        cluster_plot = "output/plots/1_{JOBID}_{kraken_level}_plot.png"
    params:
        date = date_scale,
        sample_list = config["samples"],
        out_order = plot_order
    conda:
        "envs/py3.yaml"
    shell:
        """
        python scripts/plot.py {input.stats_in} {JOBID} {params.sample_list} {params.date} -k_l {kraken_level} -s {params.out_order}
        """

rule bin_plot:
    input:
        file_out = expand("output/plots/1_{JOBID}_{kraken_level}_plot.png", JOBID = JOBID, kraken_level = kraken_level)
    output:
        contig_plot = "output/plots/{JOBID}_bin_contigs.png"
    conda:
        "envs/py3.yaml"
    shell:
        """
        cd output/results/
        cat Cluster*.fasta | awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' | sort | uniq > {JOBID}_sorted_lengths.tsv
        cd ../../
        python scripts/bin_plot.py output/results/{JOBID}_unbinned_contigs_stats.tsv output/results/{JOBID}_sorted_lengths.tsv {JOBID}
        """

rule abs_derive:
    input:
        plots = expand("output/plots/{JOBID}_bin_contigs.png", JOBID=JOBID)
    output:
        csv = "output/clustering/{JOBID}_read_counts_absolute.csv"
    params:
        thresh = CONTIG_T,
        counts = expand('output/clustering/{JOBID}_merged_counts.tsv', JOBID=JOBID)
    conda:
        "envs/py3.yaml" #change clustering (below) when add counts folder..
    shell:
        """
        python scripts/absolute_derive.py {params.counts} {JOBID} {params.thresh}
        """

rule abun_plot:
    input:
        count_in = expand("output/clustering/{JOBID}_read_counts_absolute.csv", JOBID = JOBID)
    output:
        plot_out = "output/plots/{JOBID}_{out_abun}_abun_plot.png"
    conda:
        "envs/py3.yaml"
    params:
        roa = rel_or_abs,
        top_20 = top20,
        kraken_in = expand("output/kraken/{JOBID}_{kraken_level}_top_kraken.out", JOBID = JOBID, kraken_level = kraken_level)
    shell:
        """
        cd output/results/
        for f in C*.fasta; do filename="${{f%%.*}}"; echo ">$f"; seqkit fx2tab -n $f; done > {JOBID}_binned_cluster_contig.txt
        cd ../../
        python scripts/abun_plot.py {JOBID} {input.count_in} output/results/{JOBID}_binned_cluster_contig.txt {params.roa} {params.top_20} -s {samples} Coverage -k {params.kraken_in}
        """

rule high_mags:
    input:
        expand("output/{JOBID}_cluster_summary_stats.tsv", JOBID = JOBID)
    output:
        txt = "output/{JOBID}_qual_MAGs.txt"
    params:
        checkm = expand("output/{JOBID}_checkm/{JOBID}_checkm.log", JOBID = JOBID),
        prokka = expand("output/{JOBID}_prokka/", JOBID = JOBID)
    shell:
        """
        python scripts/qual_parse.py {params.checkm} {params.prokka} > {output.txt}
        rm -r tmp/
        """
