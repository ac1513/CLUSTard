samples= 'NAB1_T0 NAB1_T1 NAB1_T2 NAB1_T3' #should be in order want output
JOBID = 'test'
RAW_SR = 'data/'
REFIN = 'data/yw_polished_anvio.fasta'
THRESH = '10000'
P_THRESH = '0.99'

#subworkflow bwa_split:
#    snakefile:
#        "bwa_Snakefile"

subworkflow para:
    snakefile:
        "para_Snakefile"

rule all:
    input:
#        "test.txt",
        expand("results/{JOBID}_summary_stats.txt", JOBID = JOBID),
        expand("plots/{JOBID}_all_plot.pdf", JOBID = JOBID)

#rule test:
#    input:
#        bwa_split(expand("inter/{JOBID}_output.txt", JOBID = JOBID))
#    output:
#        "test.txt"
#    shell:
#        """
#        echo "Done BWA" > {output}
#        """

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

rule plot:
    input:
         expand("results/{JOBID}_summary_stats.txt", JOBID = JOBID)
    output:
        "plots/{JOBID}_all_plot.pdf"
    params:
        files = "plot_in_files.txt",
        csv_loc = "results/"
    conda:
        "envs/plot2.yaml"
    shell:
        """
        ls -S results/Cluster*.csv > {params.files}
        python scripts/plot.py {params.files} {params.csv_loc}
        rm {params.files}
        """
