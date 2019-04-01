JOBID = "Methanoregulaceae_4"

subworkflow create_bin:
    snakefile:
        "bin_Snakefile"

rule all:
    input:
        expand("{jobID}_checkM", jobID=JOBID),
        expand("trees/{jobID}_clustal_concat_msa.fasta", jobID=JOBID),
        expand("trees/{jobID}_clustal_concat_msa.treefile", jobID=JOBID)

rule checkm:
    input: create_bin(expand("plots/1_{jobID}_plot.pdf", jobID=JOBID))
    output:
        expand("{jobID}_checkM", jobID=JOBID)
    log:
        expand("logs/tree/{jobID}_checkm.log", jobID=JOBID)
    threads:
        40
    conda:
        "envs/checkm.yaml"
    priority: 500
    shell:
        """
        module load bio/CheckM
        checkm lineage_wf -x fasta -t {threads} genomes/ {output}
        module unload lang/Python/2.7.15-foss-2018b
        """

rule clustalo:
    input:
        expand("{jobID}_checkM", jobID=JOBID)
    output:
        expand("trees/{jobID}_clustal_concat_msa.fasta", jobID=JOBID)
    log:
        expand("logs/tree/{jobID}_clustalo.log", jobID=JOBID)
    priority: 400
    conda:
        "envs/msa_tree.yaml"
    shell:
        """
        clustalo -i {input}/storage/tree/concatenated.fasta --dealign -o {output}
        """

rule iq_tree:
    input:
        expand("trees/{jobID}_clustal_concat_msa.fasta", jobID=JOBID)
    output:
        expand("trees/{jobID}_clustal_concat_msa.treefile", jobID=JOBID)
    threads:
        40
    priority: 300
    conda:
        "envs/msa_tree.yaml"
    params:
        pre = expand("trees/{jobID}_clustal_concat_msa", jobID=JOBID)
    shell:
        """
        iqtree -s {input} -pre {params.pre} -nt {threads}
        """
