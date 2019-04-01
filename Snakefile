JOBID = "test1"

subworkflow create_bin:
    snakefile:
        "bin_Snakefile"

rule all:
    input:
        expand("{jobID}_checkM", jobID=JOBID),
        expand("trees/{jobID}_clustal_concat_msa.fasta", jobID=JOBID),
        expand("trees/{jobID}_clustal_concat_msa.treefile", jobID=JOBID)

rule checkm:
    input: create_bin("logs/kraken_out.log")
    output:
        expand("{jobID}_checkM", jobID=JOBID)
    log:
        expand("logs/tree/{jobID}_checkm.log", jobID=JOBID)
    threads:
        40
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
    params:
        pre = expand("trees/{jobID}_clustal_concat_msa", jobID=JOBID)
    shell:
        """
        module load bio/IQ-TREE
        iqtree-mpi -s {input} -pre {params.pre} -nt {threads}
        """
