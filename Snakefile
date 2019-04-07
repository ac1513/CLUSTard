samples= 'NAB1_T0 NAB1_T1 NAB1_T2 NAB1_T3' #should be in order want output
JOBID = 'test'
RAW_SR = 'data/'
REFIN = 'data/yw_polished_anvio.fasta'
THRESH = '10000'
P_THRESH = '0.99'
PARTS = "00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 "

subworkflow bwa_split:
    snakefile:
        "bwa_Snakefile"

subworkflow para:
    snakefile:
        "para_Snakefile"

rule all:
    input:
        "test.txt",
        "test1.txt"

rule test:
    input:
        bwa_split(expand("inter/{JOBID}_output.txt", JOBID = JOBID))
    output:
        "test.txt"
    shell:
        """
        echo "Done" > {output}
        """

rule test2:
    input:
        para(expand("bins/{JOBID}_parallel_merged.out", JOBID = JOBID))
    output:
        "test1.txt"
    shell:
        """
        echo "Done" > {output}
        """
