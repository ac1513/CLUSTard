# CLUSTard

## How to run
1. Install snakemake into conda environment
2. git clone https://github.com/ac1513/CLUSTard.git
3. Create data directory (containing samples [names ending in \_R1.fastq.gz or \_R2.fastq.gz] and assembled reference *#need to change this to be more universal*)
4. Put sample prefixes into samples.tsv  (samples column) in the order that you want the final plot to be
5. Edit config.yaml file
6. Run using ``snakemake --use-conda --cluster "sbatch -t 48:00:00 --cpus-per-task={threads}" -j 6000`` (running in screen with logfile -L is best!)
* (Add --config then e.g. jobid=test to change config without having to edit the whole file - esp useful with the jobid and the kraken_level steps.)
* Currently needs to be run with Snakemake version 5.6.0 or lower
