#!/bin/bash 
#SBATCH --account=BIOL-CHONG-2019
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpjc1@york.ac.uk
#SBATCH --output=merger_%j.log
#SBATCH --mem=64G
#SBATCH --job-name=merger
#SBATCH --workdir=/mnt/lustre/groups/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/

module load lang/Python/3.6.6-foss-2018b
python parallel_merge_step2.py -i bins/parallel_sets -o bins/parallel_merged.out -n 7