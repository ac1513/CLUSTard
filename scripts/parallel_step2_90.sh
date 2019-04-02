#!/bin/bash 
#SBATCH --account=BIOL-CHONG-2019
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpjc1@york.ac.uk
#SBATCH --output=para_%j.log
#SBATCH --array=1-7
#SBATCH --job-name=para
#SBATCH --workdir=/mnt/lustre/groups/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

module load lang/Python/3.6.6-foss-2018b
python new_step2_90_parallel.py -i bins2/output.$SLURM_ARRAY_TASK_ID -o bins2/parallel_sets.$SLURM_ARRAY_TASK_ID