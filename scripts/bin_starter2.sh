#!/bin/bash 
#SBATCH --account=BIOL-CHONG-2019
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpjc1@york.ac.uk
#SBATCH --output=binfeed_%j.log
#SBATCH --array=1-7
#SBATCH --job-name=b_f_%j
#SBATCH --workdir=/mnt/lustre/groups/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/

mkdir bins2

dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

STAR=$((SLURM_ARRAY_TASK_ID-1))

START=$((STAR*10000))

END=$((START+9999))

if (($SLURM_ARRAY_TASK_ID == 7)); then
    END=$(wc -l < read_counts_derived.csv)
fi

module load lang/Python/3.6.6-foss-2018b
python bin_feeder2_90.py -i read_counts_derived.csv -o bins2/output.$SLURM_ARRAY_TASK_ID -s $START -e $END