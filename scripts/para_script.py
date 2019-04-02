# python script to generate parallel_step2 shell script
# directory for sequence files is assumed to be ..path../infiles

import os
import glob

path = os.path.dirname(os.path.realpath(__file__)) # path goes here
dir = '/' # directory containing sequence files goes here
csv_file = 'read_counts_derived.csv'

num_lines = sum(1 for line in open(path+dir+csv_file))

num_instances = num_lines//10000

f = open(path+'/parallel_step2_90.sh', 'w')

f.write('#!/bin/bash \n#SBATCH --account=BIOL-CHONG-2019\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=jpjc1@york.ac.uk\n')
f.write('#SBATCH --output=para_%j.log\n')
f.write('#SBATCH --array=1-'+str(num_instances)+'\n#SBATCH --job-name=para\n#SBATCH --workdir='+path+dir+'\n\n')
f.write('dir=$(cd -P -- \"$(dirname -- \"$0\")\" && pwd -P)\n\n')
f.write('module load lang/Python/3.6.6-foss-2018b\n')
f.write('python new_step2_90_parallel.py -i bins2/output.$SLURM_ARRAY_TASK_ID -o bins2/parallel_sets.$SLURM_ARRAY_TASK_ID')
