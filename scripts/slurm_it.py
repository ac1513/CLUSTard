# python script to generate bin_starter shell script
# directory for sequence files is assumed to be ..path../infiles

import os
import glob

path = os.path.dirname(os.path.realpath(__file__)) # path goes here
dir = '/' # directory containing sequence files goes here
csv_file = 'read_counts_derived.csv'

num_lines = sum(1 for line in open(path+dir+csv_file))

num_instances = num_lines//10000

f = open(path+'/bin_starter.sh', 'w')

f.write('#!/bin/bash \n#SBATCH --account=BIOL-CHONG-2019\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=jpjc1@york.ac.uk\n')
f.write('#SBATCH --output=binfeed_%j.log\n')
f.write('#SBATCH --array=1-'+str(num_instances)+'\n#SBATCH --job-name=b_f_%j\n#SBATCH --workdir='+path+dir+'\n\n')
f.write('mkdir bins\n\n')
f.write('dir=$(cd -P -- \"$(dirname -- \"$0\")\" && pwd -P)\n\n')
f.write('STAR=$((SLURM_ARRAY_TASK_ID-1))\n\n')
f.write('START=$((STAR*10000))\n\n')
f.write('END=$((START+9999))\n\n')
f.write('if (($SLURM_ARRAY_TASK_ID == '+str(num_instances)+')); then\n')
f.write('    END=$(wc -l < '+csv_file+')\n')
f.write('fi\n\n')
f.write('module load lang/Python/3.6.6-foss-2018b\n')
f.write('python bin_feeder2.py -i '+csv_file+' -o bins/output.$SLURM_ARRAY_TASK_ID -s $START -e $END')
f.close()
