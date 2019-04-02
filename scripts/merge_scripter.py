# python script to generate parallel_step2 shell script
# directory for sequence files is assumed to be ..path../infiles

import os
import glob

path = os.path.dirname(os.path.realpath(__file__)) # path goes here
dir = '/' # directory containing sequence files goes here
csv_file = 'read_counts_derived.csv'

num_lines = sum(1 for line in open(path+dir+csv_file))

num_instances = str(num_lines//10000)

f = open(path+'/parallel_merge.sh', 'w')

f.write('#!/bin/bash \n#SBATCH --account=BIOL-CHONG-2019\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=jpjc1@york.ac.uk\n')
f.write('#SBATCH --output=merger_%j.log\n')
f.write('#SBATCH --mem=64G\n#SBATCH --job-name=merger\n#SBATCH --workdir='+path+dir+'\n\n')
f.write('module load lang/Python/3.6.6-foss-2018b\n')
f.write('python parallel_merge_step2.py -i bins/parallel_sets -o bins/parallel_merged.out -n '+num_instances)
