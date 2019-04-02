# python script to generate the shell file to index reads to produce sample counts
# this version no longer requires manual entry of path information
# it assumes sequence files are in ..path../in_files directory

import os
import glob

path = os.path.dirname(os.path.realpath(__file__)) # path goes here
#dir = '/short_reads' # directory containing sequence files goes here

os.chdir(path) # cunningly changes working directory

a = sorted(glob.glob('*.bam'))

f = open(path+'/sorted.sh', 'w')

f.write('#!/bin/sh \n#$ -N sorter \n#$ -wd '+path+'\n\n')

for i in a:
   f.write('samtools index '+i+' \nsamtools idxstats '+i+' > counts_'+i[7:15]+'.txt\n\n')
