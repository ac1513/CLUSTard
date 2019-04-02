# python script to generate the shell script to map original reads to indexed contigs
# this version no longer requires manual entry of path information
# it assumes sequence files are in ..path../infiles directory

import os
import glob

path = '/users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/data' # path goes here
dir = '/illumina' # directory containing sequence files goes here

os.chdir(path+dir) # cunningly changes working directory

a = sorted(glob.glob('*_R1_001.fastq.gz'))
b = sorted(glob.glob('*_R2_001.fastq.gz'))
c = str(len(a))

w = sorted(glob.glob('*_1.fastq.gz'))
x = sorted(glob.glob('*_2.fastq.gz'))
y = str(len(w))

f = open('/users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapping.sh', 'w')

f.write('#!/bin/sh \n#$ -N st_mapper \n#$ -wd '+path+dir+'\n#$ -pe smp 16\n#$ -m beas\n\n')

for i in range(len(a)):
   f.write('bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta '+a[i]+' '+b[i]+' \\\n')
   f.write('| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_'+a[i][:10]+'_pe.sorted.bam -@ 16 \n')

for j in range(len(w)):
   f.write('bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta '+w[j]+' '+x[j]+' \\\n')
   f.write('| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_'+w[j][:10]+'_pe.sorted.bam -@ 16 \n')

f.close()

