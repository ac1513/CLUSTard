#!/bin/bash
#SBATCH --account=BIOL-CHONG-2019
#SBATCH --job-name=mapper
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpjc1@york.ac.uk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --output=mapper_%j.log  
 
module load bio/BWA/0.7.17-foss-2018b
module load bio/SAMtools/1.9-foss-2018b

cd /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/data/illumina

bwa mem -M -t 40 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T10_1.fastq.gz NAB_Feed_T10_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T10_pe.sorted.bam -@ 40 
bwa mem -M -t 40 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T11_1.fastq.gz NAB_Feed_T11_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T11_pe.sorted.bam -@ 40 
bwa mem -M -t 40 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T12_1.fastq.gz NAB_Feed_T12_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T12_pe.sorted.bam -@ 40 
bwa mem -M -t 40 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T13_1.fastq.gz NAB_Feed_T13_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T13_pe.sorted.bam -@ 40 
bwa mem -M -t 40 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T14_1.fastq.gz NAB_Feed_T14_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T14_pe.sorted.bam -@ 40 
bwa mem -M -t 40 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T15_1.fastq.gz NAB_Feed_T15_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T15_pe.sorted.bam -@ 40 
bwa mem -M -t 40 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T17_1.fastq.gz NAB_Feed_T17_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T17_pe.sorted.bam -@ 40 
bwa mem -M -t 40 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T9_1.fastq.gz NAB_Feed_T9_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T09_pe.sorted.bam -@ 40 
