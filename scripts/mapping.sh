#!/bin/bash
#SBATCH --account=BIOL-CHONG-2019
#SBATCH --job-name=mapper
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpjc1@york.ac.uk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --output=mapper_%j.log  
 
module load bio/BWA/0.7.17-foss-2018b
module load bio/SAMtools/1.9-foss-2018b

cd /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/data/illumina

bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_Foam_S64_L003_R1_001.fastq.gz NAB1_Foam_S64_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_Foam__pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T0_S44_L003_R1_001.fastq.gz NAB1_T0_S44_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T0_S4_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T1_S48_L003_R1_001.fastq.gz NAB1_T1_S48_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T1_S4_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T2_S52_L003_R1_001.fastq.gz NAB1_T2_S52_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T2_S5_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T3_S56_L003_R1_001.fastq.gz NAB1_T3_S56_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T3_S5_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T4_S60_L003_R1_001.fastq.gz NAB1_T4_S60_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T4_S6_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_Foam_S65_L003_R1_001.fastq.gz NAB2_Foam_S65_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_Foam__pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T0_S45_L003_R1_001.fastq.gz NAB2_T0_S45_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T0_S4_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T1_S49_L003_R1_001.fastq.gz NAB2_T1_S49_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T1_S4_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T2_S53_L003_R1_001.fastq.gz NAB2_T2_S53_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T2_S5_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T3_S57_L003_R1_001.fastq.gz NAB2_T3_S57_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T3_S5_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T4_S61_L003_R1_001.fastq.gz NAB2_T4_S61_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T4_S6_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_Foam_S66_L003_R1_001.fastq.gz NAB3_Foam_S66_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_Foam__pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T0_S46_L003_R1_001.fastq.gz NAB3_T0_S46_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T0_S4_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T1_S50_L003_R1_001.fastq.gz NAB3_T1_S50_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T1_S5_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T2_S54_L003_R1_001.fastq.gz NAB3_T2_S54_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T2_S5_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T3_S58_L003_R1_001.fastq.gz NAB3_T3_S58_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T3_S5_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T4_S62_L003_R1_001.fastq.gz NAB3_T4_S62_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T4_S6_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_Foam_S67_L003_R1_001.fastq.gz NAB4_Foam_S67_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_Foam__pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T0_S47_L003_R1_001.fastq.gz NAB4_T0_S47_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T0_S4_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T1_S51_L003_R1_001.fastq.gz NAB4_T1_S51_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T1_S5_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T2_S55_L003_R1_001.fastq.gz NAB4_T2_S55_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T2_S5_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T3_S59_L003_R1_001.fastq.gz NAB4_T3_S59_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T3_S5_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T4_S63_L003_R1_001.fastq.gz NAB4_T4_S63_L003_R2_001.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T4_S6_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T10_1.fastq.gz NAB1_T10_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T10_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T11_1.fastq.gz NAB1_T11_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T11_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T12_1.fastq.gz NAB1_T12_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T12_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T13_1.fastq.gz NAB1_T13_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T13_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T14_1.fastq.gz NAB1_T14_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T14_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T15_1.fastq.gz NAB1_T15_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T15_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T16_1.fastq.gz NAB1_T16_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T16_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T17_1.fastq.gz NAB1_T17_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T17_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T5_1.fastq.gz NAB1_T5_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T5_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T6_1.fastq.gz NAB1_T6_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T6_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T7_1.fastq.gz NAB1_T7_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T7_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T8_1.fastq.gz NAB1_T8_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T8_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB1_T9_1.fastq.gz NAB1_T9_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB1_T9_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T10_1.fastq.gz NAB2_T10_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T10_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T11_1.fastq.gz NAB2_T11_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T11_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T12_1.fastq.gz NAB2_T12_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T12_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T13_1.fastq.gz NAB2_T13_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T13_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T14_1.fastq.gz NAB2_T14_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T14_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T15_1.fastq.gz NAB2_T15_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T15_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T16_1.fastq.gz NAB2_T16_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T16_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T17_1.fastq.gz NAB2_T17_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T17_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T5_1.fastq.gz NAB2_T5_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T5_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T6_1.fastq.gz NAB2_T6_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T6_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T7_1.fastq.gz NAB2_T7_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T7_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T8_1.fastq.gz NAB2_T8_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T8_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB2_T9_1.fastq.gz NAB2_T9_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB2_T9_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T10_1.fastq.gz NAB3_T10_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T10_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T11_1.fastq.gz NAB3_T11_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T11_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T12_1.fastq.gz NAB3_T12_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T12_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T13_1.fastq.gz NAB3_T13_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T13_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T14_1.fastq.gz NAB3_T14_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T14_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T15_1.fastq.gz NAB3_T15_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T15_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T16_1.fastq.gz NAB3_T16_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T16_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T17_1.fastq.gz NAB3_T17_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T17_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T5_1.fastq.gz NAB3_T5_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T5_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T6_1.fastq.gz NAB3_T6_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T6_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T7_1.fastq.gz NAB3_T7_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T7_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T8_1.fastq.gz NAB3_T8_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T8_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB3_T9_1.fastq.gz NAB3_T9_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB3_T9_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T10_1.fastq.gz NAB4_T10_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T10_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T11_1.fastq.gz NAB4_T11_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T11_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T12_1.fastq.gz NAB4_T12_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T12_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T13_1.fastq.gz NAB4_T13_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T13_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T14_1.fastq.gz NAB4_T14_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T14_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T15_1.fastq.gz NAB4_T15_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T15_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T16_1.fastq.gz NAB4_T16_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T16_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T17_1.fastq.gz NAB4_T17_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T17_1_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T5_1.fastq.gz NAB4_T5_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T5_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T6_1.fastq.gz NAB4_T6_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T6_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T7_1.fastq.gz NAB4_T7_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T7_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T8_1.fastq.gz NAB4_T8_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T8_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB4_T9_1.fastq.gz NAB4_T9_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB4_T9_1._pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T10_1.fastq.gz NAB_Feed_T10_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T11_1.fastq.gz NAB_Feed_T11_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T12_1.fastq.gz NAB_Feed_T12_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T13_1.fastq.gz NAB_Feed_T13_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T14_1.fastq.gz NAB_Feed_T14_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T15_1.fastq.gz NAB_Feed_T15_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T17_1.fastq.gz NAB_Feed_T17_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T_pe.sorted.bam -@ 16 
bwa mem -M -t 16 /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/asm/nano/yw_polished_anvio.fasta NAB_Feed_T9_1.fastq.gz NAB_Feed_T9_2.fastq.gz \
| samtools view -buS - | samtools sort -o /users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/mapped_NAB_Feed_T_pe.sorted.bam -@ 16 
