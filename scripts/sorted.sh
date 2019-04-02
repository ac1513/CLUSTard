#!/bin/bash
#SBATCH --account=BIOL-CHONG-2019
#SBATCH --job-name=sorter
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpjc1@york.ac.uk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --output=sorter_%j.log  
 
module load bio/SAMtools/1.9-foss-2018b

cd /mnt/lustre/groups/biol-chong-2019/2019-02-05-CRISPRAD/james_madness

samtools index mapped_NAB1_T00_pe.sorted.bam 
samtools idxstats mapped_NAB1_T00_pe.sorted.bam > counts_NAB1_T00.txt

samtools index mapped_NAB1_T01_pe.sorted.bam 
samtools idxstats mapped_NAB1_T01_pe.sorted.bam > counts_NAB1_T01.txt

samtools index mapped_NAB1_T02_pe.sorted.bam 
samtools idxstats mapped_NAB1_T02_pe.sorted.bam > counts_NAB1_T02.txt

samtools index mapped_NAB1_T03_pe.sorted.bam 
samtools idxstats mapped_NAB1_T03_pe.sorted.bam > counts_NAB1_T03.txt

samtools index mapped_NAB1_T04_pe.sorted.bam 
samtools idxstats mapped_NAB1_T04_pe.sorted.bam > counts_NAB1_T04.txt

samtools index mapped_NAB1_T05_pe.sorted.bam 
samtools idxstats mapped_NAB1_T05_pe.sorted.bam > counts_NAB1_T05.txt

samtools index mapped_NAB1_T06_pe.sorted.bam 
samtools idxstats mapped_NAB1_T06_pe.sorted.bam > counts_NAB1_T06.txt

samtools index mapped_NAB1_T07_pe.sorted.bam 
samtools idxstats mapped_NAB1_T07_pe.sorted.bam > counts_NAB1_T07.txt

samtools index mapped_NAB1_T08_pe.sorted.bam 
samtools idxstats mapped_NAB1_T08_pe.sorted.bam > counts_NAB1_T08.txt

samtools index mapped_NAB1_T09_pe.sorted.bam 
samtools idxstats mapped_NAB1_T09_pe.sorted.bam > counts_NAB1_T09.txt

samtools index mapped_NAB1_T10_pe.sorted.bam 
samtools idxstats mapped_NAB1_T10_pe.sorted.bam > counts_NAB1_T10.txt

samtools index mapped_NAB1_T11_pe.sorted.bam 
samtools idxstats mapped_NAB1_T11_pe.sorted.bam > counts_NAB1_T11.txt

samtools index mapped_NAB1_T12_pe.sorted.bam 
samtools idxstats mapped_NAB1_T12_pe.sorted.bam > counts_NAB1_T12.txt

samtools index mapped_NAB1_T13_pe.sorted.bam 
samtools idxstats mapped_NAB1_T13_pe.sorted.bam > counts_NAB1_T13.txt

samtools index mapped_NAB1_T14_pe.sorted.bam 
samtools idxstats mapped_NAB1_T14_pe.sorted.bam > counts_NAB1_T14.txt

samtools index mapped_NAB1_T15_pe.sorted.bam 
samtools idxstats mapped_NAB1_T15_pe.sorted.bam > counts_NAB1_T15.txt

samtools index mapped_NAB1_T16_pe.sorted.bam 
samtools idxstats mapped_NAB1_T16_pe.sorted.bam > counts_NAB1_T16.txt

samtools index mapped_NAB1_T17_pe.sorted.bam 
samtools idxstats mapped_NAB1_T17_pe.sorted.bam > counts_NAB1_T17.txt

samtools index mapped_NAB2_T00_pe.sorted.bam 
samtools idxstats mapped_NAB2_T00_pe.sorted.bam > counts_NAB2_T00.txt

samtools index mapped_NAB2_T01_pe.sorted.bam 
samtools idxstats mapped_NAB2_T01_pe.sorted.bam > counts_NAB2_T01.txt

samtools index mapped_NAB2_T02_pe.sorted.bam 
samtools idxstats mapped_NAB2_T02_pe.sorted.bam > counts_NAB2_T02.txt

samtools index mapped_NAB2_T03_pe.sorted.bam 
samtools idxstats mapped_NAB2_T03_pe.sorted.bam > counts_NAB2_T03.txt

samtools index mapped_NAB2_T04_pe.sorted.bam 
samtools idxstats mapped_NAB2_T04_pe.sorted.bam > counts_NAB2_T04.txt

samtools index mapped_NAB2_T05_pe.sorted.bam 
samtools idxstats mapped_NAB2_T05_pe.sorted.bam > counts_NAB2_T05.txt

samtools index mapped_NAB2_T06_pe.sorted.bam 
samtools idxstats mapped_NAB2_T06_pe.sorted.bam > counts_NAB2_T06.txt

samtools index mapped_NAB2_T07_pe.sorted.bam 
samtools idxstats mapped_NAB2_T07_pe.sorted.bam > counts_NAB2_T07.txt

samtools index mapped_NAB2_T08_pe.sorted.bam 
samtools idxstats mapped_NAB2_T08_pe.sorted.bam > counts_NAB2_T08.txt

samtools index mapped_NAB2_T09_pe.sorted.bam 
samtools idxstats mapped_NAB2_T09_pe.sorted.bam > counts_NAB2_T09.txt

samtools index mapped_NAB2_T10_pe.sorted.bam 
samtools idxstats mapped_NAB2_T10_pe.sorted.bam > counts_NAB2_T10.txt

samtools index mapped_NAB2_T11_pe.sorted.bam 
samtools idxstats mapped_NAB2_T11_pe.sorted.bam > counts_NAB2_T11.txt

samtools index mapped_NAB2_T12_pe.sorted.bam 
samtools idxstats mapped_NAB2_T12_pe.sorted.bam > counts_NAB2_T12.txt

samtools index mapped_NAB2_T13_pe.sorted.bam 
samtools idxstats mapped_NAB2_T13_pe.sorted.bam > counts_NAB2_T13.txt

samtools index mapped_NAB2_T14_pe.sorted.bam 
samtools idxstats mapped_NAB2_T14_pe.sorted.bam > counts_NAB2_T14.txt

samtools index mapped_NAB2_T15_pe.sorted.bam 
samtools idxstats mapped_NAB2_T15_pe.sorted.bam > counts_NAB2_T15.txt

samtools index mapped_NAB2_T16_pe.sorted.bam 
samtools idxstats mapped_NAB2_T16_pe.sorted.bam > counts_NAB2_T16.txt

samtools index mapped_NAB2_T17_pe.sorted.bam 
samtools idxstats mapped_NAB2_T17_pe.sorted.bam > counts_NAB2_T17.txt

samtools index mapped_NAB3_T00_pe_sorted.bam 
samtools idxstats mapped_NAB3_T00_pe_sorted.bam > counts_NAB3_T00.txt

samtools index mapped_NAB3_T01_pe_sorted.bam 
samtools idxstats mapped_NAB3_T01_pe_sorted.bam > counts_NAB3_T01.txt

samtools index mapped_NAB3_T02_pe_sorted.bam 
samtools idxstats mapped_NAB3_T02_pe_sorted.bam > counts_NAB3_T02.txt

samtools index mapped_NAB3_T03_pe_sorted.bam 
samtools idxstats mapped_NAB3_T03_pe_sorted.bam > counts_NAB3_T03.txt

samtools index mapped_NAB3_T04_pe_sorted.bam 
samtools idxstats mapped_NAB3_T04_pe_sorted.bam > counts_NAB3_T04.txt

samtools index mapped_NAB3_T05_pe_sorted.bam 
samtools idxstats mapped_NAB3_T05_pe_sorted.bam > counts_NAB3_T05.txt

samtools index mapped_NAB3_T06_pe_sorted.bam 
samtools idxstats mapped_NAB3_T06_pe_sorted.bam > counts_NAB3_T06.txt

samtools index mapped_NAB3_T07_pe_sorted.bam 
samtools idxstats mapped_NAB3_T07_pe_sorted.bam > counts_NAB3_T07.txt

samtools index mapped_NAB3_T08_pe_sorted.bam 
samtools idxstats mapped_NAB3_T08_pe_sorted.bam > counts_NAB3_T08.txt

samtools index mapped_NAB3_T09_pe_sorted.bam 
samtools idxstats mapped_NAB3_T09_pe_sorted.bam > counts_NAB3_T09.txt

samtools index mapped_NAB3_T10_pe_sorted.bam 
samtools idxstats mapped_NAB3_T10_pe_sorted.bam > counts_NAB3_T10.txt

samtools index mapped_NAB3_T11_pe_sorted.bam 
samtools idxstats mapped_NAB3_T11_pe_sorted.bam > counts_NAB3_T11.txt

samtools index mapped_NAB3_T12_pe_sorted.bam 
samtools idxstats mapped_NAB3_T12_pe_sorted.bam > counts_NAB3_T12.txt

samtools index mapped_NAB3_T13_pe_sorted.bam 
samtools idxstats mapped_NAB3_T13_pe_sorted.bam > counts_NAB3_T13.txt

samtools index mapped_NAB3_T14_pe_sorted.bam 
samtools idxstats mapped_NAB3_T14_pe_sorted.bam > counts_NAB3_T14.txt

samtools index mapped_NAB3_T15_pe_sorted.bam 
samtools idxstats mapped_NAB3_T15_pe_sorted.bam > counts_NAB3_T15.txt

samtools index mapped_NAB3_T16_pe_sorted.bam 
samtools idxstats mapped_NAB3_T16_pe_sorted.bam > counts_NAB3_T16.txt

samtools index mapped_NAB3_T17_pe_sorted.bam 
samtools idxstats mapped_NAB3_T17_pe_sorted.bam > counts_NAB3_T17.txt

samtools index mapped_NAB4_T00_pe_sorted.bam 
samtools idxstats mapped_NAB4_T00_pe_sorted.bam > counts_NAB4_T00.txt

samtools index mapped_NAB4_T01_pe_sorted.bam 
samtools idxstats mapped_NAB4_T01_pe_sorted.bam > counts_NAB4_T01.txt

samtools index mapped_NAB4_T02_pe_sorted.bam 
samtools idxstats mapped_NAB4_T02_pe_sorted.bam > counts_NAB4_T02.txt

samtools index mapped_NAB4_T03_pe_sorted.bam 
samtools idxstats mapped_NAB4_T03_pe_sorted.bam > counts_NAB4_T03.txt

samtools index mapped_NAB4_T04_pe_sorted.bam 
samtools idxstats mapped_NAB4_T04_pe_sorted.bam > counts_NAB4_T04.txt

samtools index mapped_NAB4_T05_pe_sorted.bam 
samtools idxstats mapped_NAB4_T05_pe_sorted.bam > counts_NAB4_T05.txt

samtools index mapped_NAB4_T06_pe_sorted.bam 
samtools idxstats mapped_NAB4_T06_pe_sorted.bam > counts_NAB4_T06.txt

samtools index mapped_NAB4_T07_pe_sorted.bam 
samtools idxstats mapped_NAB4_T07_pe_sorted.bam > counts_NAB4_T07.txt

samtools index mapped_NAB4_T08_pe_sorted.bam 
samtools idxstats mapped_NAB4_T08_pe_sorted.bam > counts_NAB4_T08.txt

samtools index mapped_NAB4_T09_pe_sorted.bam 
samtools idxstats mapped_NAB4_T09_pe_sorted.bam > counts_NAB4_T09.txt

samtools index mapped_NAB4_T10_pe_sorted.bam 
samtools idxstats mapped_NAB4_T10_pe_sorted.bam > counts_NAB4_T10.txt

samtools index mapped_NAB4_T11_pe_sorted.bam 
samtools idxstats mapped_NAB4_T11_pe_sorted.bam > counts_NAB4_T11.txt

samtools index mapped_NAB4_T12_pe_sorted.bam 
samtools idxstats mapped_NAB4_T12_pe_sorted.bam > counts_NAB4_T12.txt

samtools index mapped_NAB4_T13_pe_sorted.bam 
samtools idxstats mapped_NAB4_T13_pe_sorted.bam > counts_NAB4_T13.txt

samtools index mapped_NAB4_T14_pe_sorted.bam 
samtools idxstats mapped_NAB4_T14_pe_sorted.bam > counts_NAB4_T14.txt

samtools index mapped_NAB4_T15_pe_sorted.bam 
samtools idxstats mapped_NAB4_T15_pe_sorted.bam > counts_NAB4_T15.txt

samtools index mapped_NAB4_T16_pe_sorted.bam 
samtools idxstats mapped_NAB4_T16_pe_sorted.bam > counts_NAB4_T16.txt

samtools index mapped_NAB4_T17_pe_sorted.bam 
samtools idxstats mapped_NAB4_T17_pe_sorted.bam > counts_NAB4_T17.txt

samtools index mapped_NAB_Feed_T09_pe.sorted.bam 
samtools idxstats mapped_NAB_Feed_T09_pe.sorted.bam > counts_NAB_Feed_T09.txt

samtools index mapped_NAB_Feed_T10_pe.sorted.bam 
samtools idxstats mapped_NAB_Feed_T10_pe.sorted.bam > counts_NAB_Feed_T10.txt

samtools index mapped_NAB_Feed_T11_pe.sorted.bam 
samtools idxstats mapped_NAB_Feed_T11_pe.sorted.bam > counts_NAB_Feed_T11.txt

samtools index mapped_NAB_Feed_T12_pe.sorted.bam 
samtools idxstats mapped_NAB_Feed_T12_pe.sorted.bam > counts_NAB_Feed_T12.txt

samtools index mapped_NAB_Feed_T13_pe.sorted.bam 
samtools idxstats mapped_NAB_Feed_T13_pe.sorted.bam > counts_NAB_Feed_T13.txt

samtools index mapped_NAB_Feed_T14_pe.sorted.bam 
samtools idxstats mapped_NAB_Feed_T14_pe.sorted.bam > counts_NAB_Feed_T14.txt

samtools index mapped_NAB_Feed_T15_pe.sorted.bam 
samtools idxstats mapped_NAB_Feed_T15_pe.sorted.bam > counts_NAB_Feed_T15.txt

samtools index mapped_NAB_Feed_T17_pe.sorted.bam 
samtools idxstats mapped_NAB_Feed_T17_pe.sorted.bam > counts_NAB_Feed_T17.txt

