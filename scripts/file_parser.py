# file_parser.py
# a modification of the
# new code used to produce output files from clustered contigs
# JC 30/12/15  modified 13/05/16
# 19/2/17 added files from 2k analysis

# code requirements

import json
import re
import csv as csv
from Bio import SeqIO
from Bio import SeqUtils as su
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# file names - change these as required

contig_file = '/users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/anvio/yw_polished_anvio.fasta'
csv_file = 'read_counts_derived.csv'
cluster_file = 'bins/total_step3_list.out' # efforts from clustering
stats_file = 'summary_stats.csv'	# this is where some results are written
wd = '/users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/results/'	# directory to save results


# dictionaries and lists

cluster_stats = []     # list of stats on cluster data for export to .csv
bun_dict = {}


# open .fasta file containing contigs and store as dict

print('Opening contig sequence file')
contig_dict = SeqIO.to_dict(SeqIO.parse(contig_file, "fasta"))
print(str(len(contig_dict))+' sequences loaded')


# open .csv file and store as list(?)

print('Loading abundance data from .csv file')

with open(csv_file, 'r') as abundance:
    bun_entry = csv.reader(abundance)
    bun_list = list(bun_entry)
    for bun_record in bun_list:
        bun_dict[bun_record[0]] = bun_record[1:]


# open .list file of contig i.d.s, loop over this list then
# generate separate .csv and .fasta files for each cluster
# .csv file should contain contig names and lengths as well as GC content
#	NOTE: the structure of these files has changed
#	      it's better, but different!

print('Loading cluster details')

with open(cluster_file, 'r') as clusters:
    working_cluster = json.load(clusters)
    for current_cluster in working_cluster:
        print('new cluster')
        cluster_filename = (wd+'Cluster_'+str(current_cluster[0])+'.csv')
        fasta_cluster_filename = (wd+'Cluster_'+str(current_cluster[0])+'.fasta')
        fasta_entry = []
        with open(cluster_filename, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE, escapechar=' ')
            csv_writer.writerow(['contig','1_00','1_01','1_02','1_03','1_04','1_05','1_06','1_07','1_08','1_09','1_10','1_11','1_12','1_13','1_14','1_15','1_16','1_17','2_00','2_01','2_02','2_03','2_04','2_05','2_06','2_07','2_08','2_09','2_10','2_11','2_12','2_13','2_14','2_15','2_16','2_17','3_00','3_01,','3_02','3_03','3_04','3_05','3_06','3_07','3_08','3_09','3_10','3_11','3_12','3_13','3_14','3_15','3_16','3_17','4_00','4_01','4_02','4_03','4_04','4_05','4_06','4_07','4_08','4_09','4_10','4_11','4_12','4_13','4_14','4_15','4_16','4_17','F_09','F_10','F_11','F_12','F_13','F_14','F_15','F_17','cover','length','GC'])
            for cluster_name in current_cluster:
                csv_string = cluster_name, ', '.join(map(str, bun_dict[cluster_name])), len(contig_dict[cluster_name]),su.GC(contig_dict[cluster_name].seq)
                csv_writer.writerow(csv_string)
                fasta_entry.append(contig_dict[cluster_name])
            SeqIO.write(fasta_entry, fasta_cluster_filename, 'fasta')
