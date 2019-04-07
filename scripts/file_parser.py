#!/usr/bin/env python3

# file_parser.py
# a modification of the
# new code used to produce output files from clustered contigs
# JC 30/12/15  modified 13/05/16
# 19/2/17 added files from 2k analysis

# code requirements

import json
import re
import csv as csv
import argparse
from Bio import SeqIO
from Bio import SeqUtils as su
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# file names - change these as required
parser = argparse.ArgumentParser(description='')
parser.add_argument('contigs', help='the files', type=str)
parser.add_argument('csv', help='read_counts_derived', type=str)
parser.add_argument('clusters', help='output from step3', type=str)
parser.add_argument('output', help='output directory', type=str)
parser.add_argument('stats', help='stats file', type=str)
parser.add_argument('-l', '--header-list', dest='header', nargs='+', default=[])
args = parser.parse_args()
contig_file = args.contigs
csv_file = args.csv
cluster_file = args.clusters
wd = args.output
stats_file = args.stats
header = args.header

#add context to header columns
header = ['contig'] + header + ['cover', 'length', 'GC']
print(header) #testing...

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
            csv_writer.writerow(header)
            for cluster_name in current_cluster:
                csv_string = cluster_name, ', '.join(map(str, bun_dict[cluster_name])), len(contig_dict[cluster_name]),su.GC(contig_dict[cluster_name].seq)
                csv_writer.writerow(csv_string)
                fasta_entry.append(contig_dict[cluster_name])
            SeqIO.write(fasta_entry, fasta_cluster_filename, 'fasta')
