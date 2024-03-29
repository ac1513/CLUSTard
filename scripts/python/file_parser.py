#!/usr/bin/env python3

# file_parser.py
# a modification of the
# new code used to produce output files from clustered contigs
# JC 30/12/15  modified 13/05/16
# Generates sequence files for each cluster

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
parser.add_argument('-l', '--header-list', dest='header', nargs='+', default=[])
args = parser.parse_args()
contig_file = args.contigs
csv_file = args.csv
cluster_file = args.clusters
wd = str("output/" + args.output)
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

# open .csv file and store as list

print('Loading abundance data from .csv file')

with open(csv_file, 'r') as abundance:
    bun_entry = csv.reader(abundance)
    bun_list = list(bun_entry)
    for bun_record in bun_list:
        bun_dict[bun_record[0]] = bun_record[1:]

print('Loading cluster details')

with open(cluster_file, 'r') as clusters:
    working_cluster = json.load(clusters)
    for current_cluster in working_cluster:
        fasta_entry = []
        top_len = 0
        for cluster_name in current_cluster:
            fasta_entry.append(contig_dict[cluster_name])
            if len(contig_dict[cluster_name]) > top_len:
                top_len = len(contig_dict[cluster_name])
                top_cluster = cluster_name
        fasta_cluster_filename = (wd+'Cluster_'+str(top_cluster)+'.fasta')
        cluster_filename = (wd+'Cluster_'+str(top_cluster)+'.csv')
        SeqIO.write(fasta_entry, fasta_cluster_filename, 'fasta')
        with open(cluster_filename, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE, escapechar=' ')
            csv_writer.writerow(header)
            for cluster_name in current_cluster:
                csv_string = cluster_name, ', '.join(map(str, bun_dict[cluster_name])), len(contig_dict[cluster_name]),su.GC(contig_dict[cluster_name].seq)
                csv_writer.writerow(csv_string)
