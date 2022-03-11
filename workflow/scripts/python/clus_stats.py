#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 19:48:34 2019

@author: ac1513
"""

import statistics
import pandas as pd
import argparse
import re
import json

# =============================================================================
# Command line parsing
# =============================================================================
parser = argparse.ArgumentParser(description='usage = python in_file prefix -k kraken_file -k_l kraken_level -cm checkm -sk seqkit_file samples_file date(y/n)')
parser.add_argument('in_file', help='the name of the file containing a list of csv files', type=str)
parser.add_argument('prefix', help='prefix of the jobs', type=str)
parser.add_argument('-cm', '--checkm_file', help = 'checkm output file - in tab format', type = str)
parser.add_argument('-sk', '--seqkit', help = 'seqkit output file - in tab format', type = str)
parser.add_argument('-k', '--kraken', help = 'kraken output file at selected genus level', type = str)
parser.add_argument('-g', '--gtdb', help = 'json file of GTDB lookup', type = str)
parser.add_argument('-q', '--qual', help = 'clustard output containing quality of mags', type = str)
parser.add_argument('-out', '--output', help = 'where to save output', type = str)

args = parser.parse_args()
in_file = args.in_file
prefix = args.prefix
output = args.output

# =============================================================================
# Read in input file(s)
# =============================================================================

set_groups = set()

with open(in_file, 'r') as text_file:
    files = text_file.read().strip().split()

gtdb_lookup = ""
if args.gtdb:
    with open(args.gtdb, 'r') as read_file:
        empty_test = read_file.read(1)
        if empty_test:
            read_file.seek(0)
            gtdb_lookup = json.load(read_file)
    if len(gtdb_lookup) == 0:
        gtdb_ncbi = " "

# =============================================================================
# Stats
# =============================================================================
stats_df = pd.DataFrame(columns=['no_seq','tot_len','av_cov','sd_cov', 'av_gc','sd_gc','n_50','comp', 'contam', 'kraken_id', 'kraken_per', 'gtdb_ncbi','qual'])


for file in files:
    with open(file, 'r') as f:
        file = str(file)
        df = pd.read_csv(f, index_col='contig')

        av_cov = str('{0:.1f}'.format(statistics.mean(df['cover'].tolist())))
        if len(df['cover']) > 1:
            sd_cov = str('{0:.1f}'.format(statistics.stdev(df['cover'].tolist())))
        else:
            sd_cov = 0
        tot_len = str(sum(df['length'].tolist()))
        av_gc =  str('{0:.1f}'.format(statistics.mean(df['GC'].tolist())))
        if len(df['GC']) > 1:
            sd_gc = str('{0:.1f}'.format(statistics.stdev(df['GC'].tolist())))
        else:
            sd_gc = 0
        na = str(file.split('/')[-1:][0][:-4])
        nu = str(len(df))

        if args.seqkit:
            seqkit_df = pd.read_csv(args.seqkit, sep = '\t', index_col =0)
            file_fa = file.replace(".csv",".fasta")
            n_50 = seqkit_df["N50"][file_fa]
            tot_len = str(seqkit_df["sum_len"][file_fa])
            nu = str(seqkit_df["num_seqs"][file_fa])

        if args.checkm_file:
            checkm_df = pd.read_csv(args.checkm_file, sep = '\t', index_col = 0)
            clus = file.split('/')[-1:][0][:-4]
            comp = checkm_df["Completeness"][clus]
            conta = checkm_df["Contamination"][clus]

        if args.kraken:
            for line in open(args.kraken, 'r'):
                if re.search(na+'_', line):
                    krak_id = line.split('\t')[-1].strip()
                    if './' in krak_id:
                        cont = ' '
                    else:
                        krak_per = line.split('\t')[1]
                        krak_per = krak_per.strip()
                        if len(gtdb_lookup) > 0:
                            if krak_id in gtdb_lookup:
                                gtdb_ncbi = gtdb_lookup[krak_id]
                            else:
                                gtdb_ncbi = "Missing"

        if args.qual:
            for line in open(args.qual, 'r'):
                if re.search(na, line):
                    qual = line.split(":")[0]
        else:
            qual = "N/A"

        stats_df.loc[na] = [nu,tot_len,av_cov,sd_cov,av_gc,sd_gc,n_50,comp,conta,krak_id,krak_per,gtdb_ncbi,qual]

stats_df.to_csv(output, sep='\t')
