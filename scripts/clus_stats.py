#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 19:48:34 2019

@author: ac1513
"""

import statistics
import pandas as pd
import argparse

# =============================================================================
# Command line parsing
# =============================================================================
parser = argparse.ArgumentParser(description='usage = python in_file prefix -k kraken_file -k_l kraken_level -cm checkm -sk seqkit_file samples_file date(y/n)')
parser.add_argument('in_file', help='the name of the file containing a list of csv files', type=str)
parser.add_argument('prefix', help='prefix of the jobs', type=str)
parser.add_argument('-cm', '--checkm_file', help = 'checkm output file - in tab format', type = str)
parser.add_argument('-sk', '--seqkit', help = 'seqkit output file - in tab format', type = str)

args = parser.parse_args()
in_file = args.in_file
prefix = args.prefix

# =============================================================================
# Read in input file(s)
# =============================================================================

set_groups = set()

with open(in_file, 'r') as text_file:
    files = text_file.read().strip().split()

# =============================================================================
# Stats
# =============================================================================
stats_df = pd.DataFrame(columns=['no_seq','tot_len','av_cov','sd_cov', 'av_gc','sd_gc','n_50','comp', 'contam'])

for file in files:
    with open(file, 'r') as f:
        file = str(file)
        df = pd.read_csv(f, index_col='contig')

        av_cov = str('{0:.1f}'.format(statistics.mean(df['cover'].tolist())))
        sd_cov = str('{0:.1f}'.format(statistics.stdev(df['cover'].tolist())))
        tot_len = str(sum(df['length'].tolist()))
        av_gc =  str('{0:.1f}'.format(statistics.mean(df['GC'].tolist())))
        sd_gc = str('{0:.1f}'.format(statistics.stdev(df['GC'].tolist())))
        na = str(file.split('/')[-1:][0].split('.')[0])
        nu = str(len(df))

        if args.seqkit:
            seqkit_df = pd.read_csv(args.seqkit, sep = '\t', index_col =0)
            file_fa = file.replace(".csv",".fasta")
            n_50 = seqkit_df["N50"][file_fa]

        if args.checkm_file:
            checkm_df = pd.read_csv(args.checkm_file, sep = '\t', index_col = 0)
            clus = file.split('/')[-1:][0][:-4]
            comp = checkm_df["Completeness"][clus]
            conta = checkm_df["Contamination"][clus]

        stats_df.loc[na] = [nu,tot_len,av_cov,sd_cov,av_gc,sd_gc,n_50,comp,conta]

stats_df.to_csv("output/" + prefix + "_cluster_summary_stats.tsv", sep='\t')
