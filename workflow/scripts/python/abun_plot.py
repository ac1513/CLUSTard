#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 19:38:16 2019

@author: ac1513
"""

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10})
import argparse
import pandas as pd
import re

# =============================================================================
# PLOTTING FUNCTION
# =============================================================================
def plotting(abun, top20, plot):
    
    if 'y' not in top20:
        col_zero = abun.columns[0]
        abun.sort_values(col_zero, inplace=True, ascending=False)
        
    abun_flip = abun.transpose()
        
    if 'y' in top20:
        plot = plot + '_top20'
        colours = ['#502db3', '#008080', '#c200f2', '#f2c200','#36a3d9', '#e6beff','#8c0025', '#f58231', '#bf0080', '#cad900','#911eb4','#e5001f','#0066bf', '#000075','#338000', '#f032e6','#1bca00','#1d4010','#9a6324','#a9a9a9']
        abun_flip.plot.bar(stacked=True, legend = None, figsize=(15,10), color=colours, width=0.9)
    else:
        abun_flip.plot.bar(stacked=True, legend = None, figsize=(15,10), width=0.9)
    
    if 'a' in plot:
        plt.margins(x = 0.01, y=0.05)
        plt.ylabel("Absolute abundance", fontsize = 18)
    if 'r' in plot:
        plt.margins(0)
        plt.ylim(0,100)
        plt.ylabel("Relative abundance (%)", fontsize = 18)

    plt.xlabel("Samples", fontsize = 18)
    
    plt.xticks(rotation='vertical', fontsize = 10, verticalalignment='center_baseline')
    plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1, frameon=False, fontsize = 10)

    plt.tight_layout()
    plt.savefig(out_loc + 'plots/' + prefix +'_' + plot +'_'+ 'abun_plot.png', bbox_inches='tight', dpi = 300)
        
# =============================================================================
# Command line parsing
# =============================================================================

parser = argparse.ArgumentParser(description='usage = python prefix csv_in binned_in plot(r/a) top20(y/n) -s sample_list -k kraken_file ')
parser.add_argument('prefix', help='prefix of the jobs', type=str)
parser.add_argument('csv_in', help='file containing absolute abundance counts for every contig', type=str)
parser.add_argument('binned_in', help='file listing contigs in each cluster', type=str)
parser.add_argument('out_loc', help = 'location of output dir', type = str)
parser.add_argument('-s', '--sample', dest='sample', nargs='+', default=[])
parser.add_argument('-k', '--kraken', dest='kraken', help = 'kraken top output file', type = str)

args = parser.parse_args()

prefix = args.prefix
csv_in = args.csv_in
binned_in = args.binned_in
sample = args.sample
out_loc = args.out_loc

if args.kraken:
    kraken = args.kraken
    taxo = "y"
else:
    taxo = "n"

# =============================================================================
# Read in files
# =============================================================================

df_abun = pd.read_csv(csv_in, index_col = 0, names = sample)
df_abun = df_abun.drop(columns='Coverage')

abun = pd.DataFrame(columns = sample).drop(columns='Coverage')

# =============================================================================
# Parse abundance and taxonomy files
# =============================================================================

with open(binned_in, 'r') as binned_list:
    for line in binned_list:
        if taxo == "y": #parse the kraken info if it's supplied
            if line.startswith(">"): #get cluster info
                cluster = line.strip('>').strip()[:-6] #strip things
                with open(kraken, 'r') as kraken_f:
                    for line2 in kraken_f:
                        if re.search(cluster+'_', line2):
                            name = line2.split('\t')[-1].strip()
                            if name =="":
                                name = 'Unclassified'
            else:
                line = line.strip().split(' ')[0] #split may not work with NAB_997 - check
                if line in df_abun.index:
                    if name in abun.index:
                        abun.loc[name] = abun.loc[name].add(df_abun.loc[line])
                    else:
                        abun.loc[name] = df_abun.loc[line]

        else:
            if line.startswith(">"): #get cluster info
                cluster = line.strip('>').strip()[:-6] #strip things
                name = cluster
            else:
                line = line.strip().split(' ')[0] #split may not work with NAB_997 - check
                if line in df_abun.index:
                    if name in abun.index:
                        abun.loc[name] = abun.loc[name].add(df_abun.loc[line])
                    else:
                        abun.loc[name] = df_abun.loc[line]

# =============================================================================
# Save absolute counts to file and plot
# =============================================================================

absolute_abun = abun.copy()
absolute_abun.to_csv(out_loc + "plots/" + prefix + "_absolute_counts.csv")
plotting(absolute_abun, "n", "a")

# =============================================================================
# Save relative counts to file
# =============================================================================

relative_abun = abun.copy()

for column in abun: #iterate over columns
    per = []
    for val in abun.loc[:, column]:
        per.append((val/abun[column].sum())*100)
    relative_abun[column] = per

relative_abun.index = abun.index.values.tolist()

relative_abun.to_csv(out_loc + "plots/" + prefix + "_relative_counts.csv")
plotting(relative_abun, "n", "r" )

# =============================================================================
# Calculate top 20 species
# =============================================================================
# =============================================================================
# Absolute
# =============================================================================

absolute_abun["sum"] = absolute_abun.sum(axis=1)
top_abs_abun = absolute_abun.sort_values('sum', axis=0, ascending=False).head(20).drop(columns = "sum")
other_abs_abun = absolute_abun.sort_values('sum', axis=0, ascending=False).iloc[20:,]
top_abs_abun.loc["Other"] = other_abs_abun.sum(axis=0).drop(columns="sum")

top_abs_abun.to_csv(out_loc + "plots/" + prefix + "_top20_absolute_counts.csv")
plotting(top_abs_abun, "y", "a")

# =============================================================================
# Relative
# =============================================================================

relative_abun["sum"] = relative_abun.sum(axis=1)
top_rel_abun = relative_abun.sort_values('sum', axis=0, ascending=False).head(20).drop(columns = "sum")
other_rel_abun = relative_abun.sort_values('sum', axis=0, ascending=False).iloc[20:,]
top_rel_abun.loc["Other"] = other_rel_abun.sum(axis=0).drop(columns="sum")

top_rel_abun.to_csv(out_loc +"plots/" + prefix + "_top20_relative_counts.csv")
plotting(top_rel_abun, "y", "r")