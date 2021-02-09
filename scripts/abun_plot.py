#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 19:38:16 2019

@author: ac1513
"""

import matplotlib.pyplot as plt
import argparse
import pandas as pd
import re

# =============================================================================
# Command line parsing
# =============================================================================

parser = argparse.ArgumentParser(description='usage = python prefix csv_in binned_in plot(r/a) top20(y/n) -s sample_list -k kraken_file ')
parser.add_argument('prefix', help='prefix of the jobs', type=str)
parser.add_argument('csv_in', help='file containing absolute abundance counts for every contig', type=str)
parser.add_argument('binned_in', help='file listing contigs in each cluster', type=str)
parser.add_argument('plot', help='relative or absolute abundance output (r/a)', type=str)
parser.add_argument('top20', help = 'output top19 and other or output all y/n', type = str)
parser.add_argument('-s', '--sample', dest='sample', nargs='+', default=[])
parser.add_argument('-k', '--kraken', dest='kraken', help = 'kraken top output file', type = str)

args = parser.parse_args()

prefix = args.prefix
csv_in = args.csv_in
plot = args.plot
binned_in = args.binned_in
sample = args.sample
top20 = args.top20

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
tot = df_abun.sum(axis = 0)
cluster_abun = pd.DataFrame(columns=df_abun.columns)
new_df_abun = pd.DataFrame(columns = sample).drop(columns='Coverage')
abun = pd.DataFrame(columns = sample).drop(columns='Coverage')

# =============================================================================
# Parse count and taxonomy files
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
                    if name in new_df_abun.index:
                        new_df_abun.loc[name] = new_df_abun.loc[name].add(df_abun.loc[line])
                    else:
                        new_df_abun.loc[name] = df_abun.loc[line]

        else:
            if line.startswith(">"): #get cluster info
                cluster = line.strip('>').strip()[:-6] #strip things
                name = cluster
            else:
                line = line.strip().split(' ')[0] #split may not work with NAB_997 - check
                if line in df_abun.index:
                    if name in new_df_abun.index:
                        new_df_abun.loc[name] = new_df_abun.loc[name].add(df_abun.loc[line])
                    else:
                        new_df_abun.loc[name] = df_abun.loc[line]

# =============================================================================
# Save absolute counts to file
# =============================================================================

absolute_abun = new_df_abun.copy()
absolute_abun.to_csv("output/plots/" + prefix + "_absolute_counts.csv")

# =============================================================================
# Save relative counts to file
# =============================================================================

relative_abun = new_df_abun.copy()

for column in new_df_abun: #iterate over columns
    per = []
    for val in new_df_abun.loc[:, column]:
        per.append((val/new_df_abun[column].sum())*100)
    relative_abun[column] = per
abun_sum = relative_abun.cumsum()
relative_abun.index = new_df_abun.index.values.tolist()
relative_abun.to_csv("output/plots/" + prefix + "_relative_counts.csv")

# =============================================================================
# Calculate top 20 species
# =============================================================================

if top20 == 'y':

    if 'a' in plot:

        absolute_abun["sum"] = absolute_abun.sum(axis=1)
        top_df_abun = absolute_abun.sort_values('sum', axis=0, ascending=False).head(19).drop(columns = "sum")
        other_df_abun = absolute_abun.sort_values('sum', axis=0, ascending=False).iloc[19:,]
        top_df_abun.loc["Other"] = other_df_abun.sum(axis=0).drop(columns="sum")
        abun = top_df_abun.copy()
        abun_sum = abun.cumsum()
        abun.to_csv("output/plots/" + prefix + "_top20_absolute_counts.csv")


    elif 'r' in plot:

        new_df_abun["sum"] = new_df_abun.sum(axis=1)
        top_df_abun = new_df_abun.sort_values('sum', axis=0, ascending=False).head(19).drop(columns = "sum")
        other_df_abun = new_df_abun.sort_values('sum', axis=0, ascending=False).iloc[19:,]
        top_df_abun.loc["Other"] = other_df_abun.sum(axis=0).drop(columns="sum")
        top_relative_abun = top_df_abun.copy()
        for column in top_df_abun: #iterate over columns
            per = []
            for val in top_df_abun.loc[:, column]:
                per.append(((val/top_df_abun[column].sum())*100))
            top_relative_abun[column] = per
        abun_sum = top_relative_abun.cumsum()
        top_relative_abun.index = top_df_abun.index.values.tolist()
        abun = top_relative_abun.copy()
        abun.to_csv("output/plots/" + prefix + "_top20_relative_counts.csv")




# =============================================================================
# Plot
# =============================================================================

prev = ""
previous = pd.Series(dtype=float)

fig = plt.figure(1)

for i in range(0, len(abun.index.values.tolist())): #for each cluster i.e list of abun index
    if not prev:
        plt.bar(abun.keys(), abun.iloc[i, :], label=abun.index.values.tolist()[i], width=0.9)
    else:
        plt.bar(abun.keys(), abun.iloc[i, :], bottom=previous, label=abun.index.values.tolist()[i],  width=0.9)
    prev = 'y'
    if i != 0:
        previous = abun.iloc[i, :] + abun_sum.iloc[i-1, :]
    else:
        previous = abun.iloc[i, :]

plt.xticks( rotation='vertical', fontsize = 4, verticalalignment='center_baseline')

plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1, frameon=False, fontsize = 5)

if 'a' in plot:
    plt.margins(x = 0.01, y=0.05)
if 'r' in plot:
    plt.margins(0)
if 'y' in top20:
    plot = plot + '_top20'

#fig.savefig(str("plots/" + prefix + '_' + plot + '_abun_plot.png'), bbox_inches='tight', dpi = 400)
plt.show()

if 'y' in top20:
    colours = ['#502db3', '#008080', '#c200f2', '#f2c200','#36a3d9', '#e6beff','#8c0025', '#f58231', '#bf0080', '#cad900','#911eb4','#e5001f','#0066bf', '#000075','#338000', '#f032e6','#1bca00','#1d4010','#9a6324','#a9a9a9']

    abun_flip = abun.transpose()
    abun_flip.plot.bar(stacked=True, legend = None, figsize=(15,10), color=colours, width=0.9)
    plt.legend(loc='center left', labelspacing=-2.5, bbox_to_anchor=(1.0, 0.5), frameon=False)
    if 'r' in plot:
        plt.ylim(0,100)
    plt.tight_layout()
    plt.savefig('output/plots/' + prefix +'_' + plot +'_'+ 'abun.png', bbox_inches='tight', dpi = 300)
