#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:52:15 2019

@author: ac1513
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import pandas as pd
import math

def bins(file, bins):
    data = open(file, 'r')
    counts = []
    for i in range(len(bins)):
        counts.append(0)
    x = data.readline()
    x = data.readline()
    while x:
       a = x.split()
       pos = 0
       for loop in bins:
          if int(a[1]) < loop*100000:
             counts[pos] += 1
             break
          else:
             pos += 1
       x = data.readline()
    data.close()
    gt = 0
    for i in counts:
       gt += i
    return counts

# =============================================================================
# Command line parsing
# =============================================================================
parser = argparse.ArgumentParser(description='usage = python entrez_down.py file_list_of_queries')
parser.add_argument('unbinned_file', help='output from checkm unbinned', type=str)
parser.add_argument('binned_file', help='output from bash script', type=str)
parser.add_argument('prefix', help='prefix of the jobs', type=str)

args = parser.parse_args()

unbinned_file = args.unbinned_file
binned_file = args.binned_file
prefix = args.prefix

#unbinned_file = "unbinned_stats.tsv"
#binned_file = "sorted_lengths.tsv"

unbin_df = pd.read_csv(unbinned_file, sep = '\t')
bin_df = pd.read_csv(binned_file, sep = '\t', names = ["Contig", "Length"])
max_len = max(math.ceil(unbin_df["Length"].max()/100000), math.ceil(bin_df["Length"].max()/100000))

groups = [0.02,0.05,0.1,0.2,0.5]
for i in range(1, max_len+1):
    groups.append(i)

unbinned = bins(unbinned_file, groups)
binned = bins(binned_file, groups)

x = []

for i in range(len(groups)):
    x.append(i)

bin_bars = []
unbin_bars = []

for i in range(len(unbinned)):
    total = unbinned[i] + binned[i]
    print(total, unbinned[i], binned[i])
    if total > 0:
        if binned[i] > 0 :
            bin_bars.append((binned[i]/total)*100)
        else:
            bin_bars.append(0)
        if unbinned[i] >0:
            unbin_bars.append((unbinned[i]/total)*100)
        else:
            unbin_bars.append(0)
    else:
        bin_bars.append(0)
        unbin_bars.append(0)

fig = plt.figure(1)
plt.bar(x, bin_bars, color='forestgreen', edgecolor='white', label='binned', width=0.9)
plt.bar(x, unbin_bars, bottom=bin_bars, color='purple', edgecolor='white', label='unbinned', width=0.9)
plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1, frameon=False)
plt.xlabel('Size (10Kb)', size=8)
plt.ylabel('% Contigs', size=8)

plt.xticks(x, groups)
plt.tick_params(labelsize = 7)
#plt.show()
fig.savefig(str('output/plots/' + prefix + '_bin_contigs.png'), bbox_inches='tight', dpi = 400)
