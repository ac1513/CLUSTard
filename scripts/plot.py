#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import statistics
import matplotlib
#matplotlib.use('pdf')
#%matplotlib inline
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import re


parser = argparse.ArgumentParser(description='usage = python entrez_down.py file_list_of_queries')
parser.add_argument('in_file', help='the name of the file containing a list of csv files', type=str)
parser.add_argument('prefix', help='prefix of the jobs', type=str)
parser.add_argument('-k', '--kraken', help = 'merged kraken input file', type=str)
parser.add_argument('-k_l', '--kraken_level', help = 'merged kraken input file', type=str)
parser.add_argument('samples', help='samples.tsv file', type=str)
args = parser.parse_args()
in_file = args.in_file
prefix = args.prefix
samples = args.samples

if args.kraken:
    kraken_file = args.kraken
else:
    print("Add a merged kraken input file if you want taxonomy info on the plot")

if args.kraken_level:
    prefix = str(prefix + "_" + args.kraken_level)
else:
    prefix = prefix

set_groups = set()
df_samples = pd.read_csv(samples, sep ='\t', index_col = 0)
groups = df_samples["group"].tolist() #get rid of header

for item in groups:
    set_groups.add(item)
dc = {}
for item in list(set_groups):
    dc[item] = -1
for i in groups:
    dc[i] +=1


matplotlib.rcParams['lines.linewidth'] = 0.5
matplotlib.rcParams['ytick.left'] = True
matplotlib.rcParams['ytick.minor.size'] = 1
matplotlib.rcParams['ytick.minor.width'] = 0.25
matplotlib.rcParams['axes.linewidth'] = 0.5
colours = ["crimson", "seagreen", "darkorange", "purple", "tab:pink", "darkslateblue", "tab:cyan"]

with open(in_file, 'r') as text_file:
    files = text_file.read().strip().split()

counter = 0
for i in range(0, len(files), 30):
    gs = gsp.GridSpec(5,6)
    gsplace = 0
    sub_files = files[i:i+30]
    counter += 1
    for file in sub_files:
        with open(file, 'r') as f:
            file = str(file)
            df = pd.read_csv(f, index_col='contig')
            av_cov = str('{0:.1f}'.format(statistics.mean(df['cover'].tolist())))
            sd_cov = str('{0:.1f}'.format(statistics.stdev(df['cover'].tolist())))
            tot_len = str('{0:.1f}'.format(sum(df['length'].tolist())/1000))
            av_gc =  str('{0:.1f}'.format(statistics.mean(df['GC'].tolist())))
            sd_gc = str('{0:.1f}'.format(statistics.stdev(df['GC'].tolist())))
            na = str(file.split('/')[-1:][0].split('.')[0][8:])
            nu = str(len(df))
            axes1 = plt.subplot(gs[gsplace])


            x_start = 0
            for item in list(set_groups):
                x_end = x_start + dc[item]
                x_range = range(x_start, x_end)
                y_mean = df.mean()[x_start:x_end]
                top = df.max()[x_start:x_end]
                bottom = df.min()[x_start:x_end]
                plt.plot(x_range, y_mean, color=colours[item])
                plt.fill_between(x_range, top, bottom, facecolor='gray', alpha=0.5)
                x_start = x_start + dc[item]+1

            plt.semilogy()

            x1,x2,y1,y2 = plt.axis()
            plt.axis((x1,x2,0.00001,10))


            plt.axhline(y=0.01, ls='--', lw = 0.25, c = 'black')

            if kraken_file:
                for line in open(kraken_file, 'r'):
                    if re.search(na, line):
                        cont = line.split('\t')[-1].strip()
                        if './' in cont:
                            cont = ' '
                        else:
                            per = line.split('\t')[1]
                            per = per.strip()
                        plt.text(0.5, 0.15, per+'%:  ' + cont, fontsize=2)
            plt.text(0.5, 4, na, fontsize = 2, fontweight='bold')
            plt.text(0.5, 1, nu+' cov:'+av_cov+'+/-'+sd_cov + ', ' + tot_len +'kb', fontsize=2)
            plt.text(0.5, 0.4, 'GC% '+ av_gc +'+/-'+ sd_gc, fontsize=2)
            plt.tick_params(axis='x', labelsize=2, pad=0, direction='out', length=1, width=0.25)
            plt.tick_params(axis = 'y', labelsize=2, pad=0, direction='out', length=1)
            plt.tick_params(right=False, top=False)
            gsplace += 1
#    plt.show()
    plt.savefig('output/plots/' + str(counter) + '_' + prefix + '_plot.pdf', type='pdf', dpi=300)
    print('Generated plot number ' + str(counter) + ' -> ' + str(counter) + '_' + prefix + '_plot.pdf')
    plt.close('all')
