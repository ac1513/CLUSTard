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
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

# =============================================================================
# Command line parsing
# =============================================================================
parser = argparse.ArgumentParser(description='usage = python entrez_down.py file_list_of_queries')
parser.add_argument('in_file', help='the name of the file containing a list of csv files', type=str)
parser.add_argument('prefix', help='prefix of the jobs', type=str)
parser.add_argument('-k', '--kraken', help = 'merged kraken input file', type=str)
parser.add_argument('-k_l', '--kraken_level', help = 'merged kraken input file', type=str)
parser.add_argument('-cm', '--checkm_file', help = 'checkm output file - in tab format', type = str)
parser.add_argument('samples', help='samples.tsv file', type=str)
parser.add_argument('dates', help='plot date scale y/n', type=str)
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

if args.checkm_file:
    checkm_file = args.checkm_file
else:
    print("Include checkm output if want comp/cont info on the plot")

# =============================================================================
# Plot global settings
# =============================================================================

matplotlib.rcParams['lines.linewidth'] = 0.5
matplotlib.rcParams['ytick.left'] = True
matplotlib.rcParams['ytick.minor.size'] = 1
matplotlib.rcParams['ytick.minor.width'] = 0.25
matplotlib.rcParams['axes.linewidth'] = 0.5
colours = ["crimson", "purple", "tab:cyan", "seagreen", "darkorange", "tab:pink", "darkslateblue"]

# =============================================================================
# Read in input file(s)
# =============================================================================

set_groups = set()
if "y" in args.dates.lower():
    df_samples = pd.read_csv(samples, sep ='\t', parse_dates = ["date"])
else:
    df_samples = pd.read_csv(samples, sep ='\t')

groups = df_samples["group"].tolist() #get rid of header

for item in groups:
    set_groups.add(item)
dc = {}
for item in list(set_groups):
    dc[item] = -1
for i in groups:
    dc[i] +=1

with open(in_file, 'r') as text_file:
    files = text_file.read().strip().split()

# =============================================================================
# Plot
# =============================================================================

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
            x_prev_start = 0
            x_prev_end = 0
            for item in list(set_groups):
                x_end = x_start + dc[item] + 1
                y_mean = df.mean()[x_start:x_end]
                top = df.max()[x_start:x_end]
                bottom = df.min()[x_start:x_end]

                if "y" in args.dates.lower():
                    df_samples["date"] = pd.to_datetime(df_samples["date"])
                    if x_start == 0:
                        x_data = df_samples["date"][x_start:x_end]
                    else:
                        x_data = df_samples["date"][x_start:x_end] + (x_data[-1:][x_start-1] - df_samples["date"][0])
                else:
                    x_data = range(x_start, x_end)
                plt.plot(x_data, y_mean, color=colours[item])
                plt.fill_between(x_data, top, bottom, facecolor='gray', alpha=0.5)
                x_prev_start = x_start
                x_prev_end = x_end
                x_start = x_end

            plt.tick_params(labelbottom=False)

            plt.semilogy()

            x1,x2,y1,y2 = plt.axis()
            plt.axis((x1,x2,0.0001,100))
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
                        if "y" in args.dates.lower():
                            plt.text(df_samples["date"][1], 0.7, per+'%:  ' + cont, fontsize=2)
                        else:
                            plt.text(0.5, 0.7, per+'%:  ' + cont, fontsize=2)
            if checkm_file:
                checkm_df = pd.read_csv(checkm_file, sep = '\t', index_col = 0)
                comp = checkm_df["Completeness"][file.split('/')[-1:][0].split('.')[0]]
                conta = checkm_df["Contamination"][file.split('/')[-1:][0].split('.')[0]]
                if "y" in args.dates.lower():
                    plt.text(df_samples["date"][1], 1.7, str(comp)+'%: Complete ' + str(conta)+'%: Contamination', fontsize=2)
                else:
                    plt.text(0.5, 1.7, str(comp) + '%:  Complete ' + str(conta) + '%:  Contamination', fontsize=2)

            if "y" in args.dates.lower():
                plt.text(df_samples["date"][1], 40, na, fontsize = 2, fontweight='bold')
                plt.text(df_samples["date"][1], 9, nu+' cov:'+av_cov+'+/-'+sd_cov + ', ' + tot_len +'kb', fontsize=2)
                plt.text(df_samples["date"][1], 4, 'GC% '+ av_gc +'+/-'+ sd_gc, fontsize=2)
            else:
                plt.text(0.5, 40, na, fontsize = 2, fontweight='bold')
                plt.text(0.5, 9, nu +', cov:'+av_cov+'+/-'+sd_cov + ', ' + tot_len +'kb', fontsize=2)
                plt.text(0.5, 4, 'GC% '+ av_gc +'+/-'+ sd_gc, fontsize=2)
            plt.tick_params(axis='x', labelsize=2, pad=0, direction='out', length=1, width=0.25)
            plt.tick_params(axis = 'y', labelsize=2, pad=0, direction='out', length=1)
            plt.tick_params(right=False, top=False)
            gsplace += 1
#    plt.show()
    plt.savefig(str(counter) + '_' + prefix + '_plot.pdf', type='pdf', dpi=300)
    print('Generated plot number ' + str(counter) + ' -> ' + str(counter) + '_' + prefix + '_plot.pdf')
    plt.close('all')
