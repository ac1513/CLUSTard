#!/usr/bin/python

# feeder3.py variant
# stand-alone version
# modified to make smaller files and hopefully run faster
# using command line arguments to calculate range of values
# so it can be run as an array job
# this version uses pre-computed values and _hopefully_
# really makes the correct comparisons!

import sys, getopt
import csv
import numpy as np
import argparse
import pandas as pd

def pcc(x, y):
	prod = np.sum(np.multiply(x,y))
	divx = np.sqrt(np.sum(np.square(x)))
	divy = np.sqrt(np.sum(np.square(y)))
	result = prod/(divx*divy)
	return result

parser = argparse.ArgumentParser(description='')
parser.add_argument('cut_diffs', help='split diffs', type=str)
parser.add_argument('all_diffs', help='all diffs', type=str)
parser.add_argument('thresh', help='pr threshold', type=float)
parser.add_argument('output', help='output', type=str)
args = parser.parse_args()
cut_diffs = args.cut_diffs
diffs = args.all_diffs
write_file = args.output
thresh = args.thresh

df_diffs_all = pd.read_csv(diffs, header=None, index_col = 0)
df_diffs_cut = pd.read_csv(cut_diffs, header=None, index_col = 0)

with open(write_file, 'w') as sender:
	line_1 = 0
	for contig_x, row in df_diffs_cut.iterrows():
		row = row.to_numpy()
		line_2 = line_1+1
		for contig_y, row1 in df_diffs_all.iloc[line_2:].iterrows():
            row1 = row1.to_numpy()
			resp_val = pcc(row, row1)
			if resp_val >= thresh:
				line_out = (contig_x, contig_y, resp_val)
				writer = csv.writer(sender)
				writer.writerow(line_out)
			line_2 += 1
		line_1 += 1
