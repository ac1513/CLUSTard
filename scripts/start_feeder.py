# start_feeder.py
# stand-alone version JC 13/06/17

import os
import csv
import numpy as np

path = os.path.dirname(os.path.realpath(__file__)) # path goes here
dir = '/in_files/' # directory containing sequence files goes here

test_file = '/read_counts_derived.csv'

print('Reading contig names')
names = []
with open(path+test_file, 'r') as f:
	for line in f.readlines():
		name, info = line.strip().split(',',1)
		names.append((name))

# determine number of columns in file
with open(path+test_file) as g:
	col_num = g.readline().count(',')

# set columns to be read
col_range = []
for loo in range(1, col_num, 1):
	col_range.append(loo)

print('Reading contig values from columns '+str(col_range))
values = np.genfromtxt(fname=path+test_file, delimiter=',', usecols=(col_range))

print('Calculating mean and standard deviation')
b = np.c_[values, np.mean(values, axis = 1), np.std(values, ddof=1, axis=1)]
diffs = np.array([col_range])

print('Calculating differences')
for i in b:  # add '[:499]' for testing
	a = np.array([[ ]])
	for v in i[:col_num-1]:
		a = np.append(a, v-i[col_num-1])
	c = np.array([a])
	diffs = np.append(diffs, c, axis=0)

print('Concatenating results')

np.save(path+'/CRISPR_names.out', names)
np.save(path+'/CRISPR_b.out', b)
np.save(path+'/CRISPR_diffs.out', diffs)
