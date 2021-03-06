#!/usr/bin/env python3

# 15/2/17 see if changing the logical flow of operations speeds things up
# make incoming file into a list
# make a set from list that matches Pcc >= 0.997
# then make a list of sets by exclusion instead of inclusion (I think)
# 18/2/17 modified to run from shell script in parallel
# from small files
# 19/2/17 now changed further to merge the resulting files

import sys, getopt
import csv
import json
import argparse


def set_default(obj):
    if isinstance(obj, set):
        return list(obj)
    raise TypeError

short_list = []
nr_list = []
final_list = []

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input-list', dest='input', nargs='+', default=[])
parser.add_argument('-o', dest='output', help='location of output', type=str)
args = parser.parse_args()
list_files = args.input
write_file = args.output

# open first file
first_file = list_files[0]

with open (first_file, 'r') as master:
#	print('opening '+str(master))
	master_list = json.load(master)
	l = master_list
	out = []
	while len(l)>0:
		first, *rest = l
		first = set(first)
		lf = -1
		while len(first)>lf:
			lf = len(first)
			rest2 = []
			for r in rest:
				if len(first.intersection(set(r)))>0:
					first |= set(r)
				else:
					rest2.append(r)
				rest = rest2
			out.append(first)
		l = rest
#	print(master_list)
#	print("has changed to:")
#	print(out)
	master_list = out

# open sequential files and merge into the master list if they match
for current_f in list_files[1:]:
	with open (current_f, 'r') as working_file:
		working_list = json.load(working_file)
		for row_1 in master_list:
			x = set(row_1)
			for row_2 in working_list:
				if x & set(row_2):
#					print('match between master and '+str(working_file))
#					print(x, row_2)
					x.update(row_2)
#					print(x)
					working_list.remove(row_2)
	for entries in working_list: # add any sets left to the end of the master list
#		print('adding unmatched entries from '+str(working_file))
		y = [list(set(entries))]
		master_list.extend(y)

l = master_list
out = []
while len(l)>0:
          first, *rest = l
          first = set(first)
          lf = -1
          while len(first)>lf:
                        lf = len(first)
                        rest2 = []
                        for r in rest:
                                if len(first.intersection(set(r)))>0:
                                        first |= set(r)
                                else:
                                        rest2.append(r)
                                rest = rest2
                        out.append(first)
          l = rest
#       print(master_list)
#       print("has changed to:")
#       print(out)
master_list = out

with open(write_file, 'w') as outgoing:
	json.dump(master_list, outgoing, default=set_default)
