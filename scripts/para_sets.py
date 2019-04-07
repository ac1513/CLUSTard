# make incoming file into a list
# make a set from list that matches Pcc >= 0.997
# then make a list of sets by exclusion instead of inclusion (I think)
# 18/2/17 modified to run from shell script in parallel

import sys, getopt
import csv
import json
import argparse

short_list = []
nr_list = []
final_list = []

parser = argparse.ArgumentParser(description='')
parser.add_argument('input', help='location of input', type=str)
parser.add_argument('output', help='location of output', type=str)
parser.add_argument('thresh', help='location of output', type=float)
args = parser.parse_args()
read_file = args.input
write_file = args.output
thresh = args.thresh

with open (read_file, 'r') as incoming:
	file_reader = csv.reader(incoming, delimiter=',')
	for row in file_reader:
		if float(row[2]) >= thresh:
			short_list.append(row[:2])

while short_list:
	top = short_list[0]
	first = set(top)
	short_list.remove(top)
	for entry in short_list:
		if first & set(entry):
			first.update(entry)
			short_list.remove(entry)
	x = [list(set(first))]	# convert set to list to make compatible with json
	final_list.extend(x)

with open(write_file, 'w') as outgoing:
	json.dump(final_list, outgoing)
