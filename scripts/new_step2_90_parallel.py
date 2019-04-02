# make incoming file into a list
# make a set from list that matches Pcc >= 0.997
# then make a list of sets by exclusion instead of inclusion (I think)
# 18/2/17 modified to run from shell script in parallel

import sys, getopt
import csv
import json

short_list = []
nr_list = []
final_list = []

def inputter(argv):
	in_file = ''
	out_file = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print('new_step2_parallel.py -i <input_file> -o <output_file>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('new_step2_parallel.py -i <input_file> -o <output_file>')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			in_file = arg
		elif opt in ("-o", "--ofile"):
			out_file = arg
	return(in_file, out_file)

variables = inputter(sys.argv[1:])
read_file = variables[0]
write_file = variables[1]

with open (read_file, 'r') as incoming:
	file_reader = csv.reader(incoming, delimiter=',')
	for row in file_reader:
		if float(row[2]) >= 0.99:
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

