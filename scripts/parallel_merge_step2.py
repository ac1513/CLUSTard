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

short_list = []
nr_list = []
final_list = []

def inputter(argv):
	in_file = ''
	out_file = ''
	file_number = 0
	try:
		opts, args = getopt.getopt(argv,"hi:o:n:",["ifile=","ofile=", "nfile="])
	except getopt.GetoptError:
		print('parallel_merge_step2.py -i <input_file_root> -o <output_file> -n <number_of_files')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('new_step2_parallel.py')
			print('     compiles multiple cluster files into a single set')
			print('     -i/--ifile     root name of input files (no number post-script)')
			print('     -o/--ofile     output file name')
			print('     -n/--num       number of files to be compiled')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			in_file = arg
		elif opt in ("-o", "--ofile"):
			out_file = arg
		elif opt in ("-n", "--num"):
			file_number = arg
	return(in_file, out_file, file_number)

variables = inputter(sys.argv[1:])
read_file = variables[0]
write_file = variables[1]
number_files = int(variables[2])

# open first file
first_file = read_file+'.1'

with open (first_file, 'r') as master:
#	print('opening '+str(master))
	master_list = json.load(master)

# open sequential files and merge into the master list if they match
for file_cycle in range(1,number_files):
	current = read_file+'.'+str(file_cycle+1)
	with open (current, 'r') as working_file:
#		print('opening '+str(working_file))
		working_list = json.load(working_file)
		for row_1 in master_list:
			x = set(row_1)
			for row_2 in working_list:
				if x & set(row_2):
#					print('match between master and '+str(working_file))
					x.update(row_2)
					working_list.remove(row_2)

# add any sets left to the end of the master list
	for entries in working_list:
#		print('adding unmatched entries from '+str(working_file))
		y = [list(set(entries))]
		master_list.extend(y)

with open(write_file, 'w') as outgoing:
	json.dump(master_list, outgoing)
