#!/usr/bin/python

# feeder3.py varient
# stand-alone version
# modified to make smaller files and hopefully run faster
# using command line arguments to calculate range of values
# so it can be run as an array job
# this version uses pre-computed values and _hopefully_
# really makes the correct comparisons!


import sys, getopt
import csv
import numpy as np


def inputter(argv):
	incoming = ''
	outgoing = ''
	starting = 0
	ending = 0
	try:
		opts, args = getopt.getopt(argv,"hi:o:s:e:",["ifile=","ofile=","st=","en="])
	except getopt.GetoptError:
		print('bin_feeder.py -i <input_file> -o <output_file> -s <start_range> -e <end_range>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('bin_feeder.py -i <input_file> -o <output_file> -s <start_range> -e <end_range>')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			incoming = arg
		elif opt in ("-o", "--ofile"):
			outgoing = arg
		elif opt in ("-s", "--start"):
			starting = arg
		elif opt in ("-e", "--end"):
			ending = arg
	return(incoming, outgoing, starting, ending)


def pcc(x, y):
	prod = np.sum(np.multiply(x,y))
	divx = np.sqrt(np.sum(np.square(x)))
	divy = np.sqrt(np.sum(np.square(y)))
	result = prod/(divx*divy)
	return result


variables = inputter(sys.argv[1:])
test_file = variables[0]
write_file = variables[1]
sstart = int(variables[2])
eend = int(variables[3])

#print('Loading contig names')
names = np.load('CRISPR_names.out.npy')

#print('Loading values')
b = np.load('CRISPR_b.out.npy')

#print('Loading differences')
diffs = np.load('CRISPR_diffs.out.npy')

#print('Preparing arrays')
name_1 = np.asarray(names[sstart:eend])
name_2 = np.asarray(names[sstart:])
b_1 = b[sstart:eend-1]
b_2 = b[sstart+1:]
diff_1 = diffs[sstart+1:eend-1]
diff_2 = diffs[sstart+2:]

total_1 = np.array((name_1, b_1, diff_1))
total_2 = np.array((name_2, b_2, diff_2))

with open(write_file, 'w') as sender:
	line_1 = 0
	for ploop1 in total_1[2]:
		contig_x = total_1[0][line_1]
		line_2 = line_1
		for ploop2 in total_2[2][line_2:]:
			contig_y = total_2[0][line_2+1]
			resp_val = pcc(ploop1, ploop2)
			if resp_val >= 0.99:
				line_out = (contig_x, contig_y, resp_val)
				writer = csv.writer(sender)
				writer.writerow(line_out)
			line_2 += 1
		line_1 += 1
