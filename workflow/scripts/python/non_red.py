#!/usr/bin/env python3

# STEP 3 STARTS HERE:
# make a non-redundant list from the sets
# generated in step 2

import json
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('input', help='location of input', type=str)
parser.add_argument('output', help='location of output', type=str)
args = parser.parse_args()
input_file = args.input
output_file = args.output


final_list = []
with open(input_file ,'r') as in_file:
	master_list = json.load(in_file)

while True:
	try:
		test = master_list[0]
	except IndexError:
		print("No clusters identified - your Pcc threshold is probably too high. Decrease it in the config file, remove the clustering directory and try again.")
		exit()
	working_list = test
	for test_list in master_list[1:]:
		if not (set(test_list).intersection(test)):
			a = False
		else:
			a = True
		if a == False:
			next
		else:
			working_list.extend(test_list)
			master_list.remove(test_list)
	master_list.remove(test)
	x = set(working_list)
	x = list(x)
	final_list.append(x)
	if master_list == []:
		break

with open(output_file, 'w') as out_file:
	json.dump(final_list, out_file)
