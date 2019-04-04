# code to select a subset of the "csv" file generated
# this works!
# JC 23/06/16
# improved a bit concerning filenames  JC 13/06/17

# function that adds all the count values together to get a total
def summer(x):
    total = 0
    for loop in x:
        total = total + int(loop)
    return(total)

# function that divides each point by total counts for this contig
def deriver(x,y):
    answers=[]
    for loop in x:
        value = int(loop)/y
        answers.append(value)
    return(answers)

import csv as csv
import os
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('loc', help='location count files are in', type=str)
parser.add_argument('jobid', help='location count files are in', type=str)
parser.add_argument('thresh', help='location count files are in', type=str)
args = parser.parse_args()
loc = args.loc
jobid = args.jobid
thresh = int(args.thresh)

dir_name = str(loc + '/')
file_name = str(jobid + '_read_counts.output')

new_record=[[]]
nr=False

with open(dir_name+file_name, 'r') as data_store:
    line = csv.reader(data_store, delimiter='\t')
    for i in line:
        if int(i[1]) >= thresh:
            counts = summer(i[2:])
            values = deriver((i[2:]), counts)
            coverage = (counts*150)/int(i[1])
            values.append(coverage)	# add coverage to the end of the entry
            values.insert(0,i[0])	# add contig name to the front of the entry
            if nr:
                new_record.append(values)
            else:			# identifies first entry in the list
                new_record = [values]
                nr = True

with open(dir_name + '/' +  jobid + '_read_counts_derived.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(new_record)
