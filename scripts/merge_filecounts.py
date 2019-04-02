# updated version using glob and other tricks to be pythonic rather than brutal

import os
import glob

def firstFile(fileName):
        bits = fileName.readline().rsplit('\t',1)
        return (bits[0])


def otherFile(fileName):
        bits = fileName.readline().split('\t')
        return (bits[2])


path = os.path.dirname(os.path.realpath(__file__))
dir = '/short_reads' # assumes files are here

os.chdir(path)

fname = sorted(glob.glob('*.txt'))

print ('files with data to be merged: '+str(fname))
# count lines in one file (they should all be the same...)

count = 0
f = open(fname[0], 'r')
for line in f:
        count +=1
f.close()

# count is now set to the length of the file

print ('number of entries per file: '+str(count))

# now need to read in all files

filedata = [open(file_name, 'r') for file_name in fname]
oo = open('read_counts.output', 'w')

for runthrough in range (0,count):
        start_marker = filedata[0]
        end_marker = filedata[-1]
        string_text = firstFile(start_marker)+'\t'
        for line_out in filedata[1:-1]:
                string_text = string_text + otherFile(line_out)+'\t'
        string_text = string_text + otherFile(end_marker)+'\n'
        oo.write(string_text)

for closer in filedata:
        closer.close()

print ('merged data saved in: read_counts.output')
