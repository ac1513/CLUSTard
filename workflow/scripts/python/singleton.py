import glob
import csv
from Bio import SeqUtils as su
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('csv', help='read_counts_derived', type=str)
parser.add_argument('loc', help='location of singleton fasta files', type=str)
parser.add_argument('-l', '--header-list', dest='header', nargs='+', default=[])
args = parser.parse_args()

header = args.header
counts_file = args.csv
loc = args.loc

singleton = loc

header = ['contig'] + header + ['cover', 'length', 'GC']


with open(singleton) as fasta:
  name = fasta.readline()
  name = name.strip(">\n")
  faseq = fasta.readline()
with open(counts_file) as counts:
  for line in counts:
    if name in line:
            count_line = line.strip('\n').split(',')
            break
filename = str(singleton[:-6] + ".csv")
with open(filename, 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(header)
    count_line.append(len(faseq))
    count_line.append(su.GC(faseq))
    writer.writerow(count_line)
