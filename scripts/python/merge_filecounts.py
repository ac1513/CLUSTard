#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 11:28:00 2020

@author: ac1513
"""

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('method', help='whether use normalised or raw counts (norm/raw)', type=str)
parser.add_argument('out', help='output file + extension', type=str)
parser.add_argument('-l', '--counts-list', dest='counts', nargs='+', default=[])
args = parser.parse_args()

counts = args.counts
method = args.method
out = args.out

if method == "norm":
    column = "normalised"
    to_drop = "raw"
elif method == "raw":
    column = "raw"
    to_drop = "normalised"
else:
    print("Error `", method, "` not valid. Please specify either `norm` or `raw`. Defaulting to normalised.")
    column = "normalised"
    to_drop = "raw"

base_df = pd.read_csv(counts[0], sep = '\t')
base_df = base_df.set_index("contigs")
initial_file = str(counts[0].replace('/','-'))
base_df = base_df.drop(to_drop, axis = 1)
base_df.columns = ['length', initial_file]

for file in counts[1:]:
    new_df = pd.read_csv(file, sep = '\t')
    new_df = new_df.set_index('contigs')
    col = str(file.replace('/','-'))
    base_df[col] = new_df[column]

base_df.to_csv(out, sep='\t', header=False)
