#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 10:55:14 2020

@author: ac1513
"""

import pandas as pd
import sys
import numpy as np

con_read_file = sys.argv[1]
read_len_file = sys.argv[2]
contigs_file = sys.argv[3]
file_out = sys.argv[4]

norm_length = 100

con_read = pd.read_csv(con_read_file, sep='\t', names=["read", "contigs"], dtype = {'read':object, 'contigs':object})
lengths = pd.read_csv(read_len_file, sep='\t', names=["read", "length"])
lengths = lengths.drop_duplicates(subset=['read'])
con_read = con_read.merge(lengths, 'left')
con_read['length'] = con_read['length'].apply(lambda x: x/norm_length).astype(np.float64)
contigs = pd.read_csv(contigs_file, sep = '\t', header = None, names = ['contigs', 'lengths'])

norm_counts = con_read.groupby('contigs',sort=False)[['length']].agg(['sum', 'count'])
norm_counts = norm_counts['length'].apply(lambda x: round(x))
norm_counts['sum'] = norm_counts['sum'].astype('int64')
norm_counts = norm_counts.reset_index()
norm_counts = norm_counts[['contigs', 'count', 'sum']]
merged = pd.merge(contigs, norm_counts, on='contigs', how ='left', sort=False)
merged = merged.fillna(value=0)
merged['count'] = merged['count'].astype(np.int64)
merged['sum'] = merged['sum'].astype(np.int64)
merged.columns = ["contigs", "length", "raw", "normalised"]

#Outputs file in format -> contig, count, length normalised count
merged.to_csv(file_out, sep='\t', index = False, header = True)
