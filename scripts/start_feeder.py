#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 13:50:48 2019

@author: ac1513
"""

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('loc', help='location count files are in', type=str)
parser.add_argument('jobid', help='location count files are in', type=str)
args = parser.parse_args()
loc = args.loc
jobid = args.jobid

df = pd.read_csv(loc + '/' + jobid + '_read_counts_derived.csv', header=None, index_col = 0)

names = df.index.values.tolist()

col_num = len(df.columns)

df2 = df #seem to need to have a copy of the df to calc mean

df = df.drop(df.columns[len(df.columns)-1], axis=1) # drop last column so don't include it in stats - is still in df2..

df2["mean"] = df.mean(axis=1)
df2["sd"] = df.std(ddof = 1, axis=1)

diffs = df.sub(df.mean(axis=1), axis=0)

diffs.to_csv(loc + '/' + jobid + '_diffs.csv') #diffs
df2.to_csv(loc + '/' + jobid + '_values.csv') #diffs
