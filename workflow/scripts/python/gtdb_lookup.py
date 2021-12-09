#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 15:53:02 2021

Script to convert GTDB names to closest LCA in NCBI taxonomy. Using GTDB metadata tsv file.

@author: ac1513
"""

import pandas as pd
import json
import argparse

# =============================================================================
# Iterates throungh tsv output going up the tree from min taxa level
# =============================================================================
def search_taxa(match, cut_taxo):
    ncbi_taxo = {"Species" : "s__", "Genus" : "g__", "Family" : "f__", "Order": "o__", "Class": "c__", "Phylum":"p__", "Domain":"d__" }
    temp_set = set()
    for taxa in cut_taxo:
        taxa_search = ncbi_taxo[taxa]
        for row in match:
            ncbi_match = row.split(";")
            for item in ncbi_match:
                if (taxa_search in item) and (len(item) > 3):
                    temp_set.add(item)
        if len(temp_set) > 0:
            break
    return [*temp_set,]

# =============================================================================
# Command line parsing
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("gtdb_file", help="kraken top file with GTDB taxonomy")
parser.add_argument("taxo_level", help="min taxo level to run at (G, S etc)")
parser.add_argument("bac_lookup", help="metadata tsv file ")
parser.add_argument("arc_lookup", help="metadata tsv file ")
parser.add_argument("output_json", help="name of json file for output e.g. .json")


args = parser.parse_args()

gtdb_file = args.gtdb_file
min_taxo= args.taxo_level
bac_lookup_file = args.bac_lookup
arc_lookup_file = args.arc_lookup
out_file = args.output_json


# =============================================================================
# Read input, most handled by parser now
# =============================================================================
#min_taxo_level = "Genus" #change this
#file = "YA_B_oct_alls1p3_gtdb_999_G_top_kraken.out"

in_list = pd.read_csv(gtdb_file, sep = "\t", names = ("Per", "1", "2", "3", "4", "ID"))
input_gtdb = in_list["ID"].str.strip().dropna().drop_duplicates().to_list()
# input_gtdb = ["UBA2224", "T78", "Methanothrix", "Sedimentibacter", "UBA1413", "Streptococcus", "UBA5389", "UBA668", "Microbacterium"]

# =============================================================================
# Reading in files and taxonomy
# =============================================================================

all_taxo = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain"]
short_taxo = {"S" : "Species", "G" : "Genus", "F" : "Family", "O": "Order", "C": "Class", "P":"Phylum", "D":"Domain" }
min_taxo_level = short_taxo[min_taxo]

cut_taxo = all_taxo[all_taxo.index(min_taxo_level):]

gtdb_lookup_bac = pd.read_csv(bac_lookup_file, sep = "\t", header = 0, low_memory=False)
gtdb_lookup_bac = gtdb_lookup_bac.loc[:, ["gtdb_taxonomy", "ncbi_taxonomy"]] #remove 108 extra columns...

gtdb_lookup_arc = pd.read_csv(arc_lookup_file, sep = "\t", header = 0)
gtdb_lookup_arc = gtdb_lookup_arc.loc[:, ["gtdb_taxonomy", "ncbi_taxonomy"]] #remove 108 extra columns...

# =============================================================================
# Searching for matches
# =============================================================================
lookup_dict = {}

for name in input_gtdb.copy():
    match = gtdb_lookup_bac[gtdb_lookup_bac["gtdb_taxonomy"].str.contains(name)]["ncbi_taxonomy"].to_list() #search Bac at minimum level
    if match == []:
        match = gtdb_lookup_arc[gtdb_lookup_arc["gtdb_taxonomy"].str.contains(name)]["ncbi_taxonomy"].to_list()
    if match != []:
        x = search_taxa(match, cut_taxo)
        lookup_dict[name] = x
        input_gtdb.remove(name)

with open(out_file, 'w', encoding = 'utf-8') as f:
    json.dump(lookup_dict, f, sort_keys = True, indent = 2)

print("Not classified: ", input_gtdb)
