#!/usr/bin/env python3

from Bio import Entrez
import argparse

parser = argparse.ArgumentParser(description='usage = python entrez_down.py file_list_of_queries')
parser.add_argument('id_file', help='the name of the file containing a list of geneIDs', type=str)
parser.add_argument('family_file', help='the name of the file containing a list of families interested in each on a new line', type=str)
parser.add_argument('outgroup_file', help='the name of the file containing a list of outgroups interested in each on a new line', type=str)
args = parser.parse_args()
id_file = args.id_file
family_file = args.family_file
outgroup_file = args.outgroup_file

with open(family_file, "r") as family:
    families = []
    for line in family:
        families.append(line.strip())
with open(outgroup_file, "r") as outgroup:
    for line in outgroup:
        families.append(line.strip())


for item in families:
    ids_oi = []
    with open(id_file, "r") as genomeID:
        line = genomeID.readline().strip()
        while line.strip():
            if line.lower() == str(">" + item).lower():
                print(line)
                line = genomeID.readline().strip()
                while not ">" in line:
                    ids_oi.append(line)
                    line = genomeID.readline().strip()
                    if line.strip() == "":
                        break
            else:
                line = genomeID.readline().strip()
    if not ids_oi:
        print(">" + item)
        print("Family not found, add it to the genomes oi file or check your spelling\n")

    if ids_oi:
        for xid in ids_oi:
            search_term = xid.strip()
            Entrez.email = "ac1513@york.ac.uk"   # required by NCBI
            search_handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(search_handle)
            search_handle.close()

            gi_list = search_results["IdList"]
            count = int(search_results["Count"])
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            batch_size = 5    # download sequences in batches so NCBI doesn't time you out

            out_file = str(search_term + ".fasta")

            with open("genomes/"+out_file, "w") as out_handle:
                for start in range(0, count, batch_size):
                    end = min(count, start+batch_size)
                    print("Going to download record %s" % (xid))
                    fetch_handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text",retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
                    data = fetch_handle.read()
                    fetch_handle.close()
                    out_handle.write(data)

        print("Download completed\n")
