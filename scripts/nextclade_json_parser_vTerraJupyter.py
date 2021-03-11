#! /usr/bin/env python

import os
import glob
import re
import shutil
import pandas as pd 
from datetime import date
import sys
import json

# note before you can use this script from the command line you must make the script executable
# locate this script in ~/scripts and be sure that ~/scripts is ammended to your $PATH variable
# to use nexclade_json_parser.py <path to json file>
# will create two output files:
# nextclade_variant_summary.csv
# nextclade_results.csv


def extract_variant_list(json_path, wd):
    # get date
#     today_date = str(date.today())

    # create pd data frame to fill
    df = pd.DataFrame()
    accession_id_list = []
    mutation_list = []
    gene_list = []
    refAA_list = []
    altAA_list = []
    codon_pos_list = []


    with open(json_path) as f:
        data = json.load(f)

    for i in range(len(data)):

    #     print(data[i]['seqName'])
        if 'aaDeletions' in data[i].keys():
            aa_deletions = data[i]['aaDeletions']
            for item in aa_deletions:
                gene=item['gene']
                refAA= item['refAA']
                altAA= 'del'
                pos=item['codon'] + 1

                mutation = '%s_%s%d%s' % (gene, refAA, pos, altAA)
                accession_id_list.append(data[i]['seqName'])
                mutation_list.append(mutation)
                gene_list.append(gene)
                refAA_list.append(refAA)
                altAA_list.append(altAA)
                codon_pos_list.append(pos)


        if 'aaSubstitutions' in data[i].keys():
            aa_subs = data[i]['aaSubstitutions']
            for item in aa_subs:
                gene = item['gene']
                refAA = item['refAA']
                if item['queryAA'] == '*':
                    altAA = 'stop'
                else:
                    altAA = item['queryAA']     
                pos = item['codon'] + 1

                mutation = '%s_%s%d%s' % (gene, refAA, pos, altAA)
                accession_id_list.append(data[i]['seqName'])
                mutation_list.append(mutation)
                gene_list.append(gene)
                refAA_list.append(refAA)
                altAA_list.append(altAA)
                codon_pos_list.append(pos)


    df['accession_id'] = accession_id_list
    df['variant_name'] = mutation_list
    df['gene'] = gene_list
    df['codon_position'] = codon_pos_list
    df['refAA'] = refAA_list
    df['altAA'] = altAA_list

    path = os.path.join(wd, 'nextclade_variant_summary.csv') 
    df.to_csv(path, index=False)
    
    
def get_nextclade(json_path, wd):
    # get date
#     today_date = str(date.today())

    # create pd data frame to fill
    accession_id_list = []
    clade_list = []
    df = pd.DataFrame()

    with open(json_path) as f:
        data = json.load(f)

    for i in range(len(data)):
        if 'clade' in data[i].keys():
            accession_id_list.append(data[i]['seqName'])
            clade_list.append(data[i]['clade'])

    df['accession_id'] = accession_id_list
    df['nextclade'] = clade_list

    path = os.path.join(wd, 'nextclade_results.csv')
    df.to_csv(path, index = False)

    
if __name__ == '__main__':
    
    if len(sys.argv) > 1:
        json_path = sys.argv[1]
    else:
        print('specify path to nextclade json file')
        
    get_nextclade(json_path)
    extract_variant_list(json_path)
        