import pickle
import glob
import re
import os
import copy
import subprocess
import timeit
from datetime import datetime, timedelta
from collections import defaultdict, Counter
from functools import reduce

import vcf
import pysam
import numpy as np
import pandas as pd
import dask.dataframe as dd
from Bio import SeqIO
from Bio.SeqUtils import seq1

from utils.file_loaders import load_metadata, load_vdb_mutation_data, merge_data, load_vdb_df, load_alias_data
from utils.quaid_func import get_counts, get_all_voc, get_recent_nt_df
from utils.quarc_func import build_count_df

import argparse

def check_dataframe_status(df):
    print(f"Shape: {df.shape}")
    print(f"Features: {list(df.columns.values)}")
    
def get_prevalence_table(vdb_df_path, genome_counts, exclude_recombinant):
    # load and build weekly vdb snp count dataframe
    vdb_df = load_vdb_df(vdb_df_path)
    count_df = build_count_df(genome_counts, vdb_df, exclude_recombinant)
    # get latest date
    date_max = count_df.index.get_level_values('Collection date').max()
    # generate start_to_date weekly prevalence rate dataframe 
    weekly_data = []
    for idx, row in count_df.iterrows():
        date = idx[0]
        lineage = idx[1]
        for mut in row['SNPs']:
            weekly_data.append((date, lineage, mut, 
                                     row['SNPs'][mut],
                                     row['Total count']))
    weekly_nt_df = pd.DataFrame(data=weekly_data, columns=['Date', 'Lineage', 'NT_mut', 'Prevalence', 'Occurence count'])
    weekly_nt_df['Prevalence rate'] = weekly_nt_df['Prevalence']/weekly_nt_df['Occurence count']
    prevalence_rate_table = weekly_nt_df.pivot(index=["Lineage", "Date"], columns="NT_mut", values="Prevalence rate").sort_index()
    prevalence_rate_table.index = [row[0] + "/"+ row[1].strftime('%Y-%m-%d') for row in prevalence_rate_table.index.values]
    
    return prevalence_rate_table.T, count_df, date_max

def get_sublineage_mutations(snp_prevalence_df, min_prevalence=0.5):
    
    snp_series = snp_prevalence_df.apply(lambda row: set(row[row > min_prevalence].index.values), axis=1)
    
    sublineage_mutations_lookup = dict()
    for index, value in snp_series.items():
        if index[-1] != '-' and value:
            sublineage_mutations_lookup[index] = value
                    
    return sublineage_mutations_lookup

def build_genome2snp_database(merged_data, quarc_db_dir):
    '''
    returns a dictionary of dictionary
    level 0:
        key: lineage/date where date is the start of the week
        value: a set of level 1 dictionary
    level 1:
        key: Accession ID of the genome 
        value: a set of SNPs that the genome contains based on VDB
    the dictionary is saved as a pkl file for future usage 
    '''
    quarc_db_df = merged_data.copy()
    quarc_db_df["Week Start"] = quarc_db_df['Collection date'].apply(lambda x: (x - timedelta(days=x.dayofweek-7)).strftime('%Y-%m-%d'))
    quarc_db_df['Lineage_Date'] = quarc_db_df['Pango lineage'] + "/" + quarc_db_df['Week Start']
    quarc_db_df.drop(['Collection date', 'Pango lineage', 'Week Start'], axis=1, inplace=True)
    quarc_db_df = quarc_db_df.reset_index().set_index(['Lineage_Date', 'Accession ID']).sort_index()
    quarc_db_df['SNPs'] = quarc_db_df.apply(lambda row: set(row.SNPs.strip().split(" ")), axis=1)
    quarc_db = quarc_db_df.groupby(level=0).apply(lambda quarc_db_df: quarc_db_df.xs(quarc_db_df.name)['SNPs'].to_dict()).to_dict()
        
    # create a binary pickle file 
    f = open(os.path.join(quarc_db_dir, "quarc_db.pkl"),"wb")
    # write the python object (dict) to pickle file
    pickle.dump(quarc_db, f)
    # close file
    f.close()
    
    return quarc_db

if __name__ == "__main__":    
    '''main function'''
    # Tested command:
    # python crykey_build_db.py -m /home/Users/ns58/MSA-01102023/metadata.tsv -l /home/Users/ns58/MSA-01102023/vdb_output/vdb_lineage_df_week.csv -v /home/Users/ns58/MSA-01102023/vdb_01102023_trimmed_nucl.txt -d 2023-01-10
    
    parser = argparse.ArgumentParser(description="Building Crykey Database based on GISAID MSA.")
    parser.add_argument("-m", "--metadata", type=str, required=True, help="Path to the metadata of the GISAID MSA in tsv format.")
    parser.add_argument("-l", "--lineages", type=str, required=True, help="Path to the vdb lineage dataframe.")
    parser.add_argument("-v", "--vdb", type=str, required=True, help="Path to the trimmed vdb file.")
    parser.add_argument("-d", "--date", type=str, required=True, help="Date of the lastest record in the GISAID MSA, with the format [YYYY-MM-DD].")
    parser.add_argument("-o", "--output", type=str, required=False, help="Output path of the Crykey Database. Default:[crykey_dbs]", default="crykey_dbs")
    parser.add_argument('--exclude-recombinant', dest='norecomb', action='store_true',
                        help="Exclude recombinant lineages. Default:[False]")
    
    parser.set_defaults(norecomb=False)
    
    args = parser.parse_args()
    metadata_f = args.metadata
    vdb_df_path = args.lineages
    trimmed_vdb_f = args.vdb
    max_date = args.date
    quarc_db_dir = args.output
    exclude_recombinant = args.norecomb
    
    # metadata_f = "/home/Users/ns58/MSA-01102023/metadata.tsv"
    # vdb_df_path = "/home/Users/ns58/MSA-01102023/vdb_output/vdb_lineage_df_week.csv"
    # trimmed_vdb_f = "/home/Users/ns58/MSA-01102023/vdb_01102023_trimmed_nucl.txt"
    # max_date = "2023-01-10"
    # quarc_db_dir = 'quarc_dbs_01102023_recombinant'
    # exclude_recombinant = False

    if not os.path.exists(quarc_db_dir):
        os.mkdir(quarc_db_dir)

    metadata = load_metadata(metadata_f)
    genome_counts = get_counts(metadata)
    metadata.reset_index(inplace=True)
    metadata.set_index("Accession ID", inplace=True)
    print("Metadata loading completed.")

    mutations_data = load_vdb_mutation_data(trimmed_vdb_f)
    mutations_data.drop(["Date", "Region"], axis=1, inplace=True)
    merged_data = merge_data(mutations_data, metadata)
    merged_data.set_index("Accession ID", inplace=True)
    genome2snp_dict = build_genome2snp_database(merged_data, quarc_db_dir)
    print("Done building Genome-to-SNP Database.")
    
    prevalence_rate_table, count_df, date_max = get_prevalence_table(vdb_df_path, genome_counts, exclude_recombinant)

    sublineage_mutations_lookup = get_sublineage_mutations(prevalence_rate_table, min_prevalence=0.5)
    f = open(os.path.join(quarc_db_dir, "mutation50_lookup.pkl"), "wb")
    pickle.dump(sublineage_mutations_lookup, f)
    f.close()

    sublineage_mutations_lookup = get_sublineage_mutations(prevalence_rate_table, min_prevalence=0)
    f = open(os.path.join(quarc_db_dir, "mutation0_lookup.pkl"), "wb")
    pickle.dump(sublineage_mutations_lookup, f)
    f.close()
    
    print("Done building SNP lookup tables.")
