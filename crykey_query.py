import pickle
import glob
import re
import os
import copy
import subprocess
import timeit
import time
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
from utils.quaid_func import get_counts, get_all_voc, build_count_df, get_recent_nt_df

import argparse

def merge_dataframe(output_dir):
    dataframe_files = sorted(os.listdir(os.path.join(output_dir, 'cryptic_dataframe')))

    df_list = []
    for filename in dataframe_files:
        temp_df = pd.read_csv(os.path.join(output_dir, 'cryptic_dataframe', filename))
        df_list.append(temp_df)

    merged_df = pd.concat(df_list)
    merged_df.to_csv(os.path.join(output_dir, 'merged_df.csv'), index=False)

    print("Total Number of Unique Mutation Sets:", len(merged_df['Nt Mutations'].unique()))
    return merged_df

def mutation_comb_msa_query(test_set, sublineage_mutations_exist_lookup, genome2snp_dict, output_acc_id=False):
    possible_voc_sets = []
    for mutation in test_set:
        try:
            possible_voc_sets.append(sublineage_mutations_exist_lookup[mutation])
        except KeyError:
            possible_voc_sets.append(set())
    lineage_date_combs = set.intersection(*possible_voc_sets)

    mutation_prevalence_dict = defaultdict(int)
    for lineage_date_comb in lineage_date_combs:
        lineage, _ = lineage_date_comb.split("/")
        try:
            for acc_id in genome2snp_dict[lineage_date_comb]:
                if test_set.issubset(genome2snp_dict[lineage_date_comb][acc_id]):
                    if output_acc_id:
                        print(acc_id)
                    mutation_prevalence_dict[lineage] += 1
        except KeyError:
            continue
    return mutation_prevalence_dict

if __name__ == "__main__":
    '''
    main function
    '''
    parser = argparse.ArgumentParser(description="Search for cryptic lineages in wastewater samples.")
    parser.add_argument("-d", "--database", type=str, required=False, help="Path to the Crykey Database Directory. Default:[crykey_dbs]", default="crykey_dbs")
    parser.add_argument("-o", "--output", type=str, required=False, help="Output directory. Default:[crykey_output]", default="crykey_output")
    
    args = parser.parse_args()
    output_dir = args.output
    quarc_db_dir = args.database
    
    merged_df = merge_dataframe(output_dir)
    print("Dataframes merged.\tMerged Dataframe Shape:", merged_df.shape)
    
    file_to_read = open(os.path.join(quarc_db_dir, "quarc_db.pkl"), "rb")
    genome2snp_dict = pickle.load(file_to_read)
    file_to_read = open(os.path.join(quarc_db_dir, "mutation0_lookup.pkl"), "rb")
    sublineage_mutations_exist_lookup = pickle.load(file_to_read)
    print("Database loaded.")
    
    query_result = dict()
    start_time = time.time()
    for i, key in enumerate(merged_df['Nt Mutations'].unique()):
        query_result[key] = mutation_comb_msa_query(set(key.split(';')), sublineage_mutations_exist_lookup, genome2snp_dict)
        if i % 100 == 99:
            print(i+1, 'Average Query Time:', (time.time()-start_time)/(i+1))

    f = open(os.path.join(output_dir, "query_result.pkl"), "wb")
    pickle.dump(query_result, f)
    f.close()
    
    records = []
    for item in query_result:
        occ_list = []
        total_occ = 0
        for l in query_result[item]:
            occ_list.append(f'{l}:{query_result[item][l]}')
            total_occ += query_result[item][l]
        #print(item, total_occ, ";".join(occ_list))
        records.append({
            'Nt Mutations': item,
            'Overall Occurrence': total_occ,
            'Lineage Occurrence': ";".join(occ_list)
        })
    query_df = pd.DataFrame(records)
    merged_df = merged_df.merge(query_df, left_on=['Nt Mutations'], right_on=['Nt Mutations'])
    merged_df.to_csv(os.path.join(output_dir, 'final_result.csv'), index=False)

    print("Done.")
    exit(0)
