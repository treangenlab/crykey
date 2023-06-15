import pandas as pd
import os
import datetime as dt
from datetime import date as date_func
import random

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

from multiprocessing import Pool, Manager
from itertools import repeat

def check_bam_index(sorted_bam_f, force_index=False):
    '''
    check whether BAM file has been indexed, if not, create index.
    '''
    if os.path.exists(sorted_bam_f):
        if (not os.path.exists(f"{sorted_bam_f}.bai")) or force_index:
            subprocess.run(["samtools", "index", sorted_bam_f], check=True)
    else:
        print(f"Error: {sorted_bam_f} does not exist.")
        
def fetch_cryptic_alignment(valid_cryptic_read_ids, sorted_bam_f, sample_id, output_dir):
    if len(valid_cryptic_read_ids) > 0:
        with open(f"{output_dir}/{sample_id}_cryptic_read_ids.txt", "w") as read_id_f:
            for read_id in valid_cryptic_read_ids:
                read_id_f.write(f"{read_id}\n")

        subprocess.run(['samtools', 'view', '-b',
                    '-N', f"{output_dir}/{sample_id}_cryptic_read_ids.txt",
                    '-o', f"{output_dir}/{sample_id}_cryptic_reads.bam",
                    sorted_bam_f],
                   check=True)

        check_bam_index(f"{output_dir}/{sample_id}_cryptic_reads.bam", force_index=True)
        
def get_reads_in_region(sorted_bam_f, target_mutation, reference, read_end_cutoff=5):
    check_bam_index(sorted_bam_f, force_index=False)
    samfile = pysam.AlignmentFile(sorted_bam_f, "rb", check_sq=False)

    read_id_dict = defaultdict(set)
    read_id_dict_fwd = defaultdict(set)
    read_id_dict_rev = defaultdict(set)

    target_mutation_set = set(target_mutation.split(";"))

    for variant_nt_label in target_mutation_set:
        record_pos = int(variant_nt_label[1:-1])
        record_ref = variant_nt_label[0]
        record_alt = variant_nt_label[-1]

        for pileupcolumn in samfile.pileup(reference.id, int(record_pos)-1, int(record_pos), min_base_quality=0, max_depth=0):
            if pileupcolumn.reference_pos == int(record_pos)-1:
                aligned_reads = pileupcolumn.pileups
                for i, base in enumerate(pileupcolumn.get_query_sequences()):
                    read_length = aligned_reads[i].alignment.query_length
                    query_pos = pileupcolumn.get_query_positions()[i] # position of the query base
                    dist_to_end = min(abs(read_length-query_pos-1), query_pos) # distance to the ends
                    if dist_to_end > read_end_cutoff and (not aligned_reads[i].alignment.is_supplementary) and (not aligned_reads[i].alignment.is_secondary):
                        if base == str(record_alt) or base == str(record_alt).lower():
                            read_id_dict[aligned_reads[i].alignment.query_name].add(variant_nt_label)
                            if not aligned_reads[i].alignment.is_reverse:
                                read_id_dict_fwd[aligned_reads[i].alignment.query_name].add(variant_nt_label)
                            else:
                                read_id_dict_rev[aligned_reads[i].alignment.query_name].add(variant_nt_label)
                        else:
                            read_id_dict[aligned_reads[i].alignment.query_name].add(record_ref+str(record_pos)+base.upper())
                            if not aligned_reads[i].alignment.is_reverse:
                                read_id_dict_fwd[aligned_reads[i].alignment.query_name].add(record_ref+str(record_pos)+base.upper())
                            else:
                                read_id_dict_rev[aligned_reads[i].alignment.query_name].add(record_ref+str(record_pos)+base.upper())
                                
    total_dp = 0
    supp_dp = 0
    weak_supp_dp = 0
    valid_cryptic_read_ids = []
    for read_id in read_id_dict:
        if len(read_id_dict[read_id]) >= len(target_mutation_set):
            total_dp += 1
            if target_mutation_set.issubset(read_id_dict[read_id]):
                supp_dp += 1
                valid_cryptic_read_ids.append(read_id)
                if target_mutation_set != read_id_dict[read_id]:
                    weak_supp_dp += 1

    dp1 = 0
    dp3 = 0
    valid_cryptic_read_ids_fwd = []
    valid_cryptic_read_ids_rev = []

    for read_id in read_id_dict_fwd:
        if len(read_id_dict_fwd[read_id]) == len(target_mutation_set):
            if read_id_dict_fwd[read_id] == target_mutation_set:
                dp3 += 1
                valid_cryptic_read_ids_fwd.append(read_id)
            else:
                dp1 += 1

    dp2 = 0
    dp4 = 0
    for read_id in read_id_dict_rev:
        if len(read_id_dict_rev[read_id]) == len(target_mutation_set):
            if read_id_dict_rev[read_id] == target_mutation_set:
                dp4 += 1
                valid_cryptic_read_ids_rev.append(read_id)
            else:
                dp2 += 1
                
    overlapping_supp_dp = len(set(valid_cryptic_read_ids_fwd).intersection(valid_cryptic_read_ids_rev))
    
    return valid_cryptic_read_ids, supp_dp, total_dp, dp1, dp2, dp3, dp4, overlapping_supp_dp, weak_supp_dp

def worker(idx, sample_id, selected_cryptic_mutations, harvest_output_dir, reference):
    sample_cryptic_df = pd.DataFrame(columns=['SRA Accession ID',
                                              'Run Index',
                                              'Nt Mutations',
                                              'Support DP',
                                              'Total DP',
                                              'DP1', 'DP2', 'DP3', 'DP4',
                                              'Overlapping Support DP',
                                              'Inconsistent Support DP',
                                              'Combined Freq'])
    sorted_bam_f = os.path.join(harvest_output_dir, f"{idx}_out", "bam_files", f"{sample_id}.sorted.bam")
    if os.path.exists(sorted_bam_f):
        cryptic_bam_files_dir = os.path.join(harvest_output_dir, f"{idx}_out", "cryptic_bam_files")
        if not os.path.exists(cryptic_bam_files_dir):
            os.mkdir(cryptic_bam_files_dir)
                
        for target_mutation in selected_cryptic_mutations:
            target_mutation_label = target_mutation.replace(";", "_")
            valid_cryptic_read_ids, supp_dp, total_dp, dp1, dp2, dp3, dp4, overlapping_supp_dp, weak_supp_dp = get_reads_in_region(sorted_bam_f, target_mutation, reference)

            output_dir = os.path.join(cryptic_bam_files_dir, target_mutation_label)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            
            if total_dp != 0:
                record_dict = {'SRA Accession ID': sample_id,
                               'Run Index': idx,
                               'Nt Mutations': target_mutation,
                               'Support DP': supp_dp,
                               'Total DP': total_dp,
                               'DP1': dp1,
                               'DP2': dp2,
                               'DP3': dp3,
                               'DP4': dp4,
                               'Overlapping Support DP': overlapping_supp_dp,
                               'Inconsistent Support DP': weak_supp_dp, 
                               'Combined Freq': supp_dp/total_dp}
            else:
                record_dict = {'SRA Accession ID': sample_id,
                               'Run Index': idx,
                               'Nt Mutations': target_mutation,
                               'Support DP': dp3+dp4,
                               'Total DP': dp1+dp2+dp3+dp4,
                               'DP1': dp1,
                               'DP2': dp2,
                               'DP3': dp3,
                               'DP4': dp4,
                               'Overlapping Support DP': overlapping_supp_dp,
                               'Inconsistent Support DP': weak_supp_dp, 
                               'Combined Freq': np.nan}

            sample_cryptic_df = sample_cryptic_df.append(record_dict, ignore_index=True)
            fetch_cryptic_alignment(valid_cryptic_read_ids, sorted_bam_f, sample_id, output_dir)
    else:
        print(idx, sample_id, 'BAM does not exists.')
        
    return sample_cryptic_df

def main():
    reference = SeqIO.read("/home/Users/yl181/wastewater/SARS-CoV-2-reference.fasta", "fasta")
    harvest_output_dir = "/home/Users/yl181/cdc_harvest_variants/Proposed_Pipeline_SRA_Outputs_Quarc"
    
    samples_in_range = []
    with open('us_sampled_sra.txt', 'r') as input_f:
        lines = input_f.readlines()
        for line in lines:
            samples_in_range.append(line.strip())

    selected_cryptic_mutations = ['A29039T;G29049A;A29301G',
                                'G29050A;A29301G',
                                'A26530G;C26577G;G26634A',
                                'T13078C;T13195C',
                                'A26530G;C26533T;C26542A;T26545G',
                                'T28823G;G28881A;G28882A;G28883C',
                                'A27259C;C27335T;A27344T;A27345T',
                                'A26530G;T26545G',
                                'C27807T;A27821C',
                                'C6402T;G6456A',
                                'G22959T;G22992A;C22995A;A23013C;A23040G;G23048A;A23055G;A23063T;T23075C',
                                'G28881A;G28882A;G28883C;G28954T',
                                'A26530G;C26577G;C26625A',
                                'C10449A;T10459C',
                                'T15682A;T15685A',
                                'T29029C;A29039T',
                                'A24966T;C25000T',
                                'A27344T;A27345T;A27354G',
                                'T23599G;C23604A;G23642T',
                                'A29039T;G29049A']
    
    clinical_cryptic_df = pd.DataFrame(columns=['SRA Accession ID',
                                                'Run Index',
                                                'Nt Mutations',
                                                'Support DP',
                                                'Total DP',
                                                'DP1', 'DP2', 'DP3', 'DP4',
                                                'Overlapping Support DP',
                                                'Inconsistent Support DP',
                                                'Combined Freq'])
    
    a_args = []
    b_args = []
    for idx, sample_id in enumerate(samples_in_range):
        a_args.append(idx)
        b_args.append(sample_id)

    pool = Pool(processes=40)
    
    res = pool.starmap(worker, zip(a_args, b_args, repeat(selected_cryptic_mutations), repeat(harvest_output_dir), repeat(reference)))
    clinical_cryptic_df = pd.concat(res, ignore_index=True)
    
    #print(clinical_cryptic_df)
    clinical_cryptic_df.to_csv('us_clinical_20_dp4.csv', index=False)
    
if __name__ == "__main__":
    main()