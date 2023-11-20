import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

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
from utils.quaid_func import get_counts, get_all_voc, build_count_df, get_recent_nt_df

import argparse

def load_sublineage_mutations_db(quarc_db_path):
    file_to_read = open(os.path.join(quarc_db_path, "mutation50_lookup.pkl"), "rb")
    sublineage_mutations_lookup = pickle.load(file_to_read)
    return sublineage_mutations_lookup

def get_target_variants(vcf_lofreq_f, target_region, min_AF=0.02, min_DP=10):
    '''
    return a list of SNPs in the format (CHROM, POS, REF, ALT)
    only SNPs with AF and DP above threshold and found in the target region are returned
    '''
    vcf_reader = vcf.Reader(open(vcf_lofreq_f, 'r'))
    target_variants = []
    for record in vcf_reader:
        if len(record.REF) == 1 and len(record.ALT[0]) == 1 and\
        record.INFO["AF"] >= min_AF and\
        record.INFO["DP"] >= min_DP and\
        target_region[0] <= record.POS <= target_region[1]:
            target_variants.append(record)
            
    return target_variants

def check_bam_index(sorted_bam_f, force_index=False):
    '''
    check whether BAM file has been indexed, if not, create index.
    '''
    if os.path.exists(sorted_bam_f):
        if (not os.path.exists(f"{sorted_bam_f}.bai")) or force_index:
            subprocess.run(["samtools", "index", sorted_bam_f], check=True)
    else:
        print(f"Error: {sorted_bam_f} does not exist.")
		
def annotate_variants(target_variants, date, site, reference, output_dir):
    '''
    return an annotation lookup table that maps from NT space to AA space
    NT mutation not in the coding region is excluded
    '''
    vcf_reader = vcf.Reader(filename='utils/dummy.vcf') 
    for record in vcf_reader:
        template_record = record
        break
    
    vcf_writer = vcf.Writer(open(os.path.join(output_dir, f'cryptic_vcfs/{date}_{site}_mutations.vcf'), 'w'), vcf_reader)
    for record in target_variants:
        vcf_writer.write_record(record)
    vcf_writer.close()
    
    os.system(f"snpEff ann NC_045512.2 -noStats {output_dir}/cryptic_vcfs/{date}_{site}_mutations.vcf > {output_dir}/cryptic_vcfs/{date}_{site}_mutations.ann.vcf")
    
    vcf_reader = vcf.Reader(filename=f'{output_dir}/cryptic_vcfs/{date}_{site}_mutations.ann.vcf')
    
    annotation_lookup = dict()
    for record in vcf_reader:
        _,_,_,gene_name,_,_,_,_,_,_,hgvs_p,_,_,_,_,_ = record.INFO['ANN'][0].split("|")
        try:
            aa_mut = hgvs_p.split(".")[1]
            ref_aa = seq1(aa_mut[:3])
            if aa_mut[-1] == "*":
                alt_aa = "*"
                pos_aa = aa_mut[3:-1]
            else:
                alt_aa = seq1(aa_mut[-3:], custom_map={"del": "-"})
                pos_aa = aa_mut[3:-3]

            if gene_name == "ORF1ab":
                try:
                    if int(pos_aa) > 4401:
                        pos_aa = str(int(pos_aa)-4401)
                        gene_name = "ORF1b"
                    else:
                        gene_name = "ORF1a"
                except ValueError:
                    print("Error on rename ORF1ab:", date, site, record)
            
            label = ref_aa+pos_aa+alt_aa
            variant_nt_label = str(record.REF)+str(record.POS)+str(record.ALT[0])
            annotation_lookup[variant_nt_label] = f'{gene_name}:{label}'

        except IndexError:
            continue
            
    return annotation_lookup

def get_allele_reads(sorted_bam_f, filtered_variants, read_end_cutoff=5):
    check_bam_index(sorted_bam_f, force_index=True)
    samfile = pysam.AlignmentFile(sorted_bam_f, "rb", check_sq=False)

    read_id_dict = defaultdict(set)
    
    for record in filtered_variants:
        variant_nt_label = str(record.REF)+str(record.POS)+str(record.ALT[0])
        for pileupcolumn in samfile.pileup(record.CHROM, int(record.POS)-1, int(record.POS), min_base_quality=0, max_depth=0):
            # zero-based
            if pileupcolumn.reference_pos == int(record.POS)-1:
                aligned_reads = pileupcolumn.pileups
                for i, base in enumerate(pileupcolumn.get_query_sequences()):
                    #print(base, record.ALT[0])
                    if base == str(record.ALT[0]) or base == str(record.ALT[0]).lower():
                        read_length = aligned_reads[i].alignment.query_length
                        query_pos = pileupcolumn.get_query_positions()[i] # position of the query base
                        dist_to_end = min(abs(read_length-query_pos-1), query_pos) # distance to the ends
                        if dist_to_end > read_end_cutoff:
                        #read_id_dict[aligned_reads[i].alignment.query_name].add(annotation_lookup[pileupcolumn.reference_pos+1])
                            read_id_dict[aligned_reads[i].alignment.query_name].add(variant_nt_label)
    
    return read_id_dict

def fetch_cryptic_alignment(date, site, valid_cryptic_read_ids, sorted_bam_f, output_dir):
    if len(valid_cryptic_read_ids) > 0:
        with open(f"{output_dir}/cryptic_reads/read_ids_{date}_{site}.txt", "w") as read_id_f:
            for read_id in valid_cryptic_read_ids:
                read_id_f.write(f"{read_id}\n")

        subprocess.run(['samtools', 'view', '-b',
                    '-N', f"{output_dir}/cryptic_reads/read_ids_{date}_{site}.txt",
                    '-o', f"{output_dir}/cryptic_alignment/cryptic_reads_{date}_{site}.bam",
                    sorted_bam_f],
                   check=True)

        check_bam_index(f"{output_dir}/cryptic_alignment/cryptic_reads_{date}_{site}.bam", force_index=True)
		
def remove_synonyms_mutation(target_variants, annotation_lookup, execute_remove=False):
    '''
    synonyms mutations and mutations not in the coding regions are removed
    '''
    filtered_variants = []
    for record in target_variants:
        variant_nt_label = str(record.REF)+str(record.POS)+str(record.ALT[0])
        try:
            gene, aa_mutation = annotation_lookup[variant_nt_label].split(":")
            if execute_remove and aa_mutation[0] != aa_mutation[-1]:
                filtered_variants.append(record)
            else:
                filtered_variants.append(record)
        except KeyError:
            continue
            
    return filtered_variants

def search_cryptic_read(read_id_dict, sorted_bam_f, reference, annotation_lookup, sublineage_mutations_lookup, min_supp_dp, min_comb_freq, read_end_cutoff=5):
    '''
    return a list of ids for valid cryptic read,
    return a list of turple containing information for each valid mutation combination, 
    in the format of ((mutations in nt space), (mutations in aa space), number of reads supporting the mutation combination, 
                      total number of reads spanning the mutation sites, combination_frequency)
    A valid mutation combination is defined as the a set of mutation that has number of read support above the user define threshold {min_supp_dp}, 
    and also has the combination frequency above the user define threshold {min_comb_freq}, the combination frequency is calculated as supp_dp/total_dp, 
    where total_dp only counts for the number of read pairs spanning all of the mutation sites in the combination.  
    '''
    cryptic_read_ids = defaultdict(list)
    
    for read_id in read_id_dict:
        if len(read_id_dict[read_id]) > 1:
            possible_voc_sets = []
            for mutation in read_id_dict[read_id]:
                try:
                    possible_voc_sets.append(sublineage_mutations_lookup[mutation])
                except KeyError:
                    possible_voc_sets.append(set())

            if len(set.intersection(*possible_voc_sets)) == 0:
                key_str = list(read_id_dict[read_id])
                key_str.sort(key=lambda test_str : test_str[1:])
                cryptic_read_ids[tuple(key_str)].append(read_id)
    
    valid_cryptic_read_ids = []
    valid_mutation_combinations = []
    for mutation_combination in cryptic_read_ids:
        supp_dp = len(cryptic_read_ids[mutation_combination])
        if supp_dp >= min_supp_dp:
            samfile = pysam.AlignmentFile(sorted_bam_f, "rb", check_sq=False)
            total_read_id_dict = defaultdict(set)
            
            for mutation in mutation_combination:
                pos = int(mutation[1:-1])
                for pileupcolumn in samfile.pileup(reference.id, pos-1, pos, min_base_quality=0, max_depth=0):
                    if pileupcolumn.reference_pos == pos-1:
                        aligned_reads = pileupcolumn.pileups
                        for i, base in enumerate(pileupcolumn.get_query_sequences()):
                            read_length = aligned_reads[i].alignment.query_length
                            query_pos = pileupcolumn.get_query_positions()[i] # position of the query base
                            dist_to_end = min(abs(read_length-query_pos-1), query_pos) # distance to the ends
                            if dist_to_end > read_end_cutoff:
                                total_read_id_dict[aligned_reads[i].alignment.query_name].add(mutation)
            
            total_dp = 0
            for read_id in total_read_id_dict:
                if len(total_read_id_dict[read_id]) == len(mutation_combination):
                    total_dp += 1
            
            comb_freq = supp_dp/total_dp
            
            if comb_freq >= min_comb_freq:
                valid_cryptic_read_ids.extend(cryptic_read_ids[mutation_combination])
                temp_aa = []
                temp_nt = []
                for mutation in mutation_combination:
                    temp_aa.append(annotation_lookup[mutation])
                    temp_nt.append(mutation)
                temp_aa = tuple(temp_aa)
                temp_nt = tuple(temp_nt)
                valid_mutation_combinations.append((temp_nt, temp_aa, supp_dp, total_dp, comb_freq))
                
    return valid_cryptic_read_ids, valid_mutation_combinations

def search_valid_mutation_combinations(date, site, sorted_bam_f, vcf_lofreq_f, cryptic_df, sublineage_mutations_lookup, reference, output_dir, target_region=(1, 29903)):
    if os.path.exists(sorted_bam_f) and os.path.exists(vcf_lofreq_f):
        target_variants = get_target_variants(vcf_lofreq_f, target_region)
        annotation_lookup = annotate_variants(target_variants, date, site, reference, output_dir)
        filtered_variants = remove_synonyms_mutation(target_variants, annotation_lookup)

        read_id_dict = get_allele_reads(sorted_bam_f, filtered_variants)
        valid_cryptic_read_ids, valid_mutation_combinations = search_cryptic_read(read_id_dict, sorted_bam_f, reference, annotation_lookup, sublineage_mutations_lookup, min_supp_dp=5, min_comb_freq=0.01)
        fetch_cryptic_alignment(date, site, valid_cryptic_read_ids, sorted_bam_f, output_dir)

        if valid_mutation_combinations != []:
            for combination in valid_mutation_combinations:
                combination_str_nt = ';'.join(combination[0])
                combination_str_aa = ';'.join(combination[1])
                supp_dp = combination[2]
                total_dp = combination[3]
                comb_freq = combination[4]
                record_dict = {'Date': date,
                               'Site': site,
                               'Nt Mutations': combination_str_nt,
                               'AA Mutations': combination_str_aa,
                               'Support DP': supp_dp,
                               'Total DP': total_dp,
                               'Combined Freq': comb_freq}
                
                #cryptic_df = cryptic_df.append(record_dict, ignore_index=True)
                cryptic_df = pd.concat([cryptic_df, pd.DataFrame([record_dict])], ignore_index=True)
    else:
        if (not os.path.exists(sorted_bam_f)) and os.path.exists(vcf_lofreq_f):
            print(date, site, 'Missing BAM File(s).')
        elif (not os.path.exists(vcf_lofreq_f)) and os.path.exists(sorted_bam_f):
            print(date, site, 'Missing VCF File(s).')
        else:
            print(date, site, 'Missing File(s).')
        
    return cryptic_df

def create_output_directory(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(os.path.join(output_dir, "cryptic_vcfs")):
        os.mkdir(os.path.join(output_dir, "cryptic_vcfs"))
    if not os.path.exists(os.path.join(output_dir, "cryptic_reads")):
        os.mkdir(os.path.join(output_dir, "cryptic_reads"))
    if not os.path.exists(os.path.join(output_dir, "cryptic_alignment")):
        os.mkdir(os.path.join(output_dir, "cryptic_alignment"))
    if not os.path.exists(os.path.join(output_dir, "cryptic_dataframe")):
        os.mkdir(os.path.join(output_dir, "cryptic_dataframe"))

def main():
    '''
    main function
    '''
    # Tested command: python crykey_wastewater.py -i test_input_metadata.tsv -r /home/Users/yl181/wastewater/SARS-CoV-2-reference.fasta -d /home/Users/yl181/wastewater/quaid/quarc_dbs_01102023_incl_recombinant
    
    parser = argparse.ArgumentParser(description="Search for cryptic lineages in wastewater samples.")
    parser.add_argument("-i", "--metadata", type=str, required=True, 
                        help="Path to the input metadata, tab seperated dataframe with the following field on each row:\
                        Sample_Collection_Date[MMDDYYYY], WWTP[str], Sorted_BAM[path-to-sorted-bam-file], VCF[path-to-VCF].")
    parser.add_argument("-r", "--reference", type=str, required=True, help="Path to the Reference Genome.")
    parser.add_argument("-d", "--database", type=str, required=False, help="Path to the Crykey Database Directory. Default:[crykey_dbs]", default="crykey_dbs")
    parser.add_argument("-o", "--output", type=str, required=False, help="Output directory. Default:[crykey_output]", default="crykey_output")
    
    args = parser.parse_args()
    
    quarc_db_path = args.database
    output_dir = args.output
    reference = SeqIO.read(args.reference, "fasta")
    metadata = pd.read_csv(args.metadata, sep='\t')
    
    create_output_directory(output_dir)
    sublineage_mutations_lookup = load_sublineage_mutations_db(quarc_db_path)
    
    for idx, row in metadata.iterrows():
        date = row['Sample_Collection_Date']
        site = row['WWTP']
        sorted_bam_f = row['Sorted_BAM']
        vcf_lofreq_f = row['VCF']

        cryptic_df = pd.DataFrame(columns=['Date', 'Site', 'Nt Mutations', 'AA Mutations', 'Support DP', 'Total DP', 'Combined Freq'])
        cryptic_df = search_valid_mutation_combinations(date, site, sorted_bam_f, vcf_lofreq_f, cryptic_df, sublineage_mutations_lookup, reference, output_dir)
        if not cryptic_df.empty:
            cryptic_df.to_csv(os.path.join(output_dir, "cryptic_dataframe", f'{date}_{site}.csv'), index=False)

        print(f'Sample {date}_{site} processed.')
				
    print('All Samples Processed.')
    
if __name__ == '__main__':
    main()
