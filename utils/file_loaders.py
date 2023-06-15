"""A set of utility file loaders for QuaID and associated pipeline(s)."""
# Standard library imports
import logging
import json
import glob
import re
from collections import Counter  # need import Counter for eval() to work correctly within the `load_vdb_df()`
# 3rd party imports
import pandas as pd

logger = logging.getLogger('Loader')

def load_metadata(path):
    logger.info(f'Loading {path}')
    _metadata = pd.read_csv(path, sep='\t', low_memory=False)
    _metadata = _metadata[(_metadata['Is complete?'] == True) & (_metadata['Host'] == 'Human')]
    _metadata = _metadata[['Pango lineage', 'Collection date', 'Accession ID']].copy()
    _metadata = _metadata[_metadata['Collection date'].str.len() == 10].copy()
    _metadata['Collection date'] = pd.to_datetime(_metadata['Collection date'])
    _metadata = _metadata[_metadata['Pango lineage'].notna()].copy()
    _metadata = _metadata.set_index(['Collection date'])
    return _metadata


def load_vcf(folder_path):
    logger.info(f'Loading VCF files from {folder_path}')
    _vcf_data = pd.DataFrame()
    data_list = []
    for fn in glob.iglob(f"{folder_path}/*/*.tsv"):
        filename = fn.split("/")[-1]
        df = pd.read_csv(fn, sep='\t')
        plant = filename.split("-")[0]
        plant = re.sub(r'([A-Za-z]+)[0-9]', r'\1', plant).upper()
        date = fn.split("/")[-2]
        date = date[3:]
        plant_date = {"Date": [date for i in range(len(df))],
                      "Plant": [plant for i in range(len(df))]}
        df = pd.concat([df, pd.DataFrame(data=plant_date)], axis=1)
        data_list.append(df)
        
    _vcf_data = pd.concat(data_list, ignore_index=True) 
    _vcf_data = _vcf_data[~(_vcf_data.ALT_AA == _vcf_data.REF_AA)].copy()
    vcf_indel_data = _vcf_data[_vcf_data.REF_LoFreq.str.len() > 1].copy()

    new_rows = []
    for i, row in vcf_indel_data.iterrows():
        for j in range(1, len(row['REF_LoFreq'])):
            new_row = row.copy()
            new_row["REF_LoFreq"] = row['REF_LoFreq'][j]
            new_row["ALT_LoFreq"] = "-"
            new_row["POS"] = row['POS'] + j
            new_rows.append(new_row)

    new_rows = pd.concat(new_rows, ignore_index=True)
    _vcf_data = pd.concat([_vcf_data, new_rows], ignore_index=True)

    _vcf_data = _vcf_data.dropna(subset=["POS"])
    _vcf_data["POS"] = _vcf_data["POS"].apply(int)
    _vcf_data["NT_mut"] = _vcf_data["REF_LoFreq"] + _vcf_data["POS"].apply(str) + _vcf_data["ALT_LoFreq"]
    _vcf_data['Variant'] = _vcf_data.Variant_LoFreq

    _vcf_data['Date'] = pd.to_datetime(_vcf_data['Date'], format='%m%d%Y')
    return _vcf_data


def load_simulated_vcf(folder_path):
    logger.info(f'Loading VCF files from {folder_path}')
    _vcf_data = pd.DataFrame()
    data_list = []
    for fn in glob.iglob(f"{folder_path}/*.tsv"):
        filename = fn.split("/")[-1]
        df = pd.read_csv(fn, sep='\t')
        date = filename.split('_')[2]
        plant_date = {"Date": [date for i in range(len(df))]}
        df = pd.concat([df, pd.DataFrame(data=plant_date)], axis=1)
        data_list.append(df)
        
    _vcf_data = pd.concat(data_list, ignore_index=True) 
    _vcf_data = _vcf_data[~(_vcf_data.ALT_AA == _vcf_data.REF_AA)].copy()
    vcf_indel_data = _vcf_data[_vcf_data.REF_LoFreq.str.len() > 1].copy()

    new_rows = []
    for i, row in vcf_indel_data.iterrows():
        for j in range(1, len(row['REF_LoFreq'])):
            new_row = row.copy()
            new_row["REF_LoFreq"] = row['REF_LoFreq'][j]
            new_row["ALT_LoFreq"] = "-"
            new_row["POS"] = row['POS'] + j
            new_rows.append(new_row)

    new_rows = pd.concat(new_rows, ignore_index=True)
    _vcf_data = pd.concat([_vcf_data, new_rows], ignore_index=True)

    _vcf_data = _vcf_data.dropna(subset=["POS"])
    _vcf_data["POS"] = _vcf_data["POS"].apply(int)
    _vcf_data["NT_mut"] = _vcf_data["REF_LoFreq"] + _vcf_data["POS"].apply(str) + _vcf_data["ALT_LoFreq"]
    _vcf_data['Variant'] = _vcf_data.Variant_LoFreq

    _vcf_data['Date'] = pd.to_datetime(_vcf_data['Date'])
    return _vcf_data


def load_alias_data():
    """Load lineage name alias data."""
    with open("alias_key.json", 'r') as inf:
        alias_data = json.load(inf)

    alias_data = {fr"^{k}.(.*)": fr"{v}.\1" for k, v in alias_data.items() if isinstance(v, str)}
    del alias_data[r'^A.(.*)']
    del alias_data[r'^B.(.*)']
    
    return alias_data


def load_vdb_mutation_data(vdb_f):
    """Load trimmed nucleotide output from the vdb processed GISAID multiple sequence alignment.

    Positional arguments:
    vdb_f -- path to the vdb output file (e.g. 'MSA_0514/vdb_051422_trimmed_nucl.txt')

    Returns:
    A pandas dataframe indexed by GISAID accession number with all SNPs expanded into their own columns.
    Reference genome (assumed to be in the first row by construction) is ommited from the dataframe.
    """
    logger.info(f'Loading {vdb_f}')
    _mutations_data = pd.read_csv(vdb_f, header=None)
    _mutations_data[["Name", "ID", "Date", "Region"]] = _mutations_data[0].str.split("|", expand=True)
    _mutations_data.rename({1: "SNPs"}, axis=1, inplace=True)
    _mutations_data.drop([0, "Name"], axis=1, inplace=True)
    _mutations_data.set_index("ID", inplace=True)
    _mutations_data = _mutations_data[1:]
    _mutations_data["SNPs"] = _mutations_data["SNPs"].str.upper()

    return _mutations_data


def load_coverage_data(coverage_folder_path):
    _coverage_data = {}
    for fn in glob.iglob(f"{coverage_folder_path}/*/*/*.txt"):
        path_split = fn.split('/')
        date = path_split[-3]
        date = date[-len("MMDDYYYY"):]
        plant = path_split[-1].split(".")[0].split("-")[0]
        plant = re.sub(r'([A-Za-z]+)[0-9]', r'\1', plant)

        try:
            df = pd.read_csv(fn, sep='\t', header=None)
            plant_date = {"Plant": [plant for i in range(len(df))],
                          "Date": [date for i in range(len(df))]}
            df = pd.concat([df, pd.DataFrame(data=plant_date)], axis=1)
            _coverage_data[f"{date}_{plant}"] = df
        except pd.errors.EmptyDataError:
            pass

    return _coverage_data


def load_vdb_df(path):
    logger.info(f'Loading {path}')
    _vdb_df = pd.read_csv(path)
    _vdb_df['Collection date'] = pd.to_datetime(_vdb_df['Collection date'])
    _vdb_df = _vdb_df.set_index(['Pango lineage', 'Collection date'])
    _vdb_df = _vdb_df.applymap(lambda x: eval(x))  # need import Counter for this to work
    return _vdb_df


def merge_data(_vdb_df, _metadata):
    _metadata = _metadata.reset_index()
    _vdb_df = _vdb_df.merge(_metadata,
                            left_on=["ID"], right_on=["Accession ID"],
                            validate="one_to_one")
    _vdb_df['Collection date'] = pd.to_datetime(_vdb_df['Collection date'])
    _vdb_df = _vdb_df.dropna()

    return _vdb_df
