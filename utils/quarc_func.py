# Standard library imports
import re
from collections import Counter, defaultdict
# 3rd party imports
import pandas as pd

def remove_ambiguous_snps(x):
    if isinstance(x, Counter):
        if '' in x:
            del x['']

        to_remove = Counter()
        for key in x:
            if re.match('[A,C,T,G][0-9]+[A,C,T,G,-]+', key) is None:
                to_remove[key] = x[key] 

        x = x - to_remove
        
    return +x

def build_count_df(counts, _vdb_df, exclude_recombinant=True):
    _count_df = pd.merge(_vdb_df, counts, left_index=True, right_index=True)
    _count_df.columns = ['SNPs', 'Total count']
    _count_df['SNPs'] = _count_df['SNPs'].map(remove_ambiguous_snps)
    _count_df = _count_df[_count_df['SNPs'] != 0].copy()
    _count_df = _count_df.reset_index()
    _count_df['Collection date'] = pd.to_datetime(_count_df['Collection date'])
    if exclude_recombinant:
        _count_df = _count_df[~(_count_df['Pango lineage'].str.startswith('X'))].copy()  # Exclude recombinant lineages
    _count_df = _count_df[~(_count_df['Pango lineage'] == 'Unassigned')].copy()      # Exclude unassigned lineages
    _count_df = _count_df[~(((_count_df['Pango lineage'].str.startswith('BA.')) | (_count_df['Pango lineage'] == 'B.1.1.529')) &\
                            (_count_df['Collection date'] < '2021-09-01'))].copy()             # Exclude Omicron prior to first sequence
    _count_df = _count_df[~(((_count_df['Pango lineage'].str.startswith('AY.')) | (_count_df['Pango lineage'] == 'B.1.617.2')) &\
                            (_count_df['Collection date'] < '2021-03-01'))].copy()             # Exclude Delta prior to first sequence
    _count_df = _count_df[~(((_count_df['Pango lineage'].str.startswith('Q.')) | (_count_df['Pango lineage'] == 'B.1.1.7')) &\
                            (_count_df['Collection date'] < '2020-09-03'))].copy()             # Exclude Alpha prior to first sequence
    _count_df = _count_df.set_index(['Collection date', 'Pango lineage'])
    _count_df = _count_df.sort_index()
    return _count_df
