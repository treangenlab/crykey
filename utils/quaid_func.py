# Standard library imports
import re
from collections import Counter, defaultdict
# 3rd party imports
import pandas as pd

def get_counts(md):
    _genome_counts = md.groupby('Pango lineage')['Accession ID'].resample('W-MON').count().copy()
    return _genome_counts


def get_dates(counts):
    return sorted(counts.loc[:, '2019-12-01':].index.get_level_values(1).unique())


def get_all_voc(level):
    WHO_variants = defaultdict(list)

    if level == 4:
        WHO_variants['Alpha'].append("B.1.1.7")
        WHO_variants['Beta'].append("B.1.351")
        WHO_variants['Epsilon'].append("B.1.427")
        WHO_variants['Epsilon'].append("B.1.429")
        WHO_variants['Eta'].append("B.1.525")
        WHO_variants['Iota'].append("B.1.526")
        WHO_variants['Kappa'].append("B.1.617.1")
        WHO_variants['B.1.617.3'].append("B.1.617.3")
        WHO_variants['Mu'].append("B.1.621")
        WHO_variants['Delta'].append("B.1.617.2")
        WHO_variants['Omicron'].append("B.1.1.529")

    elif level == 5:
        WHO_variants['Gamma'].append("B.1.1.28.1")
        WHO_variants['BA.1'].append("B.1.1.529.1")
        WHO_variants['BA.2'].append("B.1.1.529.2")
        WHO_variants['BA.3'].append("B.1.1.529.3")
        WHO_variants['BA.4'].append("B.1.1.529.4")
        WHO_variants['BA.5'].append("B.1.1.529.5")
        WHO_variants['Zeta'].append("B.1.1.28.2")
        WHO_variants['Lambda'].append("B.1.1.1.37")
    
    return WHO_variants 


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


def build_count_df(counts, _vdb_df):
    _count_df = pd.merge(_vdb_df, counts, left_index=True, right_index=True)
    _count_df.columns = ['SNPs', 'Total count']
    _count_df['SNPs'] = _count_df['SNPs'].map(remove_ambiguous_snps)
    _count_df = _count_df[_count_df['SNPs'] != 0].copy()
    _count_df = _count_df.reset_index()
    _count_df['Collection date'] = pd.to_datetime(_count_df['Collection date'])
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


def lineage_to_variant_map(lineage, WHO_variants):
    lineage = lineage.replace('.0', '')

    for variant in WHO_variants:
        if lineage in WHO_variants[variant]:
            return variant

    return lineage

def get_recent_nt_df(count_df, date, time_window):
    _date = date
    if time_window > 0:
        # Time window mode run
        try:
            recent_data = count_df.loc[pd.date_range(end=_date, periods=time_window, freq='W-MON'), :].copy()
        except KeyError:
            recent_data = count_df.loc[pd.date_range(end=count_df.index.get_level_values(0).max(), periods=time_window, freq='W-MON'), :].copy()
    else:
        # YTD run
        try:
            recent_data = count_df.loc[pd.date_range(start='2019-12-30', end=_date, freq='W-MON'), :].copy()
        except KeyError:
            recent_data = count_df.loc['2019-12-30':, :].copy()

    recent_data = recent_data.groupby('Pango lineage').sum(numeric_only=False)

    recent_data_freq = []
    for i, row in recent_data.iterrows():
        for mut in row['SNPs']:
            recent_data_freq.append((i, mut, 
                                     row['SNPs'][mut],
                                     row['Total count']))
    recent_nt_df = pd.DataFrame(data=recent_data_freq, columns=['Lineage', 'NT_mut', 'Prevalence', 'Occurence count'])
    return recent_nt_df

def get_quasiunique_df(recent_nt_df, WHO_variants, alias_data, occurence_cutoff, incl_cutoff, excl_cutoff, level=4, verbose=True):
    recent_nt_df['Lineage'] = recent_nt_df['Lineage'].replace(alias_data, regex=True)
    recent_nt_df = recent_nt_df.join(recent_nt_df['Lineage'].str.split('.', expand=True).fillna(str(0)))

    lineage_lvls = recent_nt_df.columns[4:]
    recent_nt_df = recent_nt_df.set_index(['NT_mut'] + list(lineage_lvls))

    recent_nt_df = recent_nt_df[recent_nt_df['Prevalence'] > occurence_cutoff].copy()
    recent_nt_df = recent_nt_df[(recent_nt_df['Prevalence'] / recent_nt_df['Occurence count']) >= excl_cutoff]

    to_exclude = recent_nt_df.groupby(level=[i for i in range(level+1)]).sum().groupby('NT_mut').count()['Prevalence'] > 1
    to_exclude = to_exclude.index[to_exclude]
    recent_nt_df = recent_nt_df.drop(to_exclude, axis=0, level='NT_mut')

    recent_nt_df = recent_nt_df[(recent_nt_df['Prevalence'] / recent_nt_df['Occurence count']) >= incl_cutoff]
    recent_nt_df = recent_nt_df.groupby(level=[i for i in range(level+1)]).sum()

    recent_nt_df = recent_nt_df.reset_index()

    recent_nt_df['WHO_name'] = recent_nt_df.iloc[:, 1:level+1].agg('.'.join, axis=1).apply(lambda x: lineage_to_variant_map(x, WHO_variants))

    return recent_nt_df


def get_result_table(vcf_df, quasiunique_df, coverage_data, WHO_variants, date, min_cov=0., verbose=True):
    vcf_data = vcf_df.merge(quasiunique_df, how='inner', on='NT_mut')
    
    table_data = {"Date": [], 
                  "Plant": [],
                  "Total AF (quasi-unique)": [],
                  "Total QU count": [],
                  "Number of quasi-unique mutations possible": [],
                  "Fraction of quasi-unique sites covered": [], 
                  "WHO name": []}

    plants = sorted(vcf_data.Plant.unique(), reverse=True)
    for variant in WHO_variants:
        lin_mut_nt = quasiunique_df[quasiunique_df.WHO_name == variant]
        for p in plants:
            lg = len(vcf_data[(vcf_data.Date == date) & (vcf_data.Plant == p) & (vcf_data.WHO_name == variant)].ALT_FREQ_LoFreq)
            if lg == 0:
                val = 0.
            else:
                val = sum(vcf_data[(vcf_data.Date == date) & \
                                   (vcf_data.WHO_name == variant) & \
                                   (vcf_data.Plant == p)].ALT_FREQ_LoFreq)
                val = float(val)

            try:
                cov_data = coverage_data[f"{pd.to_datetime(date).strftime('%m%d%Y')}_{p}"]
                
                for j, nt_mut in lin_mut_nt.iterrows():
                    cov = 0

                    try:
                        if nt_mut[-1] == "-":
                            cov = int(cov_data[cov_data[1] == (int(nt_mut['NT_mut'][1:-1]) - 1)][2]) + \
                                  int(cov_data[cov_data[1] == (int(nt_mut['NT_mut'][1:-1]) + 1)][2])
                        else:
                            cov = int(cov_data[cov_data[1] == int(nt_mut['NT_mut'][1:-1])][2])
                            
                    except TypeError:
                        cov = 0.

                    if cov > min_cov:
                        qu_frac_cov += 1
                
                if len(lin_mut_nt) > 0:
                    qu_frac_cov /= len(lin_mut_nt)
                else:
                    qu_frac_cov = 0.
            
            except KeyError:
                qu_frac_cov = 0.


            table_data['Date'].append(date)
            table_data['Plant'].append(p)
            table_data['Total AF (quasi-unique)'].append(val)
            table_data['Total QU count'].append(lg)
            table_data['Number of quasi-unique mutations possible'].append(len(quasiunique_df[quasiunique_df['WHO_name'] == variant]))
            table_data['Fraction of quasi-unique sites covered'].append(qu_frac_cov)
            table_data['WHO name'].append(variant)

    result_table = pd.DataFrame(data=table_data)

    return result_table


def get_result_table_sim(vcf_df, quasiunique_df, WHO_variants, date, verbose=True):
    vcf_data = vcf_df.merge(quasiunique_df, how='inner', on='NT_mut')
    
    table_data = {"Date": [], 
                  "Total AF (quasi-unique)": [],
                  "Total QU count": [],
                  "Number of quasi-unique mutations possible": [],
                  "WHO name": []}

    for variant in WHO_variants:
        lin_mut_nt = quasiunique_df[quasiunique_df.WHO_name == variant]
        lg = len(vcf_data[(vcf_data.Date == date) & (vcf_data.WHO_name == variant)].ALT_FREQ_LoFreq)
        if lg == 0:
            val = 0.
        else:
            val = sum(vcf_data[(vcf_data.Date == date) & \
                                (vcf_data.WHO_name == variant)].ALT_FREQ_LoFreq)
            val = float(val)

        table_data['Date'].append(date)
        table_data['Total AF (quasi-unique)'].append(val)
        table_data['Total QU count'].append(lg)
        table_data['Number of quasi-unique mutations possible'].append(len(quasiunique_df[quasiunique_df['WHO_name'] == variant]))
        table_data['WHO name'].append(variant)

    result_table = pd.DataFrame(data=table_data)

    return result_table