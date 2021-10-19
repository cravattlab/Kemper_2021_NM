#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 15:03:08 2021

@author: estherkemper
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:03:56 2019
##This is for individual reps of mitosis +/- lpp (no asynch, no LPP after)
##LPP1 is the experimental group, just comparing CTRL1 and LPP1 
This combines everything together, then filters based on LPP/CTRL
@author: estherkemper
"""
import numpy as np
import pandas as pd
import os

dir = os.getcwd()
#set a target cutoff defined in table
target_cutoff = 2
date = '20210830'

protein_cols = ['accession','description', 'protein', 'accession_res','gene_res']
lpp_med = 'cysteine_Mitosis_LPP(-,-)/(+,-)_combined_median'
lpp_no = 'cysteine_Mitosis_LPP(-,-)/(+,-)_combined_no'
lpp_cat = 'cysteine_Mitosis_LPP(-,-)/(+,-)_category'
basal_med = 'cysteine_Mitosis/Asynch_combined_median'
basal_cat = 'cysteine_Mitosis/Asynch_category'

numerator = 'Mitosis_CTRL1'
basal_denom = 'Asynch_CTRL2'
lpp_denom = 'Mitosis_LPP1'

def combine():
    print('Reading input datasets')
    lpp_df = pd.read_csv('{}/1_scripts/combined/{}_combined_phospho_lpp.csv'.format(dir,date))
    basal_df = pd.read_csv('{}/1_scripts/combined/{}_combined_phospho_asynch.csv'.format(dir,date))
    lpp_prot = pd.read_csv('{}/1_scripts/combined/{}_combined_protein_wp_lpp.csv'.format(dir,date))
    basal_prot = pd.read_csv('{}/1_scripts/combined/{}_combined_protein_whole_proteome.csv'.format(dir,date))
    basal_prot.columns = basal_prot.columns.str.replace(numerator+'/'+basal_denom, 'Mitosis/Asynch')
    lpp_prot.columns = lpp_prot.columns.str.replace(numerator + '/' + lpp_denom,'Mitosis_LPP(-,-)/(+,-)')
    basal_df.columns = basal_df.columns.str.replace(numerator+'/'+basal_denom, 'Mitosis/Asynch')
    lpp_df.columns = lpp_df.columns.str.replace(numerator + '/' + lpp_denom,'Mitosis_LPP(-,-)/(+,-)')
    protein = lpp_prot.merge(basal_prot, how = 'outer', on = ['accession','description','protein'], suffixes = [' protein_lpp',' protein_basal'])
    protein.columns = ['protein_' + x if x not in protein_cols else x for x in protein.columns]
    basal_df.columns = ['cysteine_' + x if x not in protein_cols else x for x in basal_df.columns]
    lpp_df.columns = ['cysteine_' + x if x not in protein_cols else x for x in lpp_df.columns]
    df = lpp_df.merge(basal_df, how = 'outer', on = protein_cols, suffixes = [' cys_lpp', ' cys_basal'])
    df = df.merge(protein, how = 'left', on = ['accession','description','protein'], suffixes = ['', ' protein'])
    rep_cols = [col for col in df.columns if '_channel_' in col]
    df.drop(rep_cols, axis = 1, inplace = True)
    df.columns = df.columns.str.replace('_ratio_','_')
    return df

def sequence_combine(row):
     reps = [col for col in row.index if 'sequence' in col]
     sequences = row[reps].dropna()
     shortest_sequence =  min(sequences, key=len)
     return pd.Series(shortest_sequence)
 
def phospho_merge():
    df = combine()
    df.rename(columns = {'sequence':'sequence cys_lpp'}, inplace = True)
    #get rid of empty columns
    df = df.dropna(axis = 1, how = 'all')
    sequence_cols = [col for col in df.columns if 'sequence' in col]
    df['sequence'] = df[sequence_cols].apply(sequence_combine, axis = 1)
    df.drop(sequence_cols, axis = 1 , inplace = True)
    df['gene'] = df.gene_res.str.split(pat = '_', expand = True)[0]
    phospho_df = pd.read_csv('{}/other_datasets/final_phospho_mann_olsen.csv'.format(dir))
    keep_cols = [col for col in phospho_df.columns if '_olsen' not in col]
    phospho_df = phospho_df[keep_cols]
    phospho_col = [col for col in phospho_df.columns if 'accession' not in col]
    df1 = df.merge(phospho_df, on = 'accession', how = 'left')
    print('Finished merging')
    return df1, phospho_col

def qc_stuff():
    df1, phospho_col = phospho_merge()
    
    comb_median_cols = sorted([col for col in df1.columns if '_combined_median' in col])
    exp_median_cols = sorted([col for col in df1.columns if 'median_TMT' in col])
    CV_cols = sorted([col for col in df1.columns if 'combined_CV' in col] )
    no_cols = sorted([col for col in df1.columns if 'combined_no' in col])
    category_cols = sorted([col for col in df1.columns if 'category' in col])
    keep_protein_cols = ['accession', 'description', 'protein', 'sequence', 'accession_res','gene_res']
    df2 = df1[keep_protein_cols + phospho_col + category_cols + comb_median_cols + CV_cols + no_cols + exp_median_cols].copy()
    keep_protein_cols = ['accession', 'description','protein','sequence','accession_res', 'gene_res']
    category_cols = [col for col in df1.columns if 'category' in col]
    comb_median_cols = sorted([col for col in df1.columns if 'combined_median' in col])
    rep_ratio_cols = sorted([col for col in df1.columns if 'median_TMT' in col])
    cv_cols = sorted([col for col in df1.columns if '_CV' in col])
    no_cols = sorted([col for col in df1.columns if '_no' in col])
    sdf = df1[keep_protein_cols + category_cols + phospho_col + comb_median_cols + cv_cols + no_cols + rep_ratio_cols ].copy()
    sdf = sdf.sort_values(lpp_med, ascending = False)
    sdf.loc[(sdf[lpp_cat].str.contains('CV too high', na = True)), lpp_med] = 'NQ'
    # sdf.loc[(sdf[lpp_no]<2) | (sdf[lpp_no].isnull()), lpp_med] = 'NQ'
    sdf.loc[(sdf[basal_cat].str.contains('CV too high', na = True)), basal_med] = 'NQ'
    # sdf.loc[(sdf[y_no]<1) | (sdf[y_no].isnull()), y_med] = 'NQ'
    df2 = sdf.loc[~(sdf[comb_median_cols] == 'NQ').all(axis = 1)]
    df2 = df2.replace('NQ', np.nan)
#     sdf = sdf.assign(accession_res = sdf.accession_res.str.rsplit('_', n = 1, expand = True)[0])
    sdf.drop(category_cols, axis =1, inplace = True)
    sdf = sdf.replace('NQ', np.nan)
    sdf[comb_median_cols + exp_median_cols] = np.log2(sdf[comb_median_cols + exp_median_cols])
    sdf_changing = sdf.loc[sdf[lpp_med] >= np.log2(target_cutoff)]
    sdf_changing[comb_median_cols] = sdf_changing[comb_median_cols].replace(np.nan, 'NQ')
    sdf[comb_median_cols] = sdf[comb_median_cols].replace(np.nan, 'NQ')
    return sdf, df2, sdf_changing

def main():
    sdf, df2, sdf_changing = qc_stuff()
    print('Writing to excel1')
    with pd.ExcelWriter('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir,date)) as writer:  # doctest: +SKIP
        df2.to_excel(writer, sheet_name = 'filtered', index = False)
        sdf.to_excel(writer, sheet_name = 'supplemental', index = False, float_format = '%.2f')
        sdf_changing.to_excel(writer, sheet_name = 'changing', index = False, float_format = '%.2f')
         
main()


