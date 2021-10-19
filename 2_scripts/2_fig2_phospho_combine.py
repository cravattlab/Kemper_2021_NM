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
original_med = 'Mitosis_LPP(-,-)/(+,-)_combined_median'
original_no = 'Mitosis_LPP(-,-)/(+,-)_combined_no'
adapted_no = 'Mitosis_LPP(-,+)/(+,+)_combined_no'
original_cat = 'Mitosis_LPP(-,-)/(+,-)_category'
adapted_med = 'Mitosis_LPP(-,+)/(+,+)_combined_median'
adapted_cat = 'Mitosis_LPP(-,+)/(+,+)_category'

lpp = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')

numerator = 'Mitosis_CTRL1'
adapted_num = 'Mitosis_CTRL3'
basal_denom = 'Asynch_CTRL2'
original_denom = 'Mitosis_LPP1'
adapted_denom = 'Mitosis_LPP3'

def combine():
    print('Reading input datasets')
    original_df = pd.read_csv('{}/1_scripts/combined/{}_phospho_combined_all_original.csv'.format(dir,date))
    adapted_df = pd.read_csv('{}/1_scripts/combined/{}_phospho_combined_adapted.csv'.format(dir,date))
    basal_df = pd.read_csv('{}/1_scripts/combined/{}_phospho_combined_basal.csv'.format(dir,date))
    basal_df.columns = ['cysteine_' + x if x not in protein_cols else x for x in basal_df.columns]
    basal_df.columns = basal_df.columns.str.replace(numerator+'/'+basal_denom, 'Mitosis/Asynch')
    original_df.columns = original_df.columns.str.replace(numerator + '/' + original_denom,'Mitosis_LPP(-,-)/(+,-)')
    adapted_df.columns = adapted_df.columns.str.replace(adapted_num + '/' + adapted_denom,'Mitosis_LPP(-,+)/(+,+)')
    print('Merging data')
    df = original_df.merge(basal_df, how = 'outer', on = protein_cols, suffixes = [' original', ' basal'])
    df = df.merge(adapted_df, how = 'outer', on = protein_cols, suffixes = ['', ' adapted'])
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
    df.rename(columns = {'sequence':'sequence original'}, inplace = True)
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
    #reorder
    df2 = df1[keep_protein_cols + phospho_col + category_cols + comb_median_cols + CV_cols + no_cols + exp_median_cols].copy()
    sdf = df1[keep_protein_cols + phospho_col + category_cols + comb_median_cols + CV_cols + no_cols + exp_median_cols].copy()
    sdf[comb_median_cols + exp_median_cols] = np.log2(sdf[comb_median_cols + exp_median_cols])
    sdf = sdf.sort_values(original_med, ascending = False)
    sdf.loc[sdf[original_no].isnull(), original_med] = 'NQ'
    sdf.loc[(sdf[original_no] < 2), original_med] = 'NQ'
    sdf.loc[(sdf[original_cat].str.contains('CV too high')) & (~sdf[original_cat].isnull()), original_med] = 'NQ'
    cv_med_cols = [col for col in comb_median_cols if 'protein' not in col]
    sdf_cols = category_cols
    for col in sdf_cols:
        name = col.rsplit(sep = '_', maxsplit = 1)[0]
        med_col = [col for col in cv_med_cols if name in col][0]
        cat_col = [col for col in sdf_cols if name in col][0]
        no_col = [col for col in no_cols if name in col][0]
        sdf.loc[(sdf[cat_col].str.contains('CV too high', na = True)), med_col] = 'NQ'
        sdf.loc[(sdf[cat_col].str.contains('CV too high', na = True)), med_col] = 'NQ'
        sdf.loc[(sdf[no_col]<2) & (sdf[original_med] == 'NQ'), med_col] = 'NQ'
    df3 = df2.loc[df2[original_no] >=2]
    df3 = df3.loc[~df3[original_cat].str.contains('CV too high')]
    for col in sdf_cols:
        name = col.rsplit(sep = '_', maxsplit = 1)[0]
        med_col = [col for col in cv_med_cols if name in col][0]
        cat_col = [col for col in sdf_cols if name in col][0]
        df3.loc[(df3[cat_col].str.contains('CV too high', na = True)), med_col] = np.nan
        df3.loc[~ (df3[cat_col].str.contains('CV too high', na = True)), med_col] = df3[med_col]
    # df3 = qc_prep(df3)
    # sdf = qc_prep(sdf)
    df3.loc[(df3[adapted_cat].str.contains('Delete', na = False)), adapted_med] = np.nan
    sdf.loc[(sdf[adapted_cat].str.contains('Delete', na = False)), adapted_med] = 'NQ'
    # sdf = sdf.loc[ (sdf[original_no].notnull()) & (sdf[adapted_no].notnull())]
    sdf[original_med] = sdf[original_med].replace('NQ', np.nan)
    sdf = sdf.sort_values(original_med, ascending = False)
    sdf[original_med] = sdf[original_med].replace(np.nan, 'NQ')
    print(sdf)
    return sdf, df3

def lpp_merge():
    sdf, df3 = qc_stuff()
    p_med = 'cysteine_' + original_med
    print(lpp.columns)
    lpp1 = lpp[['accession', 'gene_res', p_med]]
    lpp1 = lpp1.round(decimals = 1)
    lpp2 = lpp1.loc[lpp1[p_med] >= 2]
    lpp2 = pd.DataFrame(lpp2.groupby('accession')['gene_res'].apply(lambda row: ', '.join(row.astype(str)))).reset_index()
    lpp2.rename({'gene_res':'LPP-sensitive phosphosite'}, axis =1, inplace = True)
    sdf1 = sdf.merge(lpp2, on = 'accession', how = 'left')
    return sdf1, df3

def main():
    sdf, df3 = lpp_merge()
    drop_cols = [col for col in sdf.columns if 'category' in col]
    sdf.drop(drop_cols, axis =1, inplace = True)
    print('Writing to excel1')
    with pd.ExcelWriter('{}/supplemental_tables/{}_fig2_phospho_combined.xlsx'.format(dir,date)) as writer:  # doctest: +SKIP
        df3.to_excel(writer, sheet_name = 'filtered', index = False)
        sdf.to_excel(writer, sheet_name = 'supplemental', index = False, float_format = '%.2f')
         
main()


