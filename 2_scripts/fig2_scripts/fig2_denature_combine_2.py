#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 12:38:17 2021

@author: ekk
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
x_med = 'Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'Mitosis_Native/Denatured_combined_median'
x_no = 'Mitosis_LPP(-,-)/(+,-)_combined_no'
y_no = 'Mitosis_Native/Denatured_combined_no'
x_cat ='Mitosis_LPP(-,-)/(+,-)_category'
y_cat ='Mitosis_Native/Denatured_category'

numerator = 'Mitosis_CTRL1'
original_denom = 'Mitosis_LPP1'
denature_denom = 'Mitosis_urea1'

expt_table = pd.read_csv('{}/inputs/protein_inputs.csv'.format(dir), index_col = None)
lpp = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')


def combine():
    denature_df = pd.read_csv('{}/1_scripts/combined/{}_combined_denature.csv'.format(dir,date))
    original_df = pd.read_csv('{}/1_scripts/combined/{}_combined_all_original.csv'.format(dir,date))
    original_df.columns = original_df.columns.str.replace(numerator + '/' + original_denom,'Mitosis_LPP(-,-)/(+,-)')
    denature_df.columns = denature_df.columns.str.replace(numerator + '/' + denature_denom,'Mitosis_Native/Denatured')
    df = original_df.merge(denature_df, how = 'outer', on = protein_cols, suffixes = [' original', ' denature'])
    rep_cols = [col for col in df.columns if '_channel_' in col]
    df.drop(rep_cols, axis = 1, inplace = True)
    df.columns = df.columns.str.replace('_ratio_','_')
    print('Merging data')
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

def df_reorder():
    #reordering columns
    df1, phospho_col = phospho_merge()
    # sequence_cols = [col for col in df1.columns if 'sequence' in col]
    # df1['sequence'] = df1[sequence_cols].apply(sequence_combine, axis = 1)
    # df1.drop(sequence_cols, axis = 1 , inplace = True)
    keep_protein_cols = ['accession','description','protein','sequence','accession_res', 'gene_res']
    category_cols = [col for col in df1.columns if 'category' in col]
    comb_median_cols = sorted([col for col in df1.columns if 'combined_median' in col])
    rep_ratio_cols = sorted([col for col in df1.columns if 'median_TMT' in col])
    cv_cols = sorted([col for col in df1.columns if '_CV' in col])
    no_cols = sorted([col for col in df1.columns if '_no' in col])
    all_reps = sorted([col for col in df1.columns if 'channel' in col])
    sdf = df1[keep_protein_cols + phospho_col + category_cols + comb_median_cols +cv_cols + no_cols + rep_ratio_cols ].copy()
    df2 = df1[keep_protein_cols + phospho_col + category_cols + comb_median_cols +cv_cols + no_cols + rep_ratio_cols ].copy()
    df2[comb_median_cols + rep_ratio_cols] = np.log2(df2[comb_median_cols + rep_ratio_cols])
    sdf[comb_median_cols + rep_ratio_cols] = np.log2(sdf[comb_median_cols + rep_ratio_cols])
    sdf = sdf.sort_values(y_med, ascending = False)

    sdf.loc[(sdf[x_cat].str.contains('CV too high', na = True)), x_med] = 'NQ'
    sdf.loc[(sdf[x_no]<2) | (sdf[x_no].isnull()), x_med] = 'NQ'
    sdf.loc[(sdf[y_cat].str.contains('CV too high', na = True)), y_med] = 'NQ'
    sdf = sdf.loc[sdf[y_no].notnull()]
    sdf.loc[(sdf[y_no]<1) | (sdf[y_no].isnull()), y_med] = 'NQ'

    df2.loc[(df2[x_cat].str.contains('CV too high', na = True)), x_med] = np.nan
    df2.loc[(df2[x_no]<2) | (df2[x_no].isnull()), x_med] = np.nan
    df2.loc[(df2[y_cat].str.contains('CV too high', na = True)), y_med] = np.nan
    df2.loc[(df2[y_no]<1) | (df2[y_no].isnull()), y_med] = np.nan

    sdf.drop(category_cols, axis = 1, inplace = True)
    sdf = sdf.assign(accession_res = sdf.accession_res.str.rsplit('_', n =1, expand = True)[0])
    diff = df2.loc[(df2[y_med] >= np.log2(target_cutoff)) | (df2[y_med] <= np.log2(1/target_cutoff))]
    diff[comb_median_cols] = diff[comb_median_cols].replace(np.nan, 'NQ')
    
    prot_ids = list(set(diff['accession']))
    prots = df2.loc[df2.accession.isin(prot_ids)].copy()
    prots = prots.assign(accession_res = prots.accession_res.str.rsplit('_', n =1, expand = True)[0])
    prots.drop(category_cols, axis = 1, inplace = True)
    prots.loc[prots.protein == 'Uncharacterized', 'protein'] = 'Uncharacterized_' + prots.accession
    prots.sort_values('protein', ascending = True, inplace = True)
    # prots[comb_median_cols] = np.log2(prots[comb_median_cols])
    prots[comb_median_cols] = prots[comb_median_cols].replace(np.nan, 'NQ')
    
#     diff = diff[['accession', 'description', 'accession_res'] + mann_cols +  ['protein','sequence','gene_res'] + comb_median_cols + cv_cols + no_cols + rep_ratio_cols]  
    return sdf, df2, prots

def lpp_merge():
    sdf, df2, prots = df_reorder()
    p_med = 'cysteine_' + x_med
    lpp1 = lpp[['accession', 'gene_res', p_med]]
    lpp1 = lpp1.round(decimals = 1)
    lpp2 = lpp1.loc[lpp1[p_med] >= 2]
    lpp2 = pd.DataFrame(lpp2.groupby('accession')['gene_res'].apply(lambda row: ', '.join(row.astype(str)))).reset_index()
    lpp2.rename({'gene_res':'LPP-sensitive phosphosite'}, axis =1, inplace = True)
    sdf1 = sdf.merge(lpp2, on = 'accession', how = 'left')
    prots = prots.merge(lpp2, on = 'accession', how = 'left')
    mann_cols = [col for col in prots.columns if 'mann' in col]
    prots.set_index(['accession','description','protein','LPP-sensitive phosphosite']+mann_cols+['sequence','gene_res'], inplace = True)
    return sdf1, df2, prots



def main():
    sdf, df2, diff = lpp_merge()
    print('Writing to excel1')
    with pd.ExcelWriter('{}/supplemental_tables/{}_fig2_denature_combined.xlsx'.format(dir,date)) as writer:  # doctest: +SKIP
        # df2.to_excel(writer, sheet_name = 'filtered', index = False)
        # sdf.to_excel(writer, sheet_name = 'supplemental', index = False, float_format = '%.2f')
        diff.to_excel(writer, sheet_name = 'diff', float_format = '%.2f')
    # diff.to_csv('{}/supplemental_tables/{}_fig2_denature_combined.csv'.format(dir,date), float_format = '%.2f')
        
main()


