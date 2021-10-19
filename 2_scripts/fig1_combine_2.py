

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 12:12:48 2021

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
x_med = 'cysteine_Mitosis/Asynch_combined_median'
y_med = 'protein_Mitosis/Asynch_combined_median'
x_no = 'cysteine_Mitosis/Asynch_combined_no'
cat_col ='cysteine_Mitosis/Asynch_category'

expt_table = pd.read_csv('{}/inputs/protein_inputs.csv'.format(dir), index_col = None)


def combine():
    basal_df = pd.read_csv('{}/1_scripts/combined/{}_combined_basal.csv'.format(dir,date))
    protein = pd.read_csv('{}/1_scripts/combined/{}_combined_protein_whole_proteome.csv'.format(dir,date))
    protein.columns = ['protein_' + x if x not in protein_cols else x for x in protein.columns]
    basal_df.columns = ['cysteine_' + x if x not in protein_cols else x for x in basal_df.columns]
    denominator = 'Asynch_CTRL2'
    numerator = 'Mitosis_CTRL1'
    df = basal_df.merge(protein, how = 'outer', on = ['accession','description','protein'], suffixes = ['', ' protein'])
    rep_cols = [col for col in df.columns if '_channel_' in col]
    df.drop(rep_cols, axis = 1, inplace = True)
    df.columns = df.columns.str.replace('_ratio_','_')
    df.columns = df.columns.str.replace(numerator+'/'+denominator, 'Mitosis/Asynch')
    return df

# def basal_cv(row):
#     cat = row[cat_col]
#     if 'CV too high' in cat:
#         x_name = x_med.rsplit('_', maxsplit = 2)[0]
#         basal_exps = [x_name + '_' + s for s in basal_reps]
#         if row[basal_exps].count() >=2:
#             if (row[x_med] >= 2) or (row[x_med] <= 0.5):
#                 if (row[basal_exps] >= 2).all():
#                     return 'CV high but basal reps okay'
#                 elif (row[basal_exps] <= 0.5).all():
#                     return 'CV high but basal reps okay'
#                 else:
#                     return cat
#             else:
#                 return 'CV high and basal changing but together not'
#         else:
#             return cat
#     else:
#         return cat


def main_reorder():
    #reordering columns
    df1 = combine()
    keep_protein_cols = ['accession','description','protein','cysteine_sequence','accession_res', 'gene_res']
    category_cols = [col for col in df1.columns if 'category' in col]
    comb_median_cols = sorted([col for col in df1.columns if 'combined_median' in col])
    rep_ratio_cols = sorted([col for col in df1.columns if 'median_' in col])
    cv_cols = sorted([col for col in df1.columns if '_CV' in col])
    no_cols = sorted([col for col in df1.columns if '_no' in col])
    all_reps = sorted([col for col in df1.columns if 'channel' in col])
    # nunique_reps = sorted([col for col in df1.columns if 'no_unique' in col])
    # all_seq_reps = sorted([col for col in df1.columns if 'all_pep' in col])
    #select the columns we want in the order we want
    all_df = df1[keep_protein_cols + category_cols + comb_median_cols + rep_ratio_cols + cv_cols + no_cols + all_reps].copy()
    cond1 = all_df[x_no] >=2
    #quantified in 2x reps of basal
    df2 = all_df.loc[cond1]
    df2 = df2.sort_values(x_med, ascending = False)
    df3 = all_df.copy()
    df3[comb_median_cols + rep_ratio_cols + all_reps] = np.log2(df3[comb_median_cols + rep_ratio_cols + all_reps])
    df3[cat_col] = df3[cat_col].replace(np.nan, '0 or 1 rep', regex=True)
    # df3[cat_col] = df3.apply(basal_cv, axis = 1)
    df3.loc[(df3[cat_col].isnull()) | (df3[cat_col].str.contains('CV too high', na=False)), x_med] = 'NQ'
    # df3.loc[~ ((df3[cat_col].isnull()) | (df3[cat_col].str.contains('CV too high'))), x_med] = df3[x_med]
    df3.loc[(df3[x_med].isnull()) | (df3[x_no] <2), x_med] = 'NQ' 
    df3.loc[df3[y_med].isnull(), y_med] = 'NQ' 
    df3[x_med] = df3[x_med].replace('NQ', np.nan)
    df3 = df3.sort_values(x_med, ascending = False)
    df3[x_med] = df3[x_med].replace(np.nan, 'NQ')

    # df3 = df3.loc[~ ((df3[x_med] == 'NQ') & (df3[y_med] == 'NQ'))]

    # df3[x_med + '_filtered'] = df3.apply(test, axis = 1)
    # df2[cat_col] = df2.apply(basal_cv, axis = 1)
    df2 = df2.loc[~ ((df2[cat_col].isnull()) | (df2[cat_col].str.contains('CV too high', na=False)))] 

    # for col in category_cols:
    #     print(col)
    #     name = col.rsplit(sep = '_', maxsplit = 2)[0]
    #     med_col = [col for col in cv_med_cols if name in col][0]
    #     cat_col = [col for col in category_cols if name in col][0]
    #     df2.loc[(df2[cat_col].isnull()) | (df2[cat_col].str.contains('CV too high')), med_col] = np.nan
    #     df2.loc[~ ((df2[cat_col].isnull()) | (df2[cat_col].str.contains('CV too high'))), med_col] = df2[med_col]
    return all_df, df2, df3

def peptide_main():
    all_df, df2, df3 = main_reorder()
    drop_cols = [col for col in df3.columns if 'category' in col]
    df3.drop(drop_cols, axis =1, inplace = True)
    df3 = df3.assign(accession_res = df3.accession_res.str.rsplit('_', n = 1, expand = True))

    with pd.ExcelWriter('{}/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date)) as writer:  # doctest: +SKIP
        df2.to_excel(writer, sheet_name = 'filtered',  index = False)
        df3.to_excel(writer, sheet_name = 'supp_table', index = False, float_format = '%.2f')    
        # all_df.to_excel(writer, sheet_name = 'all', float_format='%.2f', index = False)
        # only_prot.to_excel(writer, sheet_name = 'only_protein', float_format='%.2f', index = False)
        # only_cys.to_excel(writer, sheet_name = 'only_cys', float_format='%.2f', index = False)
        # both.to_excel(writer, sheet_name = 'both', float_format='%.2f', index = False)
    return all_df, df2, df3
    
if __name__ == '__main__': 
    all_df, df2, df3 = peptide_main()



