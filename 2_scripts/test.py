

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

# protein_cols = ['accession','description', 'protein', 'gene_res']
protein_cols = ['unique','accession','description', 'protein', 'accession_res','gene_res']
# phospho_col = ['pos_seq_mann', 'asynch_stoichiometry_mann','mitosis_stoichiometry_mann', 'stoichiometry','pos_seq_olsen', 'asynch_stoichiometry_olsen','mitosis_stoichiometry_olsen']
x_med = 'cysteine_Mitosis/Asynch_combined_median'
y_med = 'protein_Mitosis/Asynch_combined_median'
x_no = 'cysteine_Mitosis/Asynch_combined_no'
basal_reps = ['Exp476', 'Exp483', 'Exp489', 'Exp494']
cat_col ='cysteine_Mitosis/Asynch_category'

expt_table = pd.read_csv('{}/inputs/protein_inputs.csv'.format(dir), index_col = None)


def combine():
    basal_df = pd.read_csv('{}/1_scripts/combined/{}_combined_basal.csv'.format(dir,date))
    protein = pd.read_csv('{}/1_scripts/combined/{}_combined_protein_whole_proteome.csv'.format(dir,date))
    denominator = 'Asynch_CTRL2'
    numerator = 'Mitosis_CTRL1'
    basal_keep_cols = [col for col in basal_df.columns if denominator not in col]
    basal_keep_cols = basal_keep_cols + [col for col in basal_df.columns if '/' + denominator in col]
    protein_keep_cols = [col for col in protein.columns if denominator not in col]
    protein_keep_cols = protein_keep_cols + [col for col in protein.columns if '/' + denominator in col]
    basal1_df = basal_df[basal_keep_cols]
    protein_df = protein[protein_keep_cols]
    basal1_df.columns = basal1_df.columns.str.replace('Mitosis_CTRL1/Asynch_CTRL2_ratio','cysteine_Mitosis/Asynch')
    protein_df.columns = protein_df.columns.str.replace('Mitosis_CTRL1/Asynch_CTRL2_ratio','protein_Mitosis/Asynch')
    basal1_df.columns = basal1_df.columns.str.replace(numerator,'cysteine_Mitosis/Asynch')
    protein_df.columns = protein_df.columns.str.replace(numerator,'protein_Mitosis/Asynch')
    df = basal1_df.merge(protein_df, how = 'outer', on = ['accession','description','protein'], suffixes = ['', ' protein'])
    rep_cols = [col for col in df.columns if '_rep_' in col]
    df.drop(rep_cols, axis = 1, inplace = True)
    return df

#don't think I want to merge with phospho for this dataset
# def phospho_combine(df):
#     phospho_df = pd.read_csv('{}/final_phospho_mann_olsen.csv'.format(dir))
#     phospho_col = [col for col in phospho_df.columns if 'accession' not in col]
#     df1 = df.merge(phospho_df, on = 'accession', how = 'left')
#     return df1, phospho_col

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

df3[df3.gene_res == 'RIF1_1692']
def main_reorder(df):
    #reordering columns
    df1 = combine()
    keep_protein_cols = ['accession','unique', 'description','protein','sequence','accession_res', 'gene_res']
    category_cols = [col for col in df1.columns if 'category' in col]
    comb_median_cols = sorted([col for col in df1.columns if 'combined_median' in col])
    rep_ratio_cols = sorted([col for col in df1.columns if 'median_Exp' in col])
    cv_cols = sorted([col for col in df1.columns if '_CV' in col])
    no_cols = sorted([col for col in df1.columns if '_no' in col])
    all_reps = sorted([col for col in df1.columns if 'rep' in col])
    # nunique_reps = sorted([col for col in df1.columns if 'no_unique' in col])
    # all_seq_reps = sorted([col for col in df1.columns if 'all_pep' in col])
    #select the columns we want in the order we want
    all_df = df1[keep_protein_cols + category_cols + comb_median_cols + rep_ratio_cols + cv_cols + no_cols + all_reps].copy()
    cond1 = all_df[x_no] >=2
    #quantified in 2x reps of basal
    df2 = all_df.loc[cond1]
    df2 = df2.sort_values(x_med)
    
    df3 = all_df.copy()
    df3[cat_col] = df3[cat_col].replace(np.nan, '0 or 1 rep', regex=True)
    # df3[cat_col] = df3.apply(basal_cv, axis = 1)
    df3.loc[(df3[cat_col].isnull()) | (df3[cat_col].str.contains('CV too high', na=False)), x_med] = 'NQ'
    # df3.loc[~ ((df3[cat_col].isnull()) | (df3[cat_col].str.contains('CV too high'))), x_med] = df3[x_med]
    df3.loc[df3[x_med].isnull(), x_med] = 'NQ' 
    df3.loc[df3[y_med].isnull(), y_med] = 'NQ' 
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

def peptide_main(df):
    all_df, df2, df3 = main_reorder(df)
    print(len(set(all_df.accession)))
    
    all_df.loc[all_df[x_no]<2, x_med] = np.nan
    all_df.loc[all_df[cat_col].str.contains('CV too high', na=False), x_med] = np.nan
    all_df1 = all_df.loc[~ ((all_df[x_med].isnull()) & (all_df[y_med].isnull())) ]
    
    #all data that passes criteria
    #take the proteins that aren't quantified in both, and set aside 
    both = all_df1.loc[(~all_df1[x_med].isnull()) & (~all_df1[y_med].isnull())]
    only = all_df1.loc[~(all_df1.accession.isin(set(both.accession)))]
    cys_cond = (only[y_med].isnull()) & (~only[x_med].isnull())
    prot_cond = (only[x_med].isnull()) & (~only[y_med].isnull())
    
    only_prot = only.loc[prot_cond]
    only_cys = only.loc[cys_cond]
    data = len(set(only_cys.accession)), len(set(only_prot.accession)), len(set(both.accession))
    labels = ['cys', 'prot','both'] 
    s1_venn = pd.DataFrame(data = {'labels': labels, 'cys':data})
    s1_venn.to_csv('{}/2_scripts/supplemental_tables/figs1_venn.csv'.format(dir,date), index = False)
    print(s1_venn)
    
    with pd.ExcelWriter('{}/2_scripts/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date)) as writer:  # doctest: +SKIP
        df2.to_excel(writer, sheet_name = 'filtered', float_format='%.2f', index = False)
        df3.to_excel(writer, sheet_name = 'supp_table', float_format='%.2f', index = False)    
        # all_df.to_excel(writer, sheet_name = 'all', float_format='%.2f', index = False)
        # only_prot.to_excel(writer, sheet_name = 'only_protein', float_format='%.2f', index = False)
        # only_cys.to_excel(writer, sheet_name = 'only_cys', float_format='%.2f', index = False)
        # both.to_excel(writer, sheet_name = 'both', float_format='%.2f', index = False)
    return all_df, df2, df3
    
if __name__ == '__main__': 
    df = combine()
    all_df, df2, df3 = peptide_main(df)


old = pd.read_excel('{}/2_scripts/supplemental_tables/old_supplementary_tables.xlsx'.format(dir), sheet_name = 'Table1. Asynch vs Mitosis', header =2)
old_generes = old.gene_res
new_generes = df3.gene_res

not_in_old = list(set(new_generes) - set(old_generes))
not_in_new = list(set(old_generes) - set(new_generes))

only_new_df = df3[df3.gene_res.isin(not_in_old)]
only_old_df = old[old.gene_res.isin(not_in_new)]
filtered_only_new = df2[df2.gene_res.isin(not_in_old)]
