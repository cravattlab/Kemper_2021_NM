#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os

dir = os.getcwd()
target_cutoff = 2.5
protein_cols = ['accession','description', 'protein', 'sequence', 'gene_res']
inputs = '/inputs/isotop_peptide_inputs.csv'
date = '20210830'
expt_table = pd.read_csv('{}/{}'.format(dir, inputs))
x_name = 'Mitosis_Gel-filtered/Unfiltered_'
y_name = 'Asynch_Gel-filtered/Unfiltered_'

def sequence_combine(row):
     reps = [col for col in row.index if 'sequence' in col]
     sequences = row[reps].dropna()
     shortest_sequence =  min(sequences, key=len)
     return pd.Series(shortest_sequence)

def main():
    #create zip of expt_code to names of expt_replicates
    all_groups = list(set(expt_table.reps.str.rsplit( '_', n = 1, expand = True)[0]))

    multi = pd.DataFrame(columns=protein_cols)
    for group in all_groups:
        df = pd.read_csv('{}/1_scripts/combined_isotop/{}_{}.csv'.format(dir,date, group))
        # df = df[protein_cols + ['combined_median', 'combined_CV', 'no']]
        # df.set_index(protein_cols, inplace = True)
        # data_cols = list(df.columns)
        # df.columns = [group +  '_' + s for s in data_cols]
        # df.reset_index(inplace = True)
        multi = multi.merge(df, on= protein_cols, how='outer')
        combined_median_cols = [col for col in multi.columns if 'combined_median' in col]
        std_cols = [col for col in multi.columns if 'combined_CV' in col]
        cat_cols = [col for col in multi.columns if 'category' in col]
        no_cols = [col for col in multi.columns if 'combined_no' in col]
        rep_cols = sorted([col for col in multi.columns if '_rep' in col])
        sequence_cols = [col for col in multi.columns if 'sequence' in col]
        multi['sequence'] = multi[sequence_cols].apply(sequence_combine, axis = 1)
        multi = multi[protein_cols + cat_cols + combined_median_cols + std_cols + no_cols+ rep_cols]
    multi.dropna(subset = combined_median_cols, how = 'all', axis = 0, inplace = True)
    multi.columns = multi.columns.str.replace('Gel-filtered_over_unfiltered', 'Gel-filtered/Unfiltered')
    multi.loc[multi[x_name + 'combined_no'] < 2, x_name + 'combined_median'] = np.nan
    multi.loc[multi[y_name + 'combined_no'] < 2, y_name + 'combined_median'] = np.nan
    multi.loc[multi[x_name + 'category'].str.contains('CV too high,', na=False), x_name + 'combined_median'] = np.nan
    multi.loc[multi[y_name + 'category'].str.contains('CV too high,', na=False), y_name + 'combined_median'] = np.nan
    multi = multi.sort_values(x_name+'combined_median', ascending = False)
    combined_cols = [col for col in multi.columns if 'combined_median' in col]
    cat_cols = [col for col in multi.columns if 'category' in col]
    multi.drop(cat_cols, axis = 1, inplace = True)
    diff = multi.loc[(multi[x_name + 'combined_median'] >= 2) | (multi[x_name + 'combined_median'] <= 0.5)]
    multi[combined_cols] = multi[combined_cols].replace(np.nan, 'NQ')
    diff[combined_cols] = diff[combined_cols].replace(np.nan, 'NQ')

    with pd.ExcelWriter('{}/supplemental_tables/{}_isotop_desalting_combined.xlsx'.format(dir, date)) as writer:  # doctest: +SKIP
        # protein.to_excel(writer, sheet_name = 'protein', index = False, float_format = '%.2f')
        multi.to_excel(writer, sheet_name = 'all', index = False, float_format = '%.2f')
        diff.to_excel(writer, sheet_name = 'diff', index = False, float_format='%.2f')

if __name__=="__main__":
    main()