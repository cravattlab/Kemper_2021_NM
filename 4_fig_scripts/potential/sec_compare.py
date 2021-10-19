#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 10:44:46 2021

@author: estherkemper
"""
import os
import numpy as np
import pandas as pd

lpp_cutoff = 2
asynch_cutoff = 1.6
unchanging_cutoff = 1.6
adapted_cutoff = 1.8

protein_cols = ['unique','accession','description', 'protein', 'gene_res']
phospho_col = ['asynch_stoichiometry_mann','mitosis_stoichiometry_mann', 'asynch_stoichiometry_olsen','mitosis_stoichiometry_olsen']

dir = os.getcwd()
date = '20210830'

x_med ='Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
adapted_med ='Mitosis_LPP(-,+)/(+,+)_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
# prots_df = pd.read_excel('{}/supplemental_tables/{}_fig2_proteins.xlsx'.format(dir,date), sheet_name = 'all prots with changes')
# y_df = pd.read_excel('{}/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date), sheet_name = 'filtered')
# y_exp = pd.read_excel('{}/supplemental_tables/{}_fig1_proteins.xlsx'.format(dir,date), sheet_name = 'all prots with changes')

sec = pd.read_csv('{}/other_datasets/aebersold_sec.csv'.format(dir))
sec = sec.dropna(axis = 1)

rep_cols = [col for col in df.columns if 'median_TMT' in col]
cat_cols = [col for col in df.columns if 'category' in col]
lim_cols = [col for col in df.columns if 'limited' in col]
phos_cols = [col for col in df.columns if 'stoichiometry' in col]
drop_cols = rep_cols + cat_cols + lim_cols + phos_cols

df1 = df.drop(drop_cols, axis = 1)

lpp_cond = (df1[x_med] >= lpp_cutoff) | (df1[x_med] <= 1/lpp_cutoff)
cond2 = (df1[y_med] <= (1/asynch_cutoff)) |  (df1[y_med] >= asynch_cutoff)
#peptides that are changing with LPP and asynch but not artifactual changes
changing = df1.loc[(lpp_cond) & (cond2)]
    


total_df = changing.merge(sec, how = 'inner', left_on = 'accession', right_on = 'protein_id')

sec_change = total_df.loc[total_df.SEC_change == True]
sec_same = total_df.loc[total_df.SEC_change == False]


all =df1.merge(sec, how = 'inner', left_on = 'accession', right_on = 'protein_id')
all_change = all.loc[all.SEC_change == True]
all_same = all.loc[all.SEC_change == False]

len(set(sec_change.accession_res))
len(set(sec_same.accession_res))

len(set(all_change.accession_res))
len(set(all_same.accession_res))

total_df.to_csv('{}/figures/suppfig2/potential/sec_merge.csv'.format(dir))
