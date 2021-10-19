#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:31:47 2021

@author: estherkemper
"""
import os
import numpy as np
import pandas as pd
import seaborn as sns, matplotlib.pyplot as plt
sns.set(style="ticks")
from matplotlib import rcParams
import matplotlib.colors as mc
from matplotlib_venn import venn2
from  matplotlib.ticker import PercentFormatter

plt.rcParams['pdf.fonttype'] = 42
font = {'family' : 'ARIAL',
        'weight' : 'normal',
        'size'   : 10}         
plt.rc('font', **font) #set the font style created
plt.rcParams.update({'text.color' : "black"})

target_cutoff = 2
lpp_cutoff = 2
asynch_cutoff = 1.6
unchanging_cutoff = 1.6
protein_cutoff = 1.6
adapted_cutoff = 1.8

protein_cols = ['unique','accession','description', 'protein', 'gene_res']
phospho_col = ['pos_seq_mann', 'asynch_stoichiometry_mann','mitosis_stoichiometry_mann', 'stoichiometry','pos_seq_olsen', 'asynch_stoichiometry_olsen','mitosis_stoichiometry_olsen']

dir = os.getcwd()
#set a target cutoff defined in table
date = '20210830'

x_med ='cysteine_Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
x_prot_med ='protein_Mitosis_LPP(-,-)/(+,-)_combined_median'
y_prot_med = 'protein_Mitosis/Asynch_combined_median'
adapted_med ='cysteine_Mitosis_LPP(-,+)/(+,+)_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
df.columns = df.columns.str.replace('Mitosis_LPP','cysteine_Mitosis_LPP')

lpp_df = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')

phospho_df = pd.read_csv('{}/other_datasets/final_phospho_mann_olsen.csv'.format(dir))

df1 = df[['accession','protein','accession_res','gene_res',x_med, y_med]]
lpp1 = lpp_df[['accession', 'protein','accession_res','gene_res',x_med, y_med]]

def transform(df):
    df2 = df.loc[(df[x_med] >= lpp_cutoff) | (df[x_med] <= 1/lpp_cutoff)]
    df3 = df2.loc[(df2[y_med] > 1/asynch_cutoff) & (df2[y_med] < asynch_cutoff)]
    return df3
    
df2 = transform(df1)
lpp2 = transform(lpp1)

df3 = df2.loc[df2.accession.isin(set(lpp2.accession))]
lpp3 = lpp2.loc[lpp2.accession.isin(set(df2.accession))]

foo = pd.DataFrame(lpp2.groupby(['accession', 'protein'])['gene_res'].apply(list)).reset_index()
foo = lpp2.groupby(['accession', 'protein'])[x_med].apply(list)
