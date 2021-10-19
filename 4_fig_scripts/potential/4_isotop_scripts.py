#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 11:23:18 2021

@author: estherkemper
"""
import pandas as pd
import seaborn as sns, matplotlib.pyplot as plt
import os
plt.rcParams['pdf.fonttype'] = 42
sns.set_style('white')

dir = os.getcwd()
date = '20210830'

df = pd.read_excel('{}/supplemental_tables/{}_isotop_desalting_combined.xlsx'.format(dir,date))


def data_prep(gene_res):
    rep_cols = [col for col in df.columns if 'rep' in col]
    asynch_reps = [col for col in rep_cols if 'Asynch' in col]
    mitosis_reps = [col for col in rep_cols if 'Mitosis' in col]
    df1 = df.loc[df.gene_res == gene_res, rep_cols]
    # df1 = df1[['gene_res']+rep_cols]
    rename_rep_cols = [x.split('_', maxsplit = 1)[0] for x in rep_cols]
    df1.rename(columns = dict(zip(rep_cols, rename_rep_cols)), inplace = True)
    df2 = df1.melt().dropna(axis = 0)
    gene = gene_res.split('_')[0]
    res = 'C'+ gene_res.split('_')[1]
    title = gene + ' ' + res
    return df2, title

def plot_generes(gene_res):
    df, title = data_prep(gene_res)
    fig, axs = plt.subplots(figsize = (1.5,2.4)) #figsize=(2.8, 2.8)
    sns.stripplot(x= 'variable', y = 'value', data = df, jitter = True, dodge=True, ax = axs)
    # sns.boxplot(x = 'labels', y = 'data',  data = df, boxprops=dict(alpha=.5), \
    #             linewidth = 0.75, showfliers=False, ax = axs) 
    sns.barplot(x = 'variable', y = 'value', errwidth = 1, ci = 'sd', capsize=0.2, alpha = 0.7, data = df, ax = axs)
    plt.axhline(y=2, color = 'black', linestyle = 'dotted')
    plt.title(title.format(id), family = 'ARIAL', fontsize = 10)
    plt.ylabel('Gel-filtered/Unfiltered lysate\nMS1 ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xticks(fontsize = 10, family = 'ARIAL')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.savefig('{}/figures/fig2/{}_isotop_{}_bar.pdf'.format(dir, date, gene_res), bbox_inches='tight',transparent = True)

plot_generes('MAP2K4_246')
plot_generes('DTYMK_117')
