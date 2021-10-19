#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 08:48:30 2021

@author: ekk
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns, matplotlib.pyplot as plt
sns.set(style="ticks")
from matplotlib import rcParams
import matplotlib.colors as mc
from matplotlib_venn import venn2

plt.rcParams['pdf.fonttype'] = 42
font = {'family' : 'ARIAL',
        'weight' : 'normal',
        'size'   : 10}         
plt.rc('font', **font) #set the font style created
plt.rcParams.update({'text.color' : "black"})

target_cutoff = 2
lpp_cutoff = 2
asynch_cutoff = 1.6
unchanging_cutoff = 1.5
protein_cutoff = 1.6


dir = os.getcwd()
date = '20210830'
x_med = 'cysteine_Mitosis/Asynch_combined_median'
y_med = 'protein_Mitosis/Asynch_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date), sheet_name = 'filtered')
supp_df = pd.read_excel('{}/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date), sheet_name = 'supp_table')
desalt_df = pd.read_excel('{}/supplemental_tables/{}_isotop_desalting_combined.xlsx'.format(dir,date))

def suppfig1be_boxplots(protein_name, width):
    gene1 = df.loc[df.protein == protein_name]
    gene1.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig1be_{}.csv'.format(dir, protein_name), index = False, float_format = '%.2f')
    gene1.loc[(gene1[x_med] >= 2), 'hue'] = 'blue'
    gene1.loc[(gene1[x_med] <= 0.5), 'hue'] = 'red'
    gene1.loc[(gene1[x_med] > 0.5) & (gene1[x_med] < 2), 'hue'] = '#404040'
    gene_name = gene1.iloc[0,:].gene_res.split('_')[0]
    gene1 = gene1.assign(pos = gene1.gene_res.str.split('_', expand = True)[1].astype(int))
    gene1 = gene1.sort_values('pos', ascending = True)     
    hue_order = list(gene1.hue)
    order = list('C' + gene1.pos.astype(str))     
    rep_cols = [col for col in gene1.columns if '_median_TMT' in col]
    rep_cols = [col for col in rep_cols if 'protein_' not in col]
    gene = gene1[['gene_res'] + rep_cols]
    gene = gene.melt(id_vars = 'gene_res').dropna(axis = 0)
    gene = gene.assign(variable = gene.variable.str.rsplit(pat = '_', n = 1, expand = True)[0])
    gene = gene.assign(value = np.log2(gene.value))
    gene = gene.assign(gene_res = 'C' + gene.gene_res.str.split('_', n =1, expand = True)[1])
    fig, axs = plt.subplots(figsize = (width, 2.4)) #figsize=(2.8, 2.8)
    error_kw=dict(lw=1, capsize=5, capthick=1)
    sns.boxplot(x = 'gene_res', y = 'value',  palette = hue_order, order = order, data = gene, \
                boxprops=dict(alpha=.7), showfliers=False, linewidth = 0.75,  ax = axs)
    sns.stripplot(x= 'gene_res', y = 'value', palette = hue_order, order = order, linewidth = 0.5,\
                  edgecolor = 'black', data = gene, jitter = 0.3, dodge=True,  ax = axs)
    plt.ylim(-4.5,4.5)
    plt.title(gene_name, fontsize = 10)
    handles, labels = axs.get_legend_handles_labels() 
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.ylabel('Mitosis/Asynch\nlog2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.title(gene_name, family = 'ARIAL', fontsize = 10)
    plt.savefig('{}/figures/suppfig1/{}_suppfig1be_box_{}.pdf'.format(dir, date, gene_name), bbox_inches='tight',transparent = True)
    plt.close()
suppfig1be_boxplots('GAPDH', 1.2) 
suppfig1be_boxplots('PARK7', 1.2) 
suppfig1be_boxplots('DTYMK', 1.2) 

def supp_fig1d_lpp_venn(df):
    x_no = 'cysteine_Mitosis/Asynch_combined_no'
    cat_col ='cysteine_Mitosis/Asynch_category'
    df[[x_med, y_med]] = df[[x_med, y_med]].replace('NQ', np.nan, regex = True)
    df.loc[df[x_no]<2, x_med] = np.nan
    df.loc[df[cat_col].str.contains('CV too high', na=False), x_med] = np.nan
    #all data that passes criteria
    #take the proteins that aren't quantified in both, and set aside 
    protein_cols = ['accession','protein','description','gene_res']
    keep_cols = [col for col in df.columns if 'combined_median' in col]
    sdf = df[protein_cols + keep_cols ]
    sdf = sdf.loc[~ ((sdf[x_med].isnull()) & (sdf[y_med].isnull())) ]
    both = sdf.loc[(~sdf[x_med].isnull()) & (~sdf[y_med].isnull())]
    cys = sdf.loc[(sdf[y_med].isnull()) & (~sdf[x_med].isnull())] 
    prot = sdf.loc[(sdf[x_med].isnull()) & (~sdf[y_med].isnull())] 
    sdf[[y_med]] = sdf[[y_med]].fillna('NQ')
    sdf['ABPP Gene_res: Mitosis/Asynch'] = sdf.gene_res + ': ' + sdf[x_med].astype(str)
    sdf.loc[sdf['ABPP Gene_res: Mitosis/Asynch'].str.contains('nan', regex = True, na=False)] = np.nan
    new_sdf = pd.DataFrame(sdf.groupby('accession')[['protein','description',y_med]].first())
    new_sdf['ABPP Gene_res: Mitosis/Asynch'] = pd.Series(sdf.groupby('accession')['ABPP Gene_res: Mitosis/Asynch'].apply(list))
    new_sdf = new_sdf.reset_index()
    a = new_sdf.accession.isin(list(set(cys.accession))) 
    b = new_sdf.accession.isin(list(set(prot.accession))) 
    c = new_sdf.accession.isin(list(set(both.accession))) 
    new_sdf.loc[ a,'UNENRICHED PROTEOMICS OR TMT-ABPP'] = 'TMT-ABPP'
    new_sdf.loc[b,'UNENRICHED PROTEOMICS OR TMT-ABPP'] = 'UNENRICHED PROTEOMICS'
    new_sdf.loc[((a & b) | c),'UNENRICHED PROTEOMICS OR TMT-ABPP'] = 'BOTH'
    new_sdf.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig1c.csv'.format(dir), index = False, float_format = '%.2f')
    df1 = df.loc[~ ((df[x_med].isnull()) & (df[y_med].isnull())) ]
    both = df1.loc[(~df1[x_med].isnull()) & (~df1[y_med].isnull())]
    only = df1.loc[~(df1.accession.isin(set(both.accession)))]
    cys_cond = (only[y_med].isnull()) & (~only[x_med].isnull())
    prot_cond = (only[x_med].isnull()) & (~only[y_med].isnull())
    only_prot = only.loc[prot_cond]
    only_cys = only.loc[cys_cond]
    data = len(set(only_cys.accession)), len(set(only_prot.accession)), len(set(both.accession))
    labels = ['cys', 'prot','both'] 
    s1_venn = pd.DataFrame(data = {'labels': labels, 'cys':data})
    #source data an venn datais generated in 2_fig1_combine.py
    sizes = s1_venn['cys']
    #required in 2x LPP/CTRL, 1x Asynch/Mitosis
    #for 2x for everyone, 2x rep in asynch/mitosis and 2x in lpp
    plt.figure(figsize=(3.2, 3.2))
    out = venn2(subsets = sizes, set_labels = ('Cysteine reactivity', 'Protein expression'), set_colors = ('gray','b'), alpha = 0.9 )
    plt.title('Proteins quantified in\nTMT-ABPP vs unenriched proteomics', fontsize = 10)
    for text in out.set_labels:
        text.set_fontsize(10)
    for text in out.subset_labels:
        text.set_fontsize(10)
        plt.savefig('{}/figures/suppfig1/{}_suppfig1d_venn.pdf'.format(dir, date), bbox_inches='tight')
supp_fig1d_lpp_venn(supp_df)

def supp_fig1e_GO_chart(sub, type):
    #input data generated in 3_fig1_proteins
    #then put through webgestalt, cutoff 0.05
    #then put through revigo to reduce redundancy - this is source data
    df = pd.read_csv('{}/for_figures/suppfig1/GO_analysis/suppfig1_revigo_{}_{}_simrei_05.csv'.format(dir, sub, type))
    df.columns = df.columns.str.strip(' ')
    df['Eliminated'] = df['Eliminated'].str.strip(' ')
    df['Name'] = df['Name'].str.strip(' ')
    df['Name'] = df['Name'].str.strip('"')
    df['p-value'] = -df['Value']
    rev_df = df.loc[df.Eliminated == 'False']
    rev_df = rev_df.sort_values('p-value', ascending = False)
    df1 = rev_df.iloc[0:15,:]
    df1 = df1.sort_values('p-value', ascending = True)
    x = df1['p-value']
    y = df1['Name']
    fig, axs = plt.subplots(figsize = (3, 2.4)) #figsize=(2.8, 2.8)
    plt.barh(y,x, color='black')
    plt.xticks(fontsize = 10, family = 'ARIAL')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.xlabel("Adjusted P-value -log10",fontsize = 10, family = 'ARIAL')
    plt.title('GO ' + type + ' process', fontsize = 10, family = 'ARIAL')
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig('{}/figures/suppfig1/{}_suppfig1e_{}_{}.pdf'.format(dir, date, sub, type), bbox_inches='tight')
    plt.close()
supp_fig1e_GO_chart('reactivity','cellular')

def supp_fig1h_protein(protein):
    df1 = desalt_df.loc[desalt_df.protein == protein]
    df1.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig1h_protein_{}.csv'.format(dir, protein))
    rep_cols = [col for col in df1.columns if 'rep' in col]
    rep_cols = [col for col in rep_cols if 'Mitosis' in col]
    df1 = df1[['gene_res'] + rep_cols]
    df1 = df1.melt('gene_res').dropna(axis = 0)
    df1['pos'] = df1.gene_res.str.split('_', expand = True)[1].astype(int)
    df1['res'] = 'C' + df1.pos.astype(str)
    df1 = df1.sort_values('pos')
    df2 = df1[['res'] + ['value']]
    print(df2)
    width = len(set(df2.res)) * 0.5
    fig, axs = plt.subplots(figsize = (width,2.2)) #figsize=(2.8, 2.8)
    # color = [(0.8666666666666667, 0.5176470588235295, 0.3215686274509804)] * len(set(df.gene_res))
    color = [['b'] * len(set(df2.res))][0]
    sns.swarmplot(x= 'res', y = 'value', palette = color, data = df2, linewidth = 0.5, edgecolor = 'black',  dodge=True, ax = axs)#jitter = True,
    sns.barplot(x = 'res', y = 'value', palette = color, ci = 'sd', errwidth = 1, capsize=0.2, alpha = 0.7, data = df2, ax = axs)
    plt.axhline(y=2, color = 'black', linestyle = 'dotted')
    plt.ylabel('Gel-filtered/Unfiltered lysate\nMS1 ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.title(protein, fontsize = 10)
    plt.ylim([0,7])
    plt.savefig('{}/figures/suppfig1/{}_suppfig1h_desalt_protein_{}.pdf'.format(dir, date, protein), bbox_inches='tight',transparent = True)  
    plt.close()
supp_fig1h_protein('DTYMK')

def supp_fig1i_prep(gene_res):
    rep_cols = [col for col in desalt_df.columns if 'rep' in col]
    sdf = desalt_df.loc[desalt_df.gene_res == gene_res]
    sdf.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig1i_{}.csv'.format(dir, gene_res), float_format = '%.2f')
    df1 = desalt_df.loc[desalt_df.gene_res == gene_res, rep_cols]
    rename_rep_cols = [x.split('_', maxsplit = 1)[0] for x in rep_cols]
    df1.rename(columns = dict(zip(rep_cols, rename_rep_cols)), inplace = True)
    df2 = df1.melt().dropna(axis = 0)
    gene = gene_res.split('_')[0]
    res = 'C'+ gene_res.split('_')[1]
    title = gene + ' ' + res
    return df2, title

def supp_fig1i_bar(gene_res):
    df, title = supp_fig1i_prep(gene_res)
    fig, axs = plt.subplots(figsize = (1.1,2.2)) #figsize=(2.8, 2.8)
    color = ['#404040','b']
    sns.swarmplot(x= 'variable', y = 'value', data = df, linewidth = 0.5, edgecolor = 'black', palette = color, dodge=True, ax = axs)
    sns.barplot(x = 'variable', y = 'value', errwidth = 1, ci = 'sd', capsize=0.2, alpha = 0.7, data = df, palette = color, ax = axs)
    plt.axhline(y=2, color = 'black', linestyle = 'dotted')
    plt.title(title.format(id), family = 'ARIAL', fontsize = 10)
    plt.ylabel('Gel-filtered/Unfiltered lysate\nMS1 ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.ylim([0, 7])
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.savefig('{}/figures/suppfig1/{}_suppfig1i_desalt_cys_{}_bar.pdf'.format(dir, date, gene_res), bbox_inches='tight',transparent = True)
    plt.close()
supp_fig1i_bar('DTYMK_117')
   
def supp_fig1j_protein_boxplots(protein_list, width):
    gene1 = df.loc[df.protein.isin(protein_list)]
    rep_cols = [col for col in gene1.columns if  '_median_TMT' in col]
    rep_cols = [col for col in rep_cols if 'protein_' in col]
    # gene = gene1[['protein'] + rep_cols].dropna(subset = rep_cols, how = 'all')
    gene = gene1[['protein'] + rep_cols]
    gene = gene.groupby('protein').first().reset_index()
    gene[rep_cols] = np.log2(gene[rep_cols])
    gene = gene.melt(id_vars = 'protein')
    fig, axs = plt.subplots(figsize = (width, 2.4)) #figsize=(2.8, 2.8)
    sns.boxplot(x = 'protein', y = 'value',  palette = ['green'], data = gene, \
                boxprops=dict(alpha=.7), showfliers=False, linewidth = 0.75,  ax = axs)
    sns.stripplot(x= 'protein', y = 'value', palette = ['green'], marker = '^',\
                  linewidth = 0.5, data = gene, jitter = 0.3, dodge=True,  ax = axs)
    plt.ylim(-4.5,4.5)
    handles, labels = axs.get_legend_handles_labels() 
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.ylabel('Mitosis/Asynch\nlog2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.title('Protein expression', family = 'ARIAL', fontsize = 10)
    plt.savefig('{}/figures/suppfig1/{}_suppfig1j_protein_bar.pdf'.format(dir, date), bbox_inches='tight',transparent = True)
    plt.close()
protein_list = ['DTWD1', 'NOL8', 'RRP15']
supp_fig1j_protein_boxplots(protein_list, 1.5)