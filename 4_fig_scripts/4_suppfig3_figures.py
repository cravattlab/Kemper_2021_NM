#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 11:36:26 2021

@author: ekk
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns, matplotlib.pyplot as plt
sns.set(style="ticks")
from matplotlib import rcParams
import matplotlib.colors as mc
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
adapted_cutoff = 1.8
protein_cutoff = 1.6

protein_cols = ['unique','accession','description', 'protein', 'gene_res']
phospho_col = ['pos_seq_mann', 'asynch_stoichiometry_mann','mitosis_stoichiometry_mann', 'stoichiometry','pos_seq_olsen', 'asynch_stoichiometry_olsen','mitosis_stoichiometry_olsen']

dir = os.getcwd()
date = '20210830'

x_med ='Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
adapted_med ='Mitosis_LPP(-,+)/(+,+)_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir,date), sheet_name = 'filtered')
lpp_df = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
print('Finished reading input data')

def color_prep(df):
    if (df[x_med].item() <= 1/lpp_cutoff) and (df[y_med].item() <= 1/lpp_cutoff):
        if df[adapted_med].item() <= 1/adapted_cutoff:
            colors = ['red','red','red']
            return colors
        elif df[adapted_med].item() >= adapted_cutoff:
            colors = ['red','red','blue']
            return colors
        else:
            colors = ['red','red','#404040']
            return colors
    elif (df[x_med].item() >= lpp_cutoff):
        if (df[y_med].item() >= lpp_cutoff) and (df[adapted_med].item() >= adapted_cutoff):
            colors = ['blue','blue','blue']
            return colors
    elif (df[x_med].item() < lpp_cutoff) & (df[x_med].item() > 1/lpp_cutoff):
        if (df[y_med].item() < lpp_cutoff) & (df[y_med].item() > 1/lpp_cutoff):
            if (df[adapted_med].item() >= adapted_cutoff):
                colors = ['black','black','blue']
                return colors
            
def supp_fig3ae_boxplots(gene_cys):
    gene1 = df.loc[df.gene_res == gene_cys]
    sdf = gene1.copy()
    drop_cols = [col for col in sdf.columns if 'limited' in col] + [col for col in sdf.columns if 'category' in col] + [col for col in sdf.columns if 'protein_' in col]
    sdf = sdf.drop(drop_cols, axis = 1)
    sdf.set_index(['accession', 'description', 'protein','asynch_stoichiometry_mann','mitosis_stoichiometry_mann','sequence','gene_res'], inplace = True)
    median_cols = [col for col in sdf.columns if 'median' in col]
    sdf[median_cols] = np.log2(sdf[median_cols])
    gene_res = gene_cys.split('_')[0] + ' C' + gene_cys.split('_')[1]
    sdf.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig3be_box_{}.csv'.format(dir, gene_cys), float_format = '%.2f')
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    adapted_name = adapted_med.rsplit(sep = '_', maxsplit = 2)[0]
    adapted_cols = [col for col in gene1.columns if adapted_name + '_median_TMT' in col] 
    colors = color_prep(gene1)
    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] +adapted_cols +\
        [col for col in gene1.columns if y_name + '_median_TMT' in col] 
    rep_cols = [col for col in rep_cols if 'protein_' not in col] 

    gene = gene1[['gene_res'] + rep_cols]    
    # gene[rep_cols].columns = gene[rep_cols].columns.str.split('_median_TMT_rep', expand = True)
    gene = gene.melt(id_vars = 'gene_res').dropna(axis = 0)
    gene = gene.assign(variable1 = gene.variable.str.rsplit(pat = '_', n = 3, expand = True)[0])
    gene = gene.assign(value = np.log2(gene.value))
    gene = gene.sort_values('variable1', ascending = True)
    
    fig, axs = plt.subplots( figsize=(1.2, 2.4))
    error_kw=dict(lw=1, capsize=3, capthick=1)
    sns.boxplot(x = 'variable1', y = 'value',  data = gene, palette = colors, \
                order = [y_name, x_name, adapted_name], \
                boxprops=dict(alpha=.8), linewidth = 0.75, showfliers=False) 
    sns.stripplot(x = 'variable1', y = 'value', jitter = 0.2, linewidth = 0.5,  edgecolor = 'black',data = gene, palette = colors, order =  [y_name, x_name, adapted_name])
    fig = plt.gcf()
    plt.ylim(-4.5,4.5)
    plt.ylabel('log2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.xticks(ticks = [0,1,2], labels = ['Mitosis/Asynch','LPP(-, -)/(+, -)',' LPP(-, +)/(+, +)'], rotation = 40, ha = 'right', fontsize = 10, family = 'ARIAL')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.title(gene_res, family = 'ARIAL', fontsize = 10)
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.savefig('{}/figures/suppfig3/{}_suppfig3ae_bar_{}.pdf'.format(dir,date, gene_cys), bbox_inches='tight',transparent = True)
    plt.close()
cys_list = ['EDC3_137', 'GTF2I_215', 'SLAIN2_152']
for cys in cys_list:
    supp_fig3ae_boxplots(cys)
    
def supp_fig3bf_bar(df, gene_res,sub):
    gene1 = df.loc[df.gene_res == gene_res]
    x_med1 = 'cysteine_' + x_med
    sdf = gene1.copy()
    drop_cols = [col for col in sdf.columns if 'category' in col] + [col for col in sdf.columns if 'protein_' in col]
    sdf.drop(drop_cols, axis =1, inplace = True)
    sdf.columns = sdf.columns.str.replace('cysteine_','peptide_')
    median_cols = [col for col in sdf.columns if 'median' in col]
    sdf[median_cols] = 1/sdf[median_cols]
    print(sdf[median_cols])
    colors = ['m','m','m']
    sdf.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig3bf_bar_{}.csv'.format(dir, gene_res), index = False, float_format = '%.2f')
    x_name = x_med1.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] +\
        [col for col in gene1.columns if y_name + '_median_TMT' in col] 
    rep_cols = [col for col in rep_cols if 'protein_' not in col] 
    location = (gene1['sequence'].str.find('*', 0) - 1).item()
    aa = gene1.sequence.item()[location]
    res = aa + gene1.gene_res.str.split('_', n =1, expand = True)[1]
    gene_res1 = (gene_res.split('_')[0] + ' ' + res).item()
    gene = gene1[['gene_res'] + rep_cols]    
    gene.columns = gene.columns.str.replace('_median_TMT_rep*', '')
    gene = gene.melt(id_vars = 'gene_res').dropna(axis = 0)
    gene = gene.assign(variable1 = gene.variable.str.rsplit(pat = '_', n = 1, expand = True)[0])
    len_x = len([x for x in gene.variable if x_name in x])
    len_y = len([x for x in gene.variable if y_name in x])
    count_no = np.min([len_x, len_y])
    count_append = pd.DataFrame({'gene_res': gene_res, \
                                 'variable1':'Mitosis LPP(-)', \
                                'value': 1}, index = [np.arange(0,count_no,1)])
    gene = gene.append(count_append).reset_index(drop = True)
    gene['value1'] = 1/gene['value']
    print(gene.value1)
    order = [ y_name, 'Mitosis LPP(-)', x_name]
    fig, axs = plt.subplots(figsize = (1.5,2.4)) #figsize=(2.8, 2.8)
    sns.swarmplot(x= 'variable1', y = 'value1', data = gene, linewidth = 0.5, order = order, dodge=True, ax = axs, palette = colors, edgecolor = 'black')

    sns.barplot(x = 'variable1', y = 'value1', errwidth = 1, ci = 'sd', order = order, capsize=0.3, alpha = 0.7, data = gene, ax = axs, palette = colors)
    plt.axhline(y=0.5, color = 'black', linestyle = 'dotted')
    plt.ylabel('Percent of LPP(-)\nMS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xticks(ticks = [0,1,2], labels = ['Asynch LPP(-)','Mitosis LPP(-)', 'Mitosis LPP(+)'], rotation = 40, ha = 'right', fontsize = 10, family = 'ARIAL')
    plt.yticks(ticks = np.arange(0,1.2,0.2), fontsize = 10, family = 'ARIAL')
    plt.title(gene_res1, family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.savefig('{}/figures/suppfig3/{}_suppfig3bf_phospho_{}.pdf'.format(dir, date, gene_res1, sub), bbox_inches='tight',transparent = True)
    plt.close()
for gene in ['EDC3_131', 'SLAIN2_147','GTF2I_210']:
    supp_fig3bf_bar(lpp_df, gene, 'phospho')
    
def supp_fig3d_boxplots(protein_name):
    gene1 = df.loc[df.protein == protein_name]
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    adapted_name = adapted_med.rsplit(sep = '_', maxsplit = 2)[0]
    if len(set(gene1.gene_res)) > 9:
        x_no = x_name + '_combined_no'
        y_no = y_name + '_combined_no'
        gene1 = gene1.loc[(gene1[x_no]>=8) & (gene1[y_no] >= 8)]
    else:
        pass
    sdf = gene1.copy()
    drop_cols = [col for col in sdf.columns if 'category' in col] + [col for col in sdf.columns if 'protein_' in col]
    sdf.drop(drop_cols, axis =1, inplace = True)
    sdf.set_index(['accession', 'description', 'protein','asynch_stoichiometry_mann','mitosis_stoichiometry_mann','sequence','accession_res','gene_res'], inplace = True)
    median_cols = [col for col in sdf.columns if 'median' in col]
    sdf[median_cols] = np.log2(sdf[median_cols])
    sdf.to_excel('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig3d_FLNA.xlsx'.format(dir), float_format = '%.2f')
    adapted_cols = [col for col in gene1.columns if adapted_name + '_median_TMT' in col] 
    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] +adapted_cols +\
        [col for col in gene1.columns if y_name + '_median_TMT' in col] 
    gene = gene1[['gene_res'] + rep_cols]    
    gene = gene.melt(id_vars = 'gene_res').dropna(axis = 0)
    gene = gene.assign(pos = gene.gene_res.str.split('_', n =1, expand = True)[1].astype(int))
    gene = gene.sort_values('pos', ascending = True)
    gene = gene.assign(variable = gene.variable.str.rsplit(pat = '_', n = 3, expand = True)[0])
    gene = gene.assign(value = np.log2(gene.value))
    gene = gene.assign(gene_res = 'C' + gene.gene_res.str.split('_', n =1, expand = True)[1])
    fig, axs = plt.subplots( figsize=(10, 2.4))
    order = [y_name, x_name, adapted_name]
    o = sns.color_palette()[1]
    # sns.stripplot(x= 'gene_res', y = 'value', hue = 'variable', alpha = 0.7, linewidth = 0.5, edgecolor = 'gray', hue_order = order, data = gene, jitter = True, dodge=True, palette = ['b','g','y'], ax = axs)
    sns.stripplot(x = 'gene_res', y = 'value', hue = 'variable', linewidth = 0.5, edgecolor = 'gray',\
                  data = gene, jitter = 0.2, dodge = True, palette = ['b','g','#daa520'], hue_order = order, ax = axs)
    sns.boxplot(x = 'gene_res', y = 'value', hue = 'variable', data = gene, palette = ['b','g','#daa520'], \
                hue_order = order, boxprops=dict(alpha=.7), linewidth = 0.75, showfliers=False, ax = axs) 
    plt.ylabel('log2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    handles, labels = axs.get_legend_handles_labels() 
    l = plt.legend(handles[0:3], ['Mitosis/Asynch','LPP(-, -)/LPP(+, -)','LPP(-, +)/LPP(+, +)'], bbox_to_anchor=(1, 0.75), loc=2, fontsize = 10, frameon = False)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 45, ha = 'right')
    plt.title(protein_name, family = 'ARIAL', fontsize = 10)
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.savefig('{}/figures/suppfig3/{}_suppfig3d_bar_{}.pdf'.format(dir,date, protein_name), bbox_inches='tight',transparent = True)
    plt.close()
    
supp_fig3d_boxplots('FLNA')

