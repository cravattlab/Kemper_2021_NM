#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 13:09:04 2021

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

protein_cols = ['unique','accession','description', 'protein', 'gene_res']
phospho_col = ['pos_seq_mann', 'asynch_stoichiometry_mann','mitosis_stoichiometry_mann', 'stoichiometry','pos_seq_olsen', 'asynch_stoichiometry_olsen','mitosis_stoichiometry_olsen']

dir = os.getcwd()
#set a target cutoff defined in table
date = '20210830'

x_med ='cysteine_Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
x_prot_med ='protein_Mitosis_LPP(-,-)/(+,-)_combined_median'
y_prot_med = 'protein_Mitosis/Asynch_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
df.columns = df.columns.str.replace('Mitosis_LPP','cysteine_Mitosis_LPP')

lpp_df = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
phospho_df = pd.read_csv('{}/other_datasets/final_phospho_mann_olsen.csv'.format(dir))


print('Finished reading input data')
pdf =phospho_df[['accession','mitosis_stoichiometry_mann','mitosis_stoichiometry_olsen']]
pdf = pdf.replace('; ', ',', regex = True)




def test(row):
    mann = row['mann']
    olsen = row['olsen']
    if mann is np.nan:
        if olsen:
            p = olsen
            return p
        else:
            return np.nan
        
    elif olsen is np.nan:
        if mann:
            p = mann
            return mann
        else:
            return np.nan
    else:
        p = mann + olsen
        return p

pdf['mann'] = pdf.mitosis_stoichiometry_mann.str.split(',')
pdf['olsen'] = pdf.mitosis_stoichiometry_olsen.str.split(',')
pdf1 = pdf.assign(phospho = pdf.apply(test, axis = 1))
pdf1.drop(['mann', 'olsen', 'mitosis_stoichiometry_mann', 'mitosis_stoichiometry_olsen'], axis = 1, inplace = True)
pdf2 = pdf1.explode('phospho')
pdf2 = pdf2.replace('\)', '', regex = True)
pdf2['stoich'] = pdf2.phospho.str.split('(', expand = True)[1].astype(float)
pdf2['aa'] = pdf2.phospho.str[0]
pdf2['phospho'] = pdf2.phospho.str.split('(', expand = True)[0]
pdf2['res'] = pdf2.phospho.str[1:].astype(int)
pdf2['accession_res'] = pdf2['accession'] + '_' + pdf2.res.astype(str)

def fig2_heatmap(df, x_med, y_med):
    #source data is main dataset
    print('making heatmap')
    sns.set(style="white")
    df = df.loc[(~df[x_med].isnull()) &( ~df[y_med].isnull())]
    df2 = df[['gene_res',y_med,x_med]]
    print(df2)
    df2 = df2.dropna(axis = 0, how = 'any').set_index('gene_res')
    df3 = np.log2(df2)
    df3['temp'] = df3[x_med] + df3[y_med]
    df3 = df3.sort_values('temp', ascending = False)
    df3 = df3.reset_index().drop(['gene_res', 'temp'], axis = 1)
    df3.columns = [ 'Mitosis/Asynch','Mitosis LPP(-)/(+)']
    plt.figure(figsize = (1.8,5))
    # cmap = ['#990000','#ff3333','#ffcccc','white', '#e6e6ff', '#4d4dff', '#000099']
    # cmap = sns.light_palette((240, 100, 30), n_colors = 10, as_cmap=True, input="husl")
    cmap = sns.light_palette('#000099', as_cmap = True)
    # cmap =sns.light_palette((240, 100, 30), as_cmap=True, input="husl")
    # cmap = sns.diverging_palette(15, 265, sep=40, n=256, s = 100, l = 45)
    ax = sns.heatmap(df3, cmap = cmap, vmin = 0)
    ax.figure.axes[-1].yaxis.label.set_size(10)
    ax.figure.axes[-1].tick_params(labelsize=10)
    ax.set_ylabel("Peptides", fontsize = 10, family = 'ARIAL')
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 45, ha = 'right')
    plt.yticks(np.arange(0,11000,1000).tolist(),np.arange(0,11000,1000).tolist() , fontsize = 10, family = 'ARIAL')
    plt.title('Phosphoenrichment', family = 'ARIAL', fontsize = 10)
    # ax.vlines([0.97,2], ymin = 0, ymax = 13797, color = 'black', lw = 1)
    plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrichment_heatmap.pdf'.format(dir,date), transparent = True, bbox_inches='tight')
fig2_heatmap(lpp_df, x_med, y_med)

# def white(value):
#     if value > np.log2(2.5):
#         return 3
#     elif value <= np.log2(2.5) and value > np.log2(1.8):
#         return 2
#     elif value <= np.log2(1.8) and value > np.log2(1.4):
#         return 1
#     elif value <= np.log2(1.4) and value >= np.log2(1/1.4):
#         return 0
#     elif value < np.log2(1/1.4) and value >= np.log2(1/1.8):
#         return -1
#     elif value <=np.log2(1/1.8) and value >= np.log2(1/2.5):
#         return -2
#     elif value < np.log2(1/2.5):
#         return -3
    
# def fig2_prot_heatmap(df, x_med, y_med):
#     #source data is main dataset
#     print('making heatmap')
#     df1 = df.loc[(~df[x_med].isnull()) &( ~df[y_med].isnull())]
#     df2 = df1[['accession',y_med,x_med]]
#     print(df2)
#     df2 = df2.dropna(axis = 0, how = 'any').set_index('accession')
#     df3 = np.log2(df2)
#     df3 = df3.applymap(white)
    
#     df3['temp'] = df3[x_med] + df3[y_med]
#     df3 = df3.sort_values('temp', ascending = False)
#     df3 = df3.reset_index().drop(['accession', 'temp'], axis = 1)
#     df3.columns = [ 'Mitosis/Asynch', 'Mitosis LPP(-)/(+)']
#     plt.figure(figsize = (1.8,5))
#     cmap = ['#990000','#ff3333','#ffcccc','white', '#e6e6ff', '#4d4dff', '#000099']
#     # cmap =sns.color_palette("seismic_r", as_cmap=True)
#     # cmap = sns.diverging_palette(15, 265, sep=40, n=256, s = 100, l = 30)
#     ax = sns.heatmap(df3, cmap = cmap)
#     ax.figure.axes[-1].yaxis.label.set_size(10)
#     ax.figure.axes[-1].tick_params(labelsize=10)
#     ax.set_ylabel("Proteins", fontsize = 10, family = 'ARIAL')
#     plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 45, ha = 'right')
#     plt.yticks(np.arange(0,9000,500).tolist(),np.arange(0,9000,500).tolist() , fontsize = 10, family = 'ARIAL')
#     plt.title('Protein expression', family = 'ARIAL', fontsize = 10)
#     # ax.vlines([0.97,2], ymin = 0, ymax = 13797, color = 'black', lw = 1)
#     plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrichment_protein_heatmap.pdf'.format(dir,date), transparent = True, bbox_inches='tight')
# fig2_prot_heatmap(lpp_df, x_prot_med, y_prot_med)



# def supp_fig2de_bar(df, gene_res):
#     rep_cols = [col for col in df.columns if '_TMT' in col]
#     sdf = df.loc[df.gene_res == gene_res]
#     df1 = df.loc[df.gene_res == gene_res, rep_cols]
#     rename_rep_cols = [x.split('_', maxsplit = 1)[0] for x in rep_cols]
#     df1.rename(columns = dict(zip(rep_cols, rename_rep_cols)), inplace = True)
#     df2 = df1.melt().dropna(axis = 0)
#     gene = gene_res.split('_')[0]
#     res = 'C'+ gene_res.split('_')[1]
#     title = gene + ' ' + res
    
#     fig, axs = plt.subplots(figsize = (1.5,2.4)) #figsize=(2.8, 2.8)
#     sns.stripplot(x= 'variable', y = 'value', data = df2, jitter = True, dodge=True, ax = axs)
#     # sns.boxplot(x = 'labels', y = 'data',  data = df, boxprops=dict(alpha=.5), \
#     #             linewidth = 0.75, showfliers=False, ax = axs) 
#     sns.barplot(x = 'variable', y = 'value', errwidth = 1, ci = 'sd', capsize=0.2, alpha = 0.7, data = df2, ax = axs)
#     plt.axhline(y=2, color = 'black', linestyle = 'dotted')
#     plt.title(title.format(id), family = 'ARIAL', fontsize = 10)
#     plt.ylabel('Mitosis LPP(-)/LPP(+)\nMS1 ion intensity',family = 'ARIAL', fontsize = 10)
#     plt.xticks(fontsize = 10, family = 'ARIAL')
#     plt.yticks(fontsize = 10, family = 'ARIAL')
#     # plt.ylim([0, 7])
#     plt.xlabel('')
#     axs.spines['right'].set_visible(False)
#     axs.spines['top'].set_visible(False)
#     # plt.savefig('{}/figures/suppfig2/potential/{}_phospho_{}.pdf'.format(dir, date, gene_res), bbox_inches='tight',transparent = True)
# supp_fig2de_bar(lpp_df, 'FXR1_611')

# def final_fig2ghi_boxplots(df,protein_name, sub, width):
#     gene1 = df.loc[df.protein == protein_name].reset_index()
#     # gene1.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig2ghi_{}.csv'.format(dir, protein_name))
#     gene_name = gene1.iloc[0,:].gene_res.split('_')[0]
#     gene1 = gene1.assign(loc = gene1['sequence'].str.find('*', 0) - 1)
#     gene1 = gene1.assign(aa = gene1.apply(find_aa, axis = 1))
#     x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
#     y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
#     gene1 = gene1.loc[(~gene1[x_med].isnull())& (~gene1[y_med].isnull())]
#     rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] + [col for col in gene1.columns if y_name + '_median_TMT' in col]
#     rep_cols = sorted(rep_cols)
#     gene1 = gene1.assign(gene_res1 = gene1.aa + gene1.gene_res.str.split('_', n =1, expand = True)[1])
#     gene = gene1[['gene_res','gene_res1'] + rep_cols]
#     gene = gene.melt(id_vars = ['gene_res','gene_res1']).dropna(axis = 0)
#     gene = gene.assign(pos = gene.gene_res.str.split('_', n =1, expand = True)[1].astype(int))
#     gene = gene.assign(variable = gene.variable.str.rsplit(pat = '_', n = 1, expand = True)[0])
#     gene = gene.assign(value = np.log2(gene.value))
    
#     gene = gene.sort_values('variable', ascending = True)
#     gene = gene.sort_values('pos', ascending = True)
#     base_reps = list(set([i.rsplit('_', 1)[0] for i in rep_cols]))
#     order = sorted(base_reps, reverse= False)
#     print(order, gene)
#     fig, axs = plt.subplots(figsize=(width, 2.5))
#     sns.boxplot(x = 'gene_res1', y = 'value',  hue = 'variable', data = gene, \
#                 boxprops=dict(alpha=.6), showfliers=False, hue_order = order, linewidth = 0.75, palette = ['b','g'], ax = axs)
#     sns.stripplot(x= 'gene_res1', y = 'value', hue = 'variable', data = gene, \
#                   alpha = 0.7, linewidth = 0.5, edgecolor = 'gray', hue_order = order, jitter = True, dodge=True, palette = ['b','g'], ax = axs)
#     plt.ylim(-5,5)
#     plt.title(gene_name, fontsize = 10)
#     handles, labels = axs.get_legend_handles_labels() 
#     # axs.get_legend().remove()
#     l = plt.legend(handles[0:2], ['Mitosis/Asynch','LPP(-)/LPP(+)'], bbox_to_anchor=(1, 0.625), loc=2, fontsize = 10, frameon = False) # bbox_to_anchor=(0.63, 0.95),
#     plt.axhline(y=1, color = 'black', linestyle = 'dotted')
#     plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
#     plt.ylabel('Mitosis/Asynch or LPP(-)/LPP(+)\nlog2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
#     plt.xlabel('')
#     axs.spines['right'].set_visible(False)
#     axs.spines['top'].set_visible(False)
#     plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
#     plt.yticks(fontsize = 10, family = 'ARIAL')
#     plt.title(gene_name, family = 'ARIAL', fontsize = 10)
#     plt.savefig('{}/figures/suppfig2/potential/{}_phospho_{}_{}.pdf'.format(dir, date, protein_name, sub), bbox_inches='tight',transparent = True)
# final_fig2ghi_boxplots(lpp_df,'FXR2', 'phospho', 1.2) #flnb
# final_fig2ghi_boxplots(df, 'FXR2', 'cys', 3.3) #flnb


def color_prep(df, x_med, y_med):
    if df[x_med].item() <= 1/lpp_cutoff:
        if df[y_med].item() <= 1/lpp_cutoff:
            colors = ['red','red']
            return colors
        elif df[y_med].item() >= lpp_cutoff:
            colors = ['red', 'blue']
            return colors
        else:
            colors = ['red','black']
            return colors
    elif df[x_med].item() >= lpp_cutoff:
        if  df[y_med].item() >= lpp_cutoff:
            colors = ['blue','blue']
            return colors
        elif df[y_med].item() <= 1/lpp_cutoff:
            colors = ['blue','red']
            return colors
        else: 
            colors = ['blue','black']
            return colors
    elif (df[y_med].item() <= 1/lpp_cutoff):
        colors = ['black','red']
        return colors
    elif (df[y_med].item() >= lpp_cutoff):
        colors = ['black', 'blue']
        return colors
    else:
        colors = ['black','black']
        return colors
        
            
# def fig3ghk_boxplots(df, gene_cys):
#     gene1 = df.loc[df.gene_res == gene_cys]
#     # gene1.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig3ghk_box_{}.csv'.format(dir, gene_cys))
#     x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
#     y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
#     colors = color_prep(gene1, y_med, x_med)
#     print(colors)
#     rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] +\
#         [col for col in gene1.columns if y_name + '_median_TMT' in col] 
#     rep_cols = [col for col in rep_cols if 'protein_' not in col] 
#     location = (gene1['sequence'].str.find('*', 0) - 1).item()
#     aa = gene1.sequence.item()[location]
#     cys = aa + gene1.gene_res.str.split('_', n =1, expand = True)[1]
#     gene_cys = (gene_cys.split('_')[0] + ' ' + cys).item()
#     gene = gene1[['gene_res'] + rep_cols]    
#     gene.columns = gene.columns.str.replace('_median_TMT_rep*', '')
#     # gene.columns = gene.columns.str.split('_', n = 1, expand = True)[1]
#     # print(gene.columns)
#     gene = gene.melt(id_vars = 'gene_res').dropna(axis = 0)
#     gene = gene.assign(variable1 = gene.variable.str.rsplit(pat = '_', n = 1, expand = True)[0])
#     gene = gene.assign(value = np.log2(gene.value))
#     gene = gene.sort_values('variable1', ascending = True)
#     fig, axs = plt.subplots( figsize=(1, 2.4))
#     error_kw=dict(lw=1, capsize=3, capthick=1)
#     sns.boxplot(x = 'variable1', y = 'value',  data = gene, palette = colors, \
#                 order = [y_name, x_name], \
#                 boxprops=dict(alpha=.5), linewidth = 0.75, showfliers=False) 
#     sns.stripplot(x = 'variable1', y = 'value', data = gene, palette = colors, order =  [y_name, x_name])
#     fig = plt.gcf()
#     plt.ylim(-4.5,4.5)
#     plt.ylabel('log2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
#     plt.xlabel('')
#     axs.spines['right'].set_visible(False)
#     axs.spines['top'].set_visible(False)
#     plt.xticks(ticks = [0,1], labels = ['Mitosis/Asynch','LPP(-, -)/(+, -)'], rotation = 40, ha = 'right', fontsize = 10, family = 'ARIAL')
#     plt.yticks(fontsize = 10, family = 'ARIAL')
#     plt.title(gene_cys, family = 'ARIAL', fontsize = 10)
#     plt.axhline(y=1, color = 'black', linestyle = 'dotted')
#     plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
#     plt.savefig('{}/figures/suppfig2/potential/{}_phospho_{}.pdf'.format(dir,date, gene_cys), bbox_inches='tight',transparent = True)

# fig3ghk_boxplots(lpp_df, 'FXR2_450')
# fig3ghk_boxplots(df, 'FXR2_270')


def find_aa(row):
    location = row['loc']
    aa = row.sequence[location]
    return aa

def supp_fig2de_bar(df, gene_res, sub):
    gene1 = df.loc[df.gene_res == gene_res]
    # gene1.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig3ghk_box_{}.csv'.format(dir, gene_cys))
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    gene1 = gene1.assign(loc = gene1['sequence'].str.find('*', 0) - 1)
    gene1 = gene1.assign(aa = gene1.apply(find_aa, axis = 1))
    gene1 = gene1.loc[(~gene1[x_med].isnull())& (~gene1[y_med].isnull())]
    gene1['count'] = 1
    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] + [col for col in gene1.columns if y_name + '_median_TMT' in col]
    rep_cols = sorted(rep_cols)
    gene_res1 = gene1.gene_res.str.split('_', n =1, expand = True)[0] + ' ' + gene1.aa + gene1.gene_res.str.split('_', n =1, expand = True)[1]
    gene_res1 = gene_res1.item()
    gene = gene1[['gene_res', 'count'] + rep_cols]
    gene = gene.melt(id_vars = 'gene_res').dropna(axis = 0)
    gene = gene.assign(pos = gene.gene_res.str.split('_', n =1, expand = True)[1].astype(int))
    gene = gene.assign(variable = gene.variable.str.rsplit(pat = '_', n = 1, expand = True)[0])
#     gene = gene.assign(value = np.log2(gene.value))
    fig, axs = plt.subplots(figsize = (1.8,2.4)) #figsize=(2.8, 2.8)
    sns.stripplot(x= 'variable', y = 'value', data = gene, jitter = True, dodge=True, ax = axs, palette = ['black','black', 'red'])
    # sns.boxplot(x = 'labels', y = 'data',  data = df, boxprops=dict(alpha=.5), \
    #             linewidth = 0.75, showfliers=False, ax = axs) 
    sns.barplot(x = 'variable', y = 'value', errwidth = 1, ci = 'sd', capsize=0.2, alpha = 0.7, data = gene, ax = axs, palette = ['black','black', 'red'])
    plt.axhline(y=2, color = 'black', linestyle = 'dotted')
    plt.title(gene_res1, family = 'ARIAL', fontsize = 10)
    plt.ylabel('MS3 reporter ion intensity\nnormalized to Mitosis LPP(-)',family = 'ARIAL', fontsize = 10)
    plt.xticks(ticks = [0,1, 2], labels = ['Mitosis LPP(-)','Asynch LPP(-)','LPP(+)'], rotation = 40, ha = 'right', fontsize = 10, family = 'ARIAL')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    # plt.ylim([0, 20.5])
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.savefig('{}/figures/suppfig2/potential/{}_bar_phospho_{}.pdf'.format(dir, date, gene_res1, sub), bbox_inches='tight',transparent = True)
    
supp_fig2de_bar(df, 'FXR2_270', 'cys')
supp_fig2de_bar(lpp_df, 'FXR2_450', 'phospho')

# final_fig2ghi_boxplots(lpp_df,'FXR2', 'phospho', 1.2) #flnb
# final_fig2ghi_boxplots(df, 'FXR2', 'cys', 3.3) #flnb
