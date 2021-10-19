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



def fig2e_venn_prep():
    #source data generated here
   
    pdf = phospho_df[['accession', 'mitosis_stoichiometry_mann', 'mitosis_stoichiometry_olsen']]
    pdf.loc[pdf.mitosis_stoichiometry_olsen.notnull(), 'mitosis_stoichiometry_mann'] = pdf['mitosis_stoichiometry_mann'] + '; '
    pdf['phospho'] = pdf['mitosis_stoichiometry_mann'].fillna('') + pdf['mitosis_stoichiometry_olsen'].fillna('')
    pdf['phospho'] = pdf.phospho.str.split('; ')
    drop_cols = [col for col in pdf.columns if 'stoichiometry' in col]
    pdf.drop(drop_cols, axis = 1, inplace = True)
    pdf1 = pdf.explode('phospho')
    pdf1 = pdf1.replace('\)','', regex = True)
    pdf1[['res', 'stoichiometry']] = pdf1.phospho.str.split('(', expand = True)
    pdf1['aa'] = pdf1.res.str[0]
    pdf1['accession_res'] = pdf1.accession + '_' + pdf1.res.str[1:]
    pdf2 = pdf1.drop(['phospho', 'res'], axis = 1)
    pdf3 = pdf2.loc[~pdf2.duplicated(keep = 'first')].reset_index(drop = True)
    
    
    # lpdf1 = lpp_df[['accession','accession_res', x_med, y_med, 'cysteine_Mitosis_LPP(-,-)/(+,-)_combined_no',\
    #         'cysteine_Mitosis/Asynch_combined_no','cysteine_Mitosis_LPP(-,-)/(+,-)_combined_CV',\
    #             'cysteine_Mitosis/Asynch_combined_CV']]
    # lpdf1['accession_res1'] = lpdf1.accession_res.str.split('_', n = 1, expand = True)[1]
    # lpdf1['accession_res1'] = lpdf1.accession_res1.str.rsplit('_', n = 1, expand = True)[0]
    # lpdf1 = lpdf1.sort_values('accession_res', ascending = True).reset_index(drop = True)
    # lpdf2 = lpdf1.loc[lpdf1.accession_res1.duplicated(keep = 'first')]
    
    lpp_df['accession_res1'] = lpp_df.accession_res.str.rsplit('_', n = 1, expand = True)[0]
    
    df1 = pdf3.loc[pdf3.accession.isin(list(set(lpp_df.accession)))]
    df2 = pdf3.loc[pdf3.accession_res.isin(list(set(lpp_df.accession_res1)))]
    
    #df1 the number of proteins that we quantify that mann or olsen have identified as having high stoich
    len(set(df1.accession))
    len(set(phospho_df.accession))
    len(set(lpp_df.accession))

    only_x = list(set(set(lpp_df.accession_res1) - set(pdf3.accession_res)))
    only_y = list(set(pdf3.accession_res) - set(lpp_df.accession_res1))
    total = list(set(list(lpp_df.accession_res1) + list(pdf3.accession_res)))
    both = list( set(total) - set(only_x) - set(only_y))
    sizes = [ len(only_y), len(only_x),len(both)]
    print(sizes)
    
    plt.figure(figsize=(3.2, 3.2))
    out = venn2(subsets = sizes, set_labels = ('High stoichiometry sites\nin Mann et al. or Olsen et al.', 'Any quantified sites\nthis paper'), set_colors = ('g','b'), alpha = 0.9)
    for text in out.set_labels:
        text.set_fontsize(10)
    for text in out.subset_labels:
        text.set_fontsize(10)
        plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrich_sites_venn.pdf'.format(dir,date), bbox_inches= 'tight')
        
    sizes_sub = [ len(only_y),len(both)]
    total = sum(sizes_sub)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.4,2.4)
    ax1.pie(sizes_sub, autopct = (lambda p: '{:.0f}'.format(p * total / 100)), startangle=90)
    ax1.legend(labels = ['Only in Mann et al./Olsen et al.', 'Sites also quantified\nin this paper'], loc='center left', bbox_to_anchor=(0.02, -0.3),frameon=False, handletextpad=0.3, fontsize = 10)
    plt.title('High stoichiometry phosphorylation sites in Mann et al./Olsen et al', fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrich_sites_pie.pdf'.format(dir, date), bbox_inches='tight',transparent = True)
    
    only_x = list(set(set(lpp_df.accession) - set(phospho_df.accession)))
    only_y = list(set(phospho_df.accession) - set(lpp_df.accession))
    total = list(set(list(lpp_df.accession) + list(pdf3.accession)))
    both = list( set(total) - set(only_x) - set(only_y))
    sizes1 = [ len(only_y), len(only_x),len(both)]
    print(sizes1)
    
    plt.figure(figsize=(3.2, 3.2))
    out = venn2(subsets = sizes1, set_labels = ('Proteins w/ high stoichiometry sites\nin Mann et al. or Olsen et al.', 'Proteins w/ quantified\nphosphosites in this paper'), set_colors = ('g','b'), alpha = 0.9)
    for text in out.set_labels:
        text.set_fontsize(10)
    for text in out.subset_labels:
        text.set_fontsize(10)
        plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrich_protein_venn.pdf'.format(dir,date), bbox_inches= 'tight')
    
    sizes_sub1 = [ len(only_y),len(both)]
    total1 = sum(sizes_sub1)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.4,2.4)
    ax1.pie(sizes_sub1, autopct = (lambda p: '{:.0f}'.format(p * total1 / 100)), startangle=90)
    ax1.legend(labels = ['Only found in Mann et al./Olsen et al.', 'Proteins with quantified sites\nin this paper'], loc='center left', bbox_to_anchor=(0.02, -0.3),frameon=False, handletextpad=0.3, fontsize = 10)
    plt.title('Proteins w/ high stoichiometryphosphorylation\nsite in Mann et al/Olsen et als', fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrich_protein_pie.pdf'.format(dir, date), bbox_inches='tight',transparent = True)

fig2f_pie(prots_df)


print('Finished reading input data')


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
    plt.figure(figsize = (1.8,4))
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

def white(value):
    if value > np.log2(2.5):
        return 3
    elif value <= np.log2(2.5) and value > np.log2(1.8):
        return 2
    elif value <= np.log2(1.8) and value > np.log2(1.4):
        return 1
    elif value <= np.log2(1.4) and value >= np.log2(1/1.4):
        return 0
    elif value < np.log2(1/1.4) and value >= np.log2(1/1.8):
        return -1
    elif value <=np.log2(1/1.8) and value >= np.log2(1/2.5):
        return -2
    elif value < np.log2(1/2.5):
        return -3
    
def fig2_prot_heatmap(df, x_med, y_med):
    #source data is main dataset
    print('making heatmap')
    df1 = df.loc[(~df[x_med].isnull()) &( ~df[y_med].isnull())]
    df2 = df1[['accession',y_med,x_med]]
    print(df2)
    df2 = df2.dropna(axis = 0, how = 'any').set_index('accession')
    df3 = np.log2(df2)
    df3 = df3.applymap(white)
    
    df3['temp'] = df3[x_med] + df3[y_med]
    df3 = df3.sort_values('temp', ascending = False)
    df3 = df3.reset_index().drop(['accession', 'temp'], axis = 1)
    df3.columns = [ 'Mitosis/Asynch', 'Mitosis LPP(-)/(+)']
    plt.figure(figsize = (1.8,5))
    cmap = ['#990000','#ff3333','#ffcccc','white', '#e6e6ff', '#4d4dff', '#000099']
    # cmap =sns.color_palette("seismic_r", as_cmap=True)
    # cmap = sns.diverging_palette(15, 265, sep=40, n=256, s = 100, l = 30)
    ax = sns.heatmap(df3, cmap = cmap)
    ax.figure.axes[-1].yaxis.label.set_size(10)
    ax.figure.axes[-1].tick_params(labelsize=10)
    ax.set_ylabel("Proteins", fontsize = 10, family = 'ARIAL')
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 45, ha = 'right')
    plt.yticks(np.arange(0,9000,500).tolist(),np.arange(0,9000,500).tolist() , fontsize = 10, family = 'ARIAL')
    plt.title('Protein expression', family = 'ARIAL', fontsize = 10)
    # ax.vlines([0.97,2], ymin = 0, ymax = 13797, color = 'black', lw = 1)
    plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrichment_protein_heatmap.pdf'.format(dir,date), transparent = True, bbox_inches='tight')
fig2_prot_heatmap(lpp_df, x_prot_med, y_prot_med)



def find_aa(row):
    location = row['loc']
    aa = row.sequence[location]
    return aa
def final_fig2ghi_boxplots(df,protein_name, sub):
    if sub == 'cys':
        x_med1 = x_med
    elif sub == 'phospho':
        x_med1 = 'cysteine_' + x_med
    gene1 = df.loc[df.protein == protein_name].reset_index()
    # gene1.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig2ghi_{}.csv'.format(dir, protein_name))
    gene_name = gene1.iloc[0,:].gene_res.split('_')[0]
    gene1 = gene1.assign(loc = gene1['sequence'].str.find('*', 0) - 1)
    gene1 = gene1.assign(aa = gene1.apply(find_aa, axis = 1))
    x_name = x_med1.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    gene1 = gene1.loc[(~gene1[x_med1].isnull())& (~gene1[y_med].isnull())]
    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] + [col for col in gene1.columns if y_name + '_median_TMT' in col]
    rep_cols = sorted(rep_cols)
    gene1 = gene1.assign(gene_res1 = gene1.aa + gene1.gene_res.str.split('_', n =1, expand = True)[1])
    gene = gene1[['gene_res','gene_res1'] + rep_cols]
    gene = gene.melt(id_vars = ['gene_res','gene_res1']).dropna(axis = 0)
    gene = gene.assign(pos = gene.gene_res.str.split('_', n =1, expand = True)[1].astype(int))
    gene = gene.assign(variable = gene.variable.str.rsplit(pat = '_', n = 1, expand = True)[0])
    gene = gene.assign(value = np.log2(gene.value))
    
    gene = gene.sort_values('variable', ascending = True)
    gene = gene.sort_values('pos', ascending = True)
    base_reps = list(set([i.rsplit('_', 1)[0] for i in rep_cols]))
    order = sorted(base_reps, reverse= False)
    print(order, gene)
    width = 0.5 * len(set(gene.gene_res1))
    fig, axs = plt.subplots(figsize=(width, 2.5))
    sns.boxplot(x = 'gene_res1', y = 'value',  hue = 'variable', data = gene, \
                boxprops=dict(alpha=.6), showfliers=False, hue_order = order, linewidth = 0.75, palette = ['b','g'], ax = axs)
    sns.stripplot(x= 'gene_res1', y = 'value', hue = 'variable', data = gene, \
                  alpha = 0.7, linewidth = 0.5, edgecolor = 'gray', hue_order = order, jitter = True, dodge=True, palette = ['b','g'], ax = axs)
    plt.ylim(-5,5)
    plt.title(gene_name, fontsize = 10)
    handles, labels = axs.get_legend_handles_labels() 
    # axs.get_legend().remove()
    l = plt.legend(handles[0:2], ['Mitosis/Asynch','LPP(-)/LPP(+)'], bbox_to_anchor=(1, 0.625), loc=2, fontsize = 10, frameon = False) # bbox_to_anchor=(0.63, 0.95),
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.ylabel('Mitosis/Asynch or LPP(-)/LPP(+)\nlog2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.title(gene_name, family = 'ARIAL', fontsize = 10)
    plt.savefig('{}/figures/suppfig2/potential/{}_phospho_{}_{}.pdf'.format(dir, date, protein_name, sub), bbox_inches='tight',transparent = True)
    
final_fig2ghi_boxplots(lpp_df,'BAG', 'phospho') #flnb
final_fig2ghi_boxplots(lpp_df,'FLNB', 'phospho', 1.2) #flnb
final_fig2ghi_boxplots(lpp_df,'FLNB', 'phospho', 1.2) #flnb
final_fig2ghi_boxplots(df, 'FXR2', 'cys', 3.3) #flnb



            
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
#     # plt.close()
    
# fig3ghk_boxplots(lpp_df, 'FXR2_450')
# fig3ghk_boxplots(df, 'FXR2_270')

def color_prep(df, x_med, y_med):
    if df[x_med].item() <= 1/lpp_cutoff:
        if df[y_med].item() <= 1/lpp_cutoff:
            colors = ['black','red','red']
            return colors
        elif df[y_med].item() >= lpp_cutoff:
            colors = ['black','red', 'blue']
            return colors
        else:
            colors = ['black','red','black']
            return colors
    elif df[x_med].item() >= lpp_cutoff:
        if  df[y_med].item() >= lpp_cutoff:
            colors = ['black','blue','blue']
            return colors
        elif df[y_med].item() <= 1/lpp_cutoff:
            colors = ['black','blue','red']
            return colors
        else: 
            colors = ['black','blue','black']
            return colors
    elif (df[y_med].item() <= 1/lpp_cutoff):
        colors = ['black','black','red']
        return colors
    elif (df[y_med].item() >= lpp_cutoff):
        colors = ['black','black', 'blue']
        return colors
    else:
        colors = ['black','black','black']
        return colors
        

def supp_fig2de_bar(df, gene_res,sub):
    gene1 = df.loc[df.gene_res == gene_res]
    # gene1.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig3ghk_box_{}.csv'.format(dir, gene_res))
    colors = color_prep(gene1, y_med, x_med)
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] +\
        [col for col in gene1.columns if y_name + '_median_TMT' in col] 
    rep_cols = [col for col in rep_cols if 'protein_' not in col] 
    location = (gene1['sequence'].str.find('*', 0) - 1).item()
    aa = gene1.sequence.item()[location]
    res = aa + gene1.gene_res.str.split('_', n =1, expand = True)[1]
    gene_res1 = (gene_res.split('_')[0] + ' ' + res).item()
    gene1['count'] = 1
    gene = gene1[['gene_res', 'count'] + rep_cols]    
    gene.columns = gene.columns.str.replace('_median_TMT_rep*', '')
    # gene.columns = gene.columns.str.split('_', n = 1, expand = True)[1]
    # print(gene.columns)
    gene = gene.melt(id_vars = 'gene_res').dropna(axis = 0)
    gene = gene.assign(variable1 = gene.variable.str.rsplit(pat = '_', n = 1, expand = True)[0])
    gene = gene.sort_values('variable1', ascending = True)
    gene['value1'] = 1/gene['value']

    fig, axs = plt.subplots(figsize = (1.5,2.4)) #figsize=(2.8, 2.8)
    sns.stripplot(x= 'variable1', y = 'value1', data = gene, jitter = True, dodge=True, ax = axs, palette = colors)
    # sns.boxplot(x = 'labels', y = 'data',  data = df, boxprops=dict(alpha=.5), \
    #             linewidth = 0.75, showfliers=False, ax = axs) 
    sns.barplot(x = 'variable1', y = 'value1', errwidth = 1, ci = 'sd', capsize=0.2, alpha = 0.7, data = gene, ax = axs, palette = colors)
    plt.axhline(y=0.5, color = 'black', linestyle = 'dotted')
    plt.ylabel('% of LPP(-)\nMS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xticks(ticks = [0,1,2], labels = ['Mitosis LPP(-)','Asynch LPP(-)', 'Mitosis LPP/(+)'], rotation = 40, ha = 'right', fontsize = 10, family = 'ARIAL')
    plt.yticks(ticks = np.arange(0,1.2,0.2), fontsize = 10, family = 'ARIAL')
    plt.title(gene_res1, family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.savefig('{}/figures/suppfig2/potential/{}_phospho_bar_{}.pdf'.format(dir, date, gene_res1, sub), bbox_inches='tight',transparent = True)
supp_fig2de_bar(lpp_df, 'FXR2_450', 'phospho')
supp_fig2de_bar(df, 'FXR2_270', 'cys')
supp_fig2de_bar(lpp_df, 'SLAIN2_147', 'phospho')


df1 = df.loc[(df[x_med] <= 1/lpp_cutoff) | (df[x_med] >= lpp_cutoff)]
df2 = df1.loc[(df[y_med] <= 1/asynch_cutoff) | (df1[y_med] >= asynch_cutoff)]

down = df2.loc[(df2[x_med] <= 1/lpp_cutoff)]
up = df2.loc[(df2[x_med] >= lpp_cutoff)]
authentic = down.loc[down[adapted_med] <= 1/adapted_cutoff]
artifact = down.loc[down[adapted_med] > 1/unchanging_cutoff]

all_changing = df2.loc[(df2[x_med] <= 1/lpp_cutoff) | (df2[x_med] >= lpp_cutoff)]
no_ppep = all_changing.loc[~((all_changing[x_med] <= 1/lpp_cutoff) & (all_changing[adapted_med] > 1/unchanging_cutoff))]

lpp_sens = lpp_df.loc[lpp_df[x_med] >= 2]
p_a = authentic.loc[authentic.accession.isin(list(set(lpp_sens.accession)))]
t_authentic = len(set(authentic.accession))
p_authentic = len(set(p_a.accession))

p_all_df = df.loc[df.accession.isin(list(set(lpp_sens.accession)))]
t_all = len(set(df.accession))
p_all = len(set(p_all_df.accession))

p_changing_df = all_changing.loc[all_changing.accession.isin(list(set(lpp_sens.accession)))]
t_changing = len(set(all_changing.accession))
p_changing = len(set(p_changing_df.accession))

p_artifact_df = artifact.loc[artifact.accession.isin(list(set(lpp_sens.accession)))]
t_artifact = len(set(artifact.accession))
p_artifact = len(set(p_artifact_df.accession))

p_no_ppep_df = no_ppep.loc[no_ppep.accession.isin(list(set(lpp_sens.accession)))]
t_no_ppep = len(set(no_ppep.accession))
p_no_ppep = len(set(p_no_ppep_df.accession))

p_up_df = up.loc[up.accession.isin(list(set(lpp_sens.accession)))]
t_up = len(set(up.accession))
p_up = len(set(p_up_df.accession))


def fig4a_pie(titles, sizes, sub):
    labels = 'Proteins with LPP-sensitive\nphosphorylation sites','Proteins with phosphosites\nnot quantified or unchanging'    
    total = sum(sizes)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.1,2)
    plt.title(titles, fontsize = 10)
    ax1.pie(sizes, autopct = (lambda p: '{:.0f}'.format(p * total / 100)) , startangle=90)
    # ax1.pie(sizes, autopct = (lambda p: '{:.0f}/476'.format(p * total / 100)) , startangle=90)
    ax1.legend(labels = labels, loc='center left', bbox_to_anchor=(-0.2, -0.15),frameon=False, handletextpad=0.3, fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/suppfig2/potential/{}_pie_{}.pdf'.format(dir, date, sub), bbox_inches='tight',transparent = True)
    plt.close()

def fig4a_generation():
    sizes = [p_all, t_all - p_all]
    titles = 'All quantified proteins'
    sub = 'all_proteins'
    fig4a_pie(titles, sizes, sub)
    
    titles = 'Proteins with authentic phosphorylation-\ndependent cysteine reactivity changes'
    sizes = [p_authentic, t_authentic - p_authentic]
    sub = 'authentic'
    fig4a_pie(titles, sizes, sub)
    
    titles = 'All proteins with phosphorylation-\ndependent cysteine reactivity changes'
    sizes = [p_changing, t_changing - p_changing]
    sub = 'all_changing'
    fig4a_pie(titles, sizes, sub)
    
    titles = 'Proteins with artifactual phosphorylation-\ndependent cysteine reactivity changes'
    sizes = [p_artifact, t_artifact - p_artifact]
    sub = 'artifact'
    fig4a_pie(titles, sizes, sub)
    
    titles = 'Phosphorylation-dependent cysteine reactivity\nchanges WITHOUT artifactual changes'
    sizes = [p_no_ppep, t_no_ppep - p_no_ppep]
    sub = 'no_ppep'
    fig4a_pie(titles, sizes, sub)
    
    titles = 'Phosphorylation-dependent cysteine reactivity\nchanges up in LPP(+)'
    sizes = [p_up, t_up - p_up]
    sub = 'up'
    fig4a_pie(titles, sizes, sub)
fig4a_generation()

def test(row):
    start = row['start']
    end = row['end']
    
    
lpp_sens = lpp_df.loc[lpp_df[x_med] >= 2]
artifact = down.loc[down[adapted_med] > 1/unchanging_cutoff]

# pos = lpp_sens.sequence.str.find('*') - 1 
# length = lpp_sens.sequence.str.len() -1
# res = lpp_sens.gene_res.str.split('_', expand = True)[1].astype(int)
# lpp_sens['start'] = res - pos
# lpp_sens['end'] = (length -pos) + res - 1
# lpp1 = lpp_sens[['accession','gene_res','sequence','start','end',x_med,y_med]]
# art1 = artifact[['accession','protein','accession_res','gene_res','sequence',x_med,y_med]]
# art1['res'] = art1.gene_res.str.split('_', expand = True)[1].astype(int)
# art_lpp = art1.merge(lpp1, on = ['accession'], how = 'outer', suffixes = ['_CYS', '_PHOSPHO'])

def on_tryptic(lpp_sens,artifact):
    lpp_sens['res'] = lpp_sens.gene_res.str.split('_', expand = True)[1].astype(int)
    lpp1 = lpp_sens[['accession','sequence','res',x_med,y_med]]
    art1 = artifact[['accession','protein','accession_res','gene_res','sequence',x_med,y_med]]
    pos = art1.sequence.str.find('*') - 1 
    length = art1.sequence.str.len() -1
    res = art1.gene_res.str.split('_', expand = True)[1].astype(int)
    art1['start'] = res - pos
    art1['end'] = (length -pos) + res - 1
    art_lpp = art1.merge(lpp1, on = ['accession'], how = 'outer', suffixes = ['_CYS', '_PHOSPHO'])
    on_pep = art_lpp.loc[(art_lpp.res >= art_lpp.start) & (art_lpp.res <= art_lpp.end)]
    return on_pep
on_pep = on_tryptic(lpp_sens,artifact)
on_pep1 = on_tryptic(lpp_df,artifact)
on_pep2 = on_tryptic(lpp_sens,authentic)
on_pep3 = on_tryptic(lpp_sens,df)
len(set(on_pep.accession_res))
len(set(on_pep1.accession_res))
len(set(on_pep2.accession_res))
len(set(on_pep3.accession_res))
len(set(art1.accession_res))
len(set(authentic.accession_res))
len(set(df.accession_res))

a = len(set(on_pep.accession_res))
b = len(set(on_pep1.accession_res) - set(on_pep.accession_res))
c =len(set(art1.accession_res) - set(on_pep1.accession_res))

on_pep.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_artifactual_same_pep.csv'.format(dir), index = False)

def fig4a_pie():
    sizes = [a,b,c]
    labels = 'LPP-sensitive phosphosite\non same tryptic','LPP-insensitve phosphosite\non same tryptic','No quantified phosphosite\non same tryptic'    
    total = sum(sizes)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.1,2)
    plt.title('Cysteines with artifactual phoshorylation-\ndependent cysteine reactivity changes', fontsize = 10)
    ax1.pie(sizes, autopct = (lambda p: '{:.0f}'.format(p * total / 100)) , startangle=90)
    # ax1.pie(sizes, autopct = (lambda p: '{:.0f}/476'.format(p * total / 100)) , startangle=90)
    ax1.legend(labels = labels, loc='center left', bbox_to_anchor=(-0.2, -0.3),frameon=False, handletextpad=0.3, fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/suppfig2/potential/{}_pie_same_tryptic_artifact.pdf'.format(dir, date), bbox_inches='tight',transparent = True)
    plt.close()