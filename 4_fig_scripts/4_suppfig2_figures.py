#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 13:51:09 2021

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
import matplotlib.ticker as ticker

target_cutoff = 2
lpp_cutoff = 2
asynch_cutoff = 1.6
unchanging_cutoff = 1.5
protein_cutoff = 1.6

protein_cols = ['accession','accession_res','description', 'protein', 'gene_res']
phospho_col = ['pos_seq_mann', 'asynch_stoichiometry_mann','mitosis_stoichiometry_mann', 'stoichiometry','pos_seq_olsen', 'asynch_stoichiometry_olsen','mitosis_stoichiometry_olsen']

dir = os.getcwd()
#set a target cutoff defined in table
date = '20210830'

x_med ='Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
adapted_med ='Mitosis_LPP(-,+)/(+,+)_combined_median'
denature_med = 'Mitosis_Native/Denatured_combined_median'

# prot_med = 'protein_Mitosis/Asynch_ratio_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
y_df = pd.read_excel('{}/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date), sheet_name = 'filtered')
y_exp = pd.read_excel('{}/supplemental_tables/{}_fig1_proteins.xlsx'.format(dir,date), sheet_name = 'all prots with changes')
desalt_df = pd.read_excel('{}/supplemental_tables/{}_isotop_desalting_combined.xlsx'.format(dir,date))
denature_df = pd.read_excel('{}/supplemental_tables/{}_fig2_denature_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
dis_map = pd.read_csv('{}/other_datasets/IUPRED/20210303_all_fasta.csv'.format(dir))

print('Finished reading input data')

def supp_fig2b_waterfall(df, x_med, label):
    #source data main dataset
    plt.subplots(figsize=(3,2.7))
    df1 = df[['gene_res', 'accession', x_med]] 
    df1 = df1.loc[~(df1[x_med] == 'NQ')]
    df1.loc[:,x_med] = np.log2(df1.loc[:,x_med].astype('float'))
    df1 = df1.sort_values(x_med, ascending = True)
    y = df1[[x_med]]
    x = np.arange(0,len(df1))
    up_cond = (df1[x_med] > np.log2(lpp_cutoff))
    down_cond = (df1[x_med] <= np.log2(1/lpp_cutoff))
    both = up_cond| down_cond
    plt.scatter(x[~both], y[~both], lw = 0, c = '#404040', marker = 'o',s = 8, alpha = 0.6)
    plt.scatter(x[up_cond], y[up_cond], lw = 0, c = '#00b7eb', marker = 'o', s = 8, alpha = 0.6)
    plt.scatter(x[down_cond], y[down_cond], lw = 0, c = '#F28C28', marker = 'o', s = 8, alpha = 0.6)
    plt.axhline(y = 1, c = 'black', linestyle = 'dotted')
    plt.axhline(y = -1, c = 'black', linestyle = 'dotted')
    plt.xlabel('Peptides', fontsize =10, family = 'ARIAL')
    plt.ylabel(label + '\nlog2 MS3 reporter ion intensity',fontsize =10, family = 'ARIAL')
    plt.xticks(fontsize =10, family = 'ARIAL')
    plt.yticks(fontsize =10, family = 'ARIAL', ticks= [-4,-2,0,2,4])
    title = x_med.split('/')[0]
    plt.savefig('{}/figures/suppfig2/{}_suppfig2b_waterfall.pdf'.format(dir,date), bbox_inches='tight')
    plt.close()
label = 'Native/Denatured'
supp_fig2b_waterfall(denature_df, denature_med, label)


def supp_fig2c_prep():
    # dis_map = fig4d_disorder_map()
    # lig_map = fig4d_lig_map()
    dis_map = pd.read_csv('{}/other_datasets/IUPRED/20210303_all_fasta.csv'.format(dir))
    dis_map = dis_map.assign(accession_res = (dis_map.accession + '_C' + dis_map.accession_res.str.rsplit('_',expand = True)[1]))
    dis_map = dis_map[['accession_res', 'IUPRED SCORE']]
    df1 = denature_df.loc[~denature_df[denature_med].isnull()]
    # lpp_cond = (df1[denature_med] >= lpp_cutoff) | (df1[denature_med] <= 1/lpp_cutoff)
    # am_diff = df1.loc[(lpp_cond)]
    drop_cols = [col for col in df1.columns if 'median_TMT' in col]#+ [col for col in df1.columns if '_CV' in col]+ [col for col in df1.columns if '_no' in col]
    df1 = df1.drop(drop_cols, axis = 1)
    df1 = df1.assign(accession_res = (df1.accession + '_C' + df1.gene_res.str.split('_',expand = True)[1]))
    # am_diff = am_diff.drop(drop_cols, axis = 1)
    # am_diff = am_diff.assign(accession_res = (am_diff.accession + '_C' + am_diff.gene_res.str.split('_',expand = True)[1])) 
    # am_diff1 = am_diff.merge(dis_map, how = 'left', on = 'accession_res')
    df2 = df1.merge(dis_map, how = 'left', on = 'accession_res')
    
    down =  df2.loc[(df2[denature_med] <= 1/lpp_cutoff)]
    up = df2.loc[(df2[denature_med] >= lpp_cutoff)]
    keep_cols = [col for col in df2.columns if 'combined_median' in col]
    sdf = df2[protein_cols +['sequence']+ keep_cols + ['IUPRED SCORE']]
    sdf[keep_cols] = np.log2(sdf[keep_cols])
    sdf['Reactivity in Native proteome'] =  'Unchanged'
    sdf.loc[sdf.accession_res.isin(list(set(down.accession_res))), 'Reactivity in Native proteome'] = 'Decreased'
    sdf.loc[sdf.accession_res.isin(list(set(up.accession_res))), 'Reactivity in Native proteome'] = 'Increased'
    sdf['Reside in disordered region IUPRED > 0.5?'] = False
    sdf.loc[(sdf['IUPRED SCORE']> 0.5), 'Reside in disordered region IUPRED > 0.5?'] = True
    sdf.to_excel('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig2c.xlsx'.format(dir), index = False, float_format = '%.2f')
    df3 = df2[~df2.accession_res.isin(set(down.accession_res))]
    df3 = df3[~df3.accession_res.isin(set(up.accession_res))]
    # df3.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_supp_fig2c.csv'.format(dir), index = False, float_format = '%.2f')
    return up, down,  df2, df3

def supp_fig2c_disorder_prep():
    up, down, df2, df3 = supp_fig2c_prep()
    labels = ['Lower','Unchanging','Higher']
    dis_number = [len(set(down.loc[down['IUPRED SCORE'] >= 0.5, 'accession_res'])),\
        len(set(df3.loc[df3['IUPRED SCORE'] >= 0.5, 'accession_res'])),\
                 len(set(up.loc[up['IUPRED SCORE'] >= 0.5, 'accession_res']))\
                ]
    total_number = [len(set(down.accession_res)),len(set(df3.accession_res)),len(set(up.accession_res))]
    print(labels, dis_number, total_number)
    return labels, dis_number, total_number

def fig4c_plot_disorder():   
    labels, dis_number, total_number = supp_fig2c_disorder_prep()
    y_value =  [i / j for i, j in zip(dis_number, total_number)] 
    percent_labels = [str(i)+ '/\n' + str(j) for i,j in zip(dis_number, total_number)]    
    plt.subplots(figsize=(2.3,2.6))
    ax = sns.barplot(x = labels, y = y_value, palette = ('#F28C28','#404040', '#00b7eb'))
    plt.ylabel('Predicted instrinsic disorder\n(Percent of total)', fontsize =10, family = 'ARIAL')
    plt.xticks(fontsize =10, family = 'ARIAL', rotation = 45, ha = 'right')
    plt.yticks(fontsize =10, family = 'ARIAL')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.ylim(0,0.25)
    plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(0.05))
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1, decimals = 0))
    for i in range(len(percent_labels)):
        plt.annotate(str(percent_labels[i]), xy=(i,y_value[i]), ha='center', va='bottom', fontsize = 10, family = 'ARIAL')
    plt.savefig('{}/figures/suppfig2/{}_suppfig2c_disorder.pdf'.format(dir, date), bbox_inches='tight')
    plt.close()
fig4c_plot_disorder()

def supp_fig2d_correlation_plot(df, x_med, y_med):
    #source data in main dataset
    print('Running Supp Fig2 correlation plot')
    plt.subplots(figsize=(2.8,2.8)) 
    df1= df.loc[~(df[y_med] == 'NQ')]
    df1 = df1.loc[~(df1[x_med] == 'NQ')]
    df1 = df1[['gene_res', 'accession', x_med,y_med]]
    df1[[x_med,y_med]] = np.log2(df1[[x_med,y_med]].astype('float') )
    x = df1[y_med]
    y = df1[x_med]
    cond1 = (df1[x_med] >= np.log2(lpp_cutoff))  
    cond2 = (df1[x_med] <= np.log2(1/lpp_cutoff))
    cond3 =(df1[y_med] >= np.log2(lpp_cutoff)) & (~cond1) & (~cond2)
    cond4 = (df1[y_med] <= np.log2(1/lpp_cutoff))& (~cond1) & (~cond2)
    all = (cond1) | (cond2) | (cond3) | (cond4)
    plt.scatter(x[cond3], y[cond3], lw = 0.1, edgecolors = 'black', alpha = 0.5,c = '#9bddff', marker = '.')
    plt.scatter(x[cond4], y[cond4], lw = 0.1, edgecolors = 'black',alpha = 0.5,c = '#FFAC1C', marker = '.')
    plt.scatter(x[cond1], y[cond1], lw = 0.1, edgecolors = 'black',alpha = 0.5,c = 'blue', marker = '.')
    plt.scatter(x[cond2], y[cond2], lw = 0.1, edgecolors = 'black',alpha = 0.5,c = 'red', marker = '.')
    plt.scatter(x[~all], y[~all], lw = 0.1, edgecolors = 'black',alpha= 0.3, c = '#404040', marker = '.')
    plt.xlabel('Mitosis Native/Denatured\nlog2 MS3 reporter ion intensity', fontsize =10, family = 'ARIAL')
    plt.ylabel('Mitosis LPP(-)/LPP(+)\nlog2 MS3 reporter ion intensity',fontsize =10, family = 'ARIAL')
    plt.axhline(y = 1, c = 'black', linestyle = 'dotted')
    plt.axhline(y = -1, c = 'black', linestyle = 'dotted')
    plt.axvline(x = 1, c = 'black', linestyle = 'dotted')
    plt.axvline(x = -1, c = 'black', linestyle = 'dotted')
    plt.xlim(-4.5,4.5) 
    plt.ylim(-4.5,4.5) 
    plt.xticks(fontsize =10, family = 'ARIAL',ticks= [-4,-2,0,2,4])
    plt.yticks(fontsize =10, family = 'ARIAL',ticks= [-4,-2,0,2,4])
    plt.savefig('{}/figures/suppfig2/{}_suppfig2d_correlation.pdf'.format(dir,date), bbox_inches='tight')
    plt.close()
supp_fig2d_correlation_plot(denature_df, x_med, denature_med)

def supp_fig2efg_boxplots(df, x_med, y_med, protein):
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    df[[x_med, y_med]] = df[[x_med, y_med]].replace('NQ', np.nan)
    gene1 = df.loc[df.protein == protein]
    keep_cols = [col for col in gene1.columns if x_name in col] + [col for col in gene1.columns if y_name in col]
    median_cols = [col for col in gene1.columns if 'median' in col]
    sdf = gene1.copy()
    sdf[median_cols] = np.log2(sdf[median_cols])
    sdf = sdf.sort_values('protein', ascending = True)
    sdf.set_index(['accession', 'description', 'protein','sequence','accession_res','gene_res'], inplace = True)
    sdf = sdf[keep_cols]
    sdf.to_excel('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig2efg_{}.xlsx'.format(dir, protein))
    gene_name = gene1.iloc[0,:].gene_res.split('_')[0]
    gene1 = gene1.loc[~((gene1[x_med].isnull()) | (gene1[y_med].isnull()))]
    if len(set(gene1.gene_res)) > 12:
        x_no = x_name + '_combined_no'
        y_no = y_name + '_combined_no'
        gene1 = gene1.loc[(gene1[x_no]>=8) & (gene1[y_no] >= 2)]
    # elif len(set(gene1.gene_res)) > 9:
    #     x_no = x_name + '_combined_no'
    #     y_no = y_name + '_combined_no'
    #     gene1 = gene1.loc[(gene1[x_no]>=7) & (gene1[y_no] >= 2)]
    else:
        pass

    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] + [col for col in gene1.columns if y_name + '_median_TMT' in col]
    rep_cols = sorted(rep_cols)
    gene = gene1[['gene_res'] + rep_cols]
    gene = gene.melt(id_vars = 'gene_res').dropna(axis = 0)
    gene = gene.assign(pos = gene.gene_res.str.split('_', n =1, expand = True)[1].astype(int))
    gene = gene.assign(variable = gene.variable.str.rsplit(pat = '_', n = 1, expand = True)[0])
    gene = gene.assign(value = np.log2(gene.value))
    gene = gene.assign(gene_res = 'C' + gene.gene_res.str.split('_', n =1, expand = True)[1])
    gene = gene.sort_values('variable', ascending = True)
    gene = gene.sort_values('pos', ascending = True)
    base_reps = list(set([i.rsplit('_', 1)[0] for i in rep_cols]))
    order = sorted(base_reps, reverse= True)
    width = len(set(gene.pos)) * 0.45
    fig, axs = plt.subplots(figsize = (width, 2.6)) #figsize=(2.8, 2.8)
    sns.boxplot(x = 'gene_res', y = 'value',  hue = 'variable', data = gene, \
                boxprops=dict(alpha=.7), showfliers=False, hue_order = order, linewidth = 0.75, palette = ['#F28C28','g'], ax = axs)
    sns.stripplot(x= 'gene_res', y = 'value',hue = 'variable',\
                  linewidth = 0.5, edgecolor = 'gray', hue_order = order, data = gene, \
                  jitter = 0.2, dodge=True, palette = ['#F28C28','g'], ax = axs)
    # plt.ylim(-4.5,4.5)
    plt.title(gene_name, fontsize = 10)
    handles, labels = axs.get_legend_handles_labels() 
    axs.get_legend().remove()
    # l = plt.legend(handles[0:2], ['Native/Denature','LPP(-)/LPP(+)'], bbox_to_anchor=(1, 0.625), loc=2, fontsize = 10, frameon = False) # bbox_to_anchor=(0.63, 0.95),
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.ylim(-4.6,4.6)
    plt.ylabel('Native/Denature or LPP(-)/LPP(+)\nlog2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.title(gene_name, family = 'ARIAL', fontsize = 10)
    plt.savefig('{}/figures/suppfig2/{}_suppfig2efg_denature_bar_{}.pdf'.format(dir, date, gene_name), bbox_inches='tight',transparent = True)
    plt.close()
supp_fig2efg_boxplots(denature_df, x_med, denature_med, 'FLNB')
supp_fig2efg_boxplots(denature_df, x_med, denature_med, 'NUMA1')
supp_fig2efg_boxplots(denature_df, x_med, denature_med, 'BAG3')

def supp_fig2h_protein(protein):
    df1 = desalt_df.loc[desalt_df.protein == protein]
    desalt_no =  'Mitosis_Gel-filtered/Unfiltered_combined_no'
    df1 = df1.loc[df1[desalt_no] >= 2]
    df1.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig2h_protein_{}.csv'.format(dir, protein), float_format = '%.2f')
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
    fig, axs = plt.subplots(figsize = (width,2.6)) #figsize=(2.8, 2.8)
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
    plt.savefig('{}/figures/suppfig2/{}_suppfig2h_desalt_protein_{}.pdf'.format(dir, date, protein), bbox_inches='tight',transparent = True) 
    plt.close()
supp_fig2h_protein('MAP2K4')

def supp_fig2i_prep(gene_res):
    rep_cols = [col for col in desalt_df.columns if 'rep' in col]
    sdf = desalt_df.loc[desalt_df.gene_res == gene_res]
    sdf.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig2i_{}.csv'.format(dir, gene_res), float_format = '%.2f')
    df1 = desalt_df.loc[desalt_df.gene_res == gene_res, rep_cols]
    rename_rep_cols = [x.split('_', maxsplit = 1)[0] for x in rep_cols]
    df1.rename(columns = dict(zip(rep_cols, rename_rep_cols)), inplace = True)
    df2 = df1.melt().dropna(axis = 0)
    gene = gene_res.split('_')[0]
    res = 'C'+ gene_res.split('_')[1]
    title = gene + ' ' + res
    return df2, title

def supp_fig2i_bar(gene_res):
    df, title = supp_fig2i_prep(gene_res)
    fig, axs = plt.subplots(figsize = (1.1,2.6)) #figsize=(2.8, 2.8)
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
    plt.savefig('{}/figures/suppfig2/{}_suppfig2i_desalt_cys_{}_bar.pdf'.format(dir, date, gene_res), bbox_inches='tight',transparent = True)
    plt.close()
supp_fig2i_bar('MAP2K4_246')

def supp_fig2j_boxplots(cys_list, width):
    gene1 = df.loc[df.gene_res.isin(cys_list)]
    sdf = gene1.copy()
    adapted_name = adapted_med.rsplit(sep = '_', maxsplit = 2)[0]
    drop_cols = [col for col in sdf.columns if 'protein_' in col] + [col for col in sdf.columns if adapted_name in col] + [col for col in sdf.columns if 'limited' in col]
    drop_cols = drop_cols + [col for col in sdf.columns if 'category' in col]
    sdf.drop(drop_cols, axis =1, inplace = True)
    medians = [col for col in sdf.columns if 'median' in col]
    sdf[medians] = np.log2(sdf[medians])
    sdf.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig2j_MAP2K.csv'.format(dir))
    gene_name = gene1.iloc[0,:].gene_res.split('_')[0]
    gene_name = gene_name[0:-1]
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] + [col for col in gene1.columns if y_name + '_median_TMT' in col]
    rep_cols = [col for col in rep_cols if 'protein_' not in col]
    gene = gene1[['gene_res'] + rep_cols]
    gene = gene.melt(id_vars = 'gene_res').dropna(axis = 0)
    gene = gene.assign(pos = gene.gene_res.str.split('_', n =1, expand = True)[1].astype(int))
    gene = gene.assign(value = np.log2(gene.value))
    gene = gene.sort_values('gene_res', ascending = True)
    gene = gene.assign(variable = gene.variable.str.rsplit(pat = '_', n = 1, expand = True)[0])
    map2k7 = pd.DataFrame({'gene_res':['MAP2K7_260','MAP2K7_260'], \
                           'variable':['Mitosis_LPP(-,-)/(+,-)_median_TMT_rep','cysteine_Mitosis/Asynch_median_TMT_rep'],\
                               'value': [np.nan, np.nan], \
                                   'pos':[260,260]})
    gene = gene.append(map2k7)
    order = sorted(list(set([i.rsplit('_', 1)[0] for i in rep_cols])), reverse = True)
    fig, axs = plt.subplots(figsize = (width, 2.6)) 
    sns.boxplot(x = 'gene_res', y = 'value',  hue = 'variable', data = gene, boxprops=dict(alpha=.7), showfliers=False, hue_order = order, linewidth = 0.75, palette = ['b','g'], ax = axs)
    sns.stripplot(x= 'gene_res', y = 'value', hue = 'variable',linewidth = 0.5, edgecolor = 'gray', hue_order = order, data = gene, jitter = 0.2, dodge=True, palette = ['b','g'], ax = axs)
    plt.ylim(-4.5,4.5)
    handles, labels = axs.get_legend_handles_labels() 
    l = plt.legend(handles[0:2], ['Mitosis/Asynch','LPP(-)/LPP(+)'], bbox_to_anchor=(1, 0.625), loc=2, fontsize = 10, frameon = False)
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.ylabel('Mitosis/Asynch or LPP(-)/LPP(+)\nlog2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.xticks(ticks = [0,1,2,3,4], labels = ['MAP2K1/2 207/211','MAP2K3/6 207/196','MAP2K4 246', 'MAP2K5 300', 'MAP2K7 260'], fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    axs.get_legend().remove()
    plt.title(gene_name + ' aligned cysteines', family = 'ARIAL', fontsize = 10)
    plt.savefig('{}/figures/suppfig2/{}_suppfig2j_MAP2K_bar.pdf'.format(dir, date), bbox_inches='tight',transparent = True)
    plt.close()
cys_list = ['MAP2K2_211', 'MAP2K3_207','MAP2K4_246','MAP2K5_300','MAP2K7_260'] #'MAP2K1_207','MAP2K6_196'
# cys_list = ['MAP2K1_207','MAP2K2_211', 'MAP2K3_207','MAP2K4_246','MAP2K5_300','MAP2K6_196','MAP2K7_260'] #'MAP2K1_207','MAP2K6_196'
supp_fig2j_boxplots(cys_list, 2.4)




