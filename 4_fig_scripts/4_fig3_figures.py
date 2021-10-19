#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 11:47:33 2021

@author: estherkemper
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

protein_cols = ['accession','description', 'protein', 'gene_res']
phospho_col = ['asynch_stoichiometry_mann','mitosis_stoichiometry_mann']

dir = os.getcwd()
date = '20210830'

x_med ='Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
adapted_med ='Mitosis_LPP(-,+)/(+,+)_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir,date), sheet_name = 'filtered')
print('Finished reading input data')

def fig3a_STP_bar():
    motif = 'SP|TP'
    df1 = df.loc[~df[y_med].isnull()].reset_index(drop = True)
    low = df1.loc[(df1[x_med] <= 1/lpp_cutoff) & ((df1[y_med] >= asynch_cutoff) | (df1[y_med] <= (1/asynch_cutoff)))]
    high = df1.loc[(df1[x_med] >= lpp_cutoff) & ((df1[y_med] >= asynch_cutoff) | (df1[y_med] <= (1/asynch_cutoff)))]
    labels = ['More','Less', 'All']
    total = [len(set(high.accession_res)), len(set(low.accession_res)), len(set(df1.accession_res))]
    high_motif = high.loc[high.sequence.str.contains(motif)]
    low_motif = low.loc[low.sequence.str.contains(motif)]
    all_motif = df1.loc[df1.sequence.str.contains(motif)]
    counts = [len(set(high_motif.accession_res)), len(set(low_motif.accession_res)), len(set(all_motif.accession_res))]
    y_value =  [i / j for i, j in zip(counts, total)] 
    percent_labels = [str(i)+ '/\n' + str(j) for i,j in zip(counts, total)]    
    plt.subplots(figsize=(2.2,2.5))
    ax = sns.barplot(x = labels, y = y_value, palette = ('blue','red','gray')) ##daa520
    plt.ylabel('Percent peptides containing\nS/T-P substrate motif', fontsize =10, family = 'ARIAL')
    plt.xticks(fontsize =10, family = 'ARIAL')
    plt.yticks(fontsize =10, family = 'ARIAL')
    plt.ylim(0,1)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    for i in range(len(percent_labels)):
        plt.annotate(str(percent_labels[i]), xy=(i,y_value[i]), ha='center', va='bottom', fontsize = 10, family = 'ARIAL')
    plt.savefig('{}/figures/fig3/{}_fig3a_STP_bar.pdf'.format(dir, date), bbox_inches='tight')
    keep_cols = [col for col in df1.columns if 'combined_median' in col]
    df2 = df1[protein_cols +['accession_res','sequence'] + keep_cols]
    df2.loc[df2.accession_res.isin(list(set(low.accession_res))), 'Reactivity in LPP(-)'] = 'Down'
    df2.loc[df2.accession_res.isin(list(set(high.accession_res))), 'Reactivity in LPP(-)'] = 'Up'
    df2.loc[df2.sequence.str.contains(motif), 'Contain S/T-P on '] = 'True'
    df2.loc[~df2.sequence.str.contains(motif),'Contain S/T-P on '] = 'False'
    median_cols = [col for col in df2.columns if 'median' in col]
    df2[median_cols] = np.log2(df2[median_cols])
    df2.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig3a_stp.csv'.format(dir), float_format = '%.2f', index = False)
    plt.close()
fig3a_STP_bar()


def fig3b_protein_boxplots(protein_name):
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    gene1 = df.loc[df.protein == protein_name]
    sdf = gene1.copy()
    sdf.set_index(['accession', 'description', 'protein','asynch_stoichiometry_mann','mitosis_stoichiometry_mann','sequence','accession_res','gene_res'], inplace = True)
    keep_cols = [col for col in sdf.columns if x_name in col]
    sdf = sdf[keep_cols]
    median_cols = [col for col in sdf.columns if 'median' in col]
    cat_cols = [col for col in sdf.columns if 'category' in col]
    sdf[median_cols] = np.log2(sdf[median_cols])
    sdf.drop(cat_cols, axis =1 , inplace = True)
    sdf = sdf.sort_values('protein', ascending = True)
    
    sdf.to_excel('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig3b_{}.xlsx'.format(dir, protein_name),  float_format = '%.2f')
    gene1.loc[(gene1[x_med] >= 2), 'hue'] = 'blue'
    gene1.loc[(gene1[x_med] <= 0.5), 'hue'] = 'red'
    gene1.loc[(gene1[x_med] > 0.5) & (gene1[x_med] < 2), 'hue'] = '#404040'
    gene_name = gene1.iloc[0,:].gene_res.split('_')[0]
    gene1 = gene1.assign(pos = gene1.gene_res.str.split('_', expand = True)[1].astype(int))
    gene1 = gene1.sort_values('pos', ascending = True)     
    hue_order = list(gene1.hue)
    order = list('C' + gene1.pos.astype(str))     
    gene_name = gene1.iloc[0,:].gene_res.split('_')[0]
    
    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col]
    gene = gene1[['protein', 'gene_res'] + rep_cols]
    gene = gene.melt(id_vars = ['protein', 'gene_res']).dropna(axis = 0)
    gene = gene.assign(variable = gene.variable.str.rsplit(pat = '_', n = 1, expand = True)[0])
    gene = gene.assign(value = np.log2(gene.value))
    gene =gene.sort_values('gene_res', ascending = False)
    gene = gene.assign(gene_res = 'C' + gene.gene_res.str.split('_', n =1, expand = True)[1])
    fig, axs = plt.subplots(figsize=(3, 2))
    error_kw=dict(lw=1, capsize=5, capthick=1)
    sns.boxplot(x = 'gene_res', y = 'value',  data = gene, width = 0.5, palette = hue_order, order = order, boxprops=dict(alpha=.7), showfliers=False, linewidth = 0.75)
    sns.stripplot(x = 'gene_res', y = 'value', jitter = 0.2, linewidth = 0.5, edgecolor = 'black', data = gene, palette = hue_order, order = order)
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.ylim(-5, 2)
    plt.ylabel('log2 LPP(-)/LPP(+)\nMS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xlabel('',family = 'ARIAL', fontsize = 10)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 45, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.title(gene_name, family = 'ARIAL', fontsize = 10)
    plt.savefig('{}/figures/fig3/{}_fig3b_bar_{}.pdf'.format(dir,date, gene_name), bbox_inches='tight',transparent = True)
    plt.close()
fig3b_protein_boxplots('EDC3')
fig3b_protein_boxplots('GTF2I')
    
def fig3e_correlation_plot():
    #source data is main dataset
    df1= df.loc[~df[adapted_med].isnull()]
    df1 = df1.loc[~df1[x_med].isnull()]
    df1 = df1[['gene_res', 'accession', x_med,adapted_med]]
    df1[[x_med,adapted_med]] = np.log2(df1[[x_med,adapted_med]] )
    df1 = df1.sort_values(adapted_med, ascending = True)
    plt.subplots(figsize=(2.8,2.8))
    y = df1[adapted_med]
    x = df1[x_med]
    cond1 = (df1[adapted_med] >= np.log2(adapted_cutoff)) & (df1[x_med] >= np.log2(lpp_cutoff))
    cond2 = (df1[adapted_med] <= np.log2(1/adapted_cutoff)) & (df1[x_med] <= np.log2(1/lpp_cutoff))
    cond4 = (df1[adapted_med] >= np.log2(adapted_cutoff)) & (df1[x_med] <= np.log2(1/lpp_cutoff)) 
    cond5 = (df1[adapted_med] > np.log2(1/adapted_cutoff)) & (df1[x_med] <= np.log2(1/lpp_cutoff))
    up_mask = (cond1) #| (cond3) 
    down_mask = (cond2)# | (cond4)
    both = (up_mask) | (down_mask) | (cond5)
    plt.scatter(x[up_mask], y[up_mask], lw = 0.1, alpha = 0.4,edgecolors = 'black',c = 'blue', marker = '.',s=80)
    plt.scatter(x[down_mask], y[down_mask], lw = 0.1, alpha = 0.4,edgecolors = 'black',c = 'red', marker = '.',s=80)
    plt.scatter(x[cond5], y[cond5], lw = 0.1, alpha = 0.4,edgecolors = 'black',c = '#daa520', marker = '.',s=80)
    plt.scatter(x[~both], y[~both], lw = 0.1, alpha= 0.4, edgecolors = 'black',c = '#404040', marker = '.',s=80)
    plt.axhline(y = np.log2(adapted_cutoff), c = 'black', linestyle = 'dotted')
    plt.axhline(y = np.log2(1/adapted_cutoff), c = 'black', linestyle = 'dotted')
    plt.axvline(x = 1, c = 'black', linestyle = 'dotted')
    plt.axvline(x = -1, c = 'black', linestyle = 'dotted')
    plt.ylabel('Adapted protocol\nlog2 LPP(-, +)/LPP(+, +)', fontsize =10, family = 'ARIAL')
    plt.xlabel('Original protocol\nlog2 LPP(-, -)/LPP(+, -)',fontsize =10, family = 'ARIAL')
    plt.ylim(-4.5,4.5)
    plt.xlim(-4.5,4.5) 
    plt.xticks(fontsize =10, family = 'ARIAL',ticks= [-4,-2,0,2,4])
    plt.yticks(fontsize =10, family = 'ARIAL',ticks= [-4,-2,0,2,4])
    plt.savefig('{}/figures/fig3/{}_fig3e_correlation1.pdf'.format(dir,date), bbox_inches='tight')
    plt.close()
fig3e_correlation_plot()

def fig3f_pie(df):
    #cysteines that still have LPP/CTRL changes after LPP-after treatment
    df1 = df.loc[(df[y_med] >= asynch_cutoff) | (df[y_med] <= 1/asynch_cutoff)]
    up = df1.loc[df1[x_med] <= 1/lpp_cutoff]
    up = up.loc[~(up[adapted_med].isnull())]
    up['category'] = 'Uncategorized'
    up.loc[up[adapted_med] <= 1/adapted_cutoff, 'category'] = 'Authentic'
    up.loc[up[adapted_med] > 1/unchanging_cutoff,'category'] = 'Artifact'
    #write to source data
    keep_cols = [col for col in up.columns if 'combined_median' in col]
    source_up = up[protein_cols +['accession_res','sequence'] + keep_cols + ['category']]
    median_cols = [col for col in source_up.columns if 'median' in col]
    source_up[median_cols] = np.log2(source_up[median_cols])
    drop_cols = [col for col in source_up.columns if 'limited' in col] + [col for col in source_up.columns if 'protein_' in col]
    source_up.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig3f_pie.csv'.format(dir), index = False, float_format = '%.2f')
    authentic = up.loc[up[adapted_med] <= 1/adapted_cutoff]
    artifact = up.loc[up[adapted_med] > 1/unchanging_cutoff]  
    uncategorized = up.loc[(up[adapted_med] > 1/adapted_cutoff) & (up[adapted_med] <= 1/unchanging_cutoff)]
    labels = 'Authentic', 'Artifact','Uncategorized'
    sizes = [len(set(authentic.accession_res)), len(set(artifact.accession_res)), len(set(uncategorized.accession_res))]
    total = sum(sizes)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.3,2.3)
    ax1.pie(sizes, autopct = (lambda p: '{:.0f}'.format(p * total / 100)) , colors = ['red','#daa520','lightgray'], startangle=24)
    ax1.legend( labels = labels, loc='center left', bbox_to_anchor=(0.1, -0.3),frameon=False, handletextpad=0.3, fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/fig3/{}_fig3f_pie.pdf'.format(dir, date), bbox_inches='tight',transparent = True)
    plt.close()
fig3f_pie(df)
    
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
                colors = ['#404040','#404040','blue']
                return colors
            
def fig3ghk_boxplots(gene_cys):
    gene1 = df.loc[df.gene_res == gene_cys]
    sdf = gene1.copy()
    drop_cols = [col for col in sdf.columns if 'limited' in col] + [col for col in sdf.columns if 'category' in col] + [col for col in sdf.columns if 'protein_' in col]
    sdf = sdf.drop(drop_cols, axis = 1)
    sdf.set_index(['accession', 'description', 'protein','asynch_stoichiometry_mann','mitosis_stoichiometry_mann','sequence','gene_res'], inplace = True)
    median_cols = [col for col in sdf.columns if 'median' in col]
    sdf[median_cols] = np.log2(sdf[median_cols])
    sdf.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig3ghk_box_{}.csv'.format(dir, gene_cys), float_format = '%.2f')
    gene_res = gene_cys.split('_')[0] + ' C' + gene_cys.split('_')[1]
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    adapted_name = adapted_med.rsplit(sep = '_', maxsplit = 2)[0]
    adapted_cols = [col for col in gene1.columns if adapted_name + '_median_TMT' in col] 
    colors = color_prep(gene1)
    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] +adapted_cols +\
        [col for col in gene1.columns if y_name + '_median_TMT' in col] 
    rep_cols = [col for col in rep_cols if 'protein_' not in col] 
    gene = gene1[['gene_res'] + rep_cols]    
    gene = gene.melt(id_vars = 'gene_res').dropna(axis = 0)
    gene = gene.assign(variable1 = gene.variable.str.rsplit(pat = '_', n = 3, expand = True)[0])
    gene = gene.assign(value = np.log2(gene.value))
    gene = gene.sort_values('variable1', ascending = True)
    fig, axs = plt.subplots( figsize=(1.2, 2.4))
    error_kw=dict(lw=1, capsize=3, capthick=1)
    sns.boxplot(x = 'variable1', y = 'value',  data = gene, palette = colors, \
                order = [y_name, x_name, adapted_name], \
                boxprops=dict(alpha=.7), linewidth = 0.75, showfliers=False) 
    sns.stripplot(x = 'variable1', y = 'value', data = gene, linewidth = 0.5, edgecolor = 'black', jitter = 0.2, palette = colors, order =  [y_name, x_name, adapted_name])
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
    plt.savefig('{}/figures/fig3/{}_fig3ghk_bar_{}.pdf'.format(dir,date, gene_cys), bbox_inches='tight',transparent = True)
    plt.close()
cys_list = ['MAP2K1_277', 'PIAS1_481', 'FLNA_1453']
for cys in cys_list:
    fig3ghk_boxplots(cys)
    
    
def fig3j_STP_bar():
    motif = 'SP|TP'
    df1 = df.loc[~df[y_med].isnull()]
    low = df1.loc[(df1[x_med] <= 1/lpp_cutoff) & ((df1[y_med] >= asynch_cutoff) | (df1[y_med] <= (1/asynch_cutoff)))]
    #source data
    keep_cols = [col for col in low.columns if '_combined' in col]
    sdf = low[protein_cols +['sequence'] + keep_cols]
    sdf = sdf.loc[sdf[adapted_med].notnull()]
    sdf.loc[sdf.sequence.str.contains(motif), 'Contains S/T-P on peptide'] = 'True'
    sdf.loc[~sdf.sequence.str.contains(motif), 'Contains S/T-P on peptide'] = 'False'
    sdf.loc[sdf[adapted_med].notnull(),'Authentic or Artifact'] = 'Unassigned'
    sdf.loc[sdf[adapted_med] > 1/unchanging_cutoff,'Authentic or Artifact'] = 'Artifact'
    sdf.loc[sdf[adapted_med] <= 1/adapted_cutoff, 'Authentic or Artifact'] = 'Authentic'
    keep_cols = [col for col in sdf.columns if 'combined_median' in col]
    median_cols = [col for col in sdf.columns if 'median' in col]
    sdf[median_cols] = np.log2(sdf[median_cols])
    drop_cols = [col for col in sdf.columns if 'limited' in col] + [col for col in sdf.columns if 'protein_' in col]
    sdf.drop(drop_cols, axis = 1, inplace = True)
    sdf.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig3j_stp.csv'.format(dir), index = False, float_format = '%.2f')
    ppep = low.loc[low[adapted_med] > 1/unchanging_cutoff]
    changing = low.loc[low[adapted_med] <= 1/adapted_cutoff]
    ppep_motif = ppep.loc[ppep.sequence.str.contains(motif)]
    changing_motif = changing.loc[changing.sequence.str.contains(motif)]
    labels = ['Authentic','Artifact']
    counts = [ len(set(changing_motif.accession_res)),len(set(ppep_motif.accession_res))]
    total = [len(set(changing.accession_res)), len(set(ppep.accession_res))]
    y_value =  [i / j for i, j in zip(counts, total)] 
    percent_labels = [str(i)+ '/\n' + str(j) for i,j in zip(counts, total)]    
    plt.subplots(figsize=(1.6,2.5))
    ax = sns.barplot(x = labels, y = y_value, palette = ('red','#daa520'))
    plt.ylabel('Percent peptides containing\nS/T-P substrate motif', fontsize =10, family = 'ARIAL')
    plt.xticks(fontsize =10, family = 'ARIAL')
    plt.yticks(fontsize =10, family = 'ARIAL')
    plt.ylim(0,1)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    for i in range(len(percent_labels)):
        plt.annotate(str(percent_labels[i]), xy=(i,y_value[i]), ha='center', va='bottom', fontsize = 10, family = 'ARIAL')
    plt.savefig('{}/figures/fig3/{}_fig3j_STP_bar.pdf'.format(dir, date), bbox_inches='tight')
    plt.close()
fig3j_STP_bar()



    
    
    
    