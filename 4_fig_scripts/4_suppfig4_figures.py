#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 13:35:39 2021

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
from  matplotlib.ticker import PercentFormatter
import re

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

x_med ='Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
adapted_med ='Mitosis_LPP(-,+)/(+,+)_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
prots_df = pd.read_excel('{}/supplemental_tables/{}_fig2_proteins.xlsx'.format(dir,date), sheet_name = 'all prots with changes')
y_df = pd.read_excel('{}/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date), sheet_name = 'filtered')
y_exp = pd.read_excel('{}/supplemental_tables/{}_fig1_proteins.xlsx'.format(dir,date), sheet_name = 'all prots with changes')
lpp_df = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')

def supp_fig4b_pie(titles, sizes, sub):
    labels = 'Proteins with LPP-sensitive\nphosphorylation sites','Proteins with phosphosites\nnot quantified or unchanging'    
    total = sum(sizes)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.1,2)
    plt.title(titles, fontsize = 10)
    ax1.pie(sizes, autopct = (lambda p: '{:.0f}'.format(p * total / 100)) , startangle=90, colors = ['m','darkgray'])
    # ax1.pie(sizes, autopct = (lambda p: '{:.0f}/476'.format(p * total / 100)) , startangle=90)
    ax1.legend(labels = labels, loc='center left', bbox_to_anchor=(-0.2, -0.2),frameon=False, handletextpad=0.3, fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/suppfig4/{}_suppfig4b_pie_{}.pdf'.format(dir, date, sub), bbox_inches='tight',transparent = True)
    # plt.close()

def supp_fig4_pie():
    phospho_med = 'cysteine_' + x_med
    lpp_sens = lpp_df.loc[lpp_df[phospho_med] >= 2]
    df1 = df.loc[(df[x_med] <= 1/lpp_cutoff) | (df[x_med] >= lpp_cutoff)]
    df2 = df1.loc[(df[y_med] <= 1/asynch_cutoff) | (df1[y_med] >= asynch_cutoff)]
    down = df2.loc[(df2[x_med] <= 1/lpp_cutoff)]
    artifact = down.loc[down[adapted_med] > 1/unchanging_cutoff]
    
    all_changing = df2.loc[(df2[x_med] <= 1/lpp_cutoff) | (df2[x_med] >= lpp_cutoff)]
    no_ppep = all_changing.loc[~all_changing.accession_res.isin(artifact.accession_res)]

    p_all_df = df.loc[df.accession.isin(list(set(lpp_sens.accession)))]
    t_all = len(set(df.accession))
    p_all = len(set(p_all_df.accession))
    
    p_no_ppep_df = no_ppep.loc[no_ppep.accession.isin(list(set(lpp_sens.accession)))]
    t_no_ppep = len(set(no_ppep.accession))
    p_no_ppep = len(set(p_no_ppep_df.accession))
    
    sizes = [p_all, t_all - p_all]
    titles = 'All quantified proteins'
    sub = 'all_proteins'
    supp_fig4b_pie(titles, sizes, sub)

    titles1 = 'All proteins with phosphorylation-\ndependent cysteine reactivity changes'
    sizes1 = [p_no_ppep, t_no_ppep - p_no_ppep]
    sub1 = 'no_ppep'
    supp_fig4b_pie(titles1, sizes1, sub1)
supp_fig4_pie()
    
def supp_fig4a_prep():
    #proteins with high stoich phosphorylation, fig5 NUMBER IS DIFF FROM FIG 3 BC PROTEINS NOT sites
    df1 = df.loc[(df[y_med] >= asynch_cutoff) | (df[y_med] <= 1/asynch_cutoff)]
    up = df1.loc[df[x_med] <= 1/lpp_cutoff]
    # up = up.loc[~(up[adapted_med].isnull())]
    #write to source data
    # up.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig3f_pie.csv'.format(dir))
    authentic = up.loc[up[adapted_med] <= 1/adapted_cutoff]
    artifact = up.loc[up[adapted_med] > 1/unchanging_cutoff] 
    changing = df1.loc[(df1[x_med] >= lpp_cutoff) | (df[x_med] <= 1/lpp_cutoff)]
    am_diff = changing.loc[~changing.gene_res.isin(artifact)]
    #no phospho
    mann = 'mitosis_stoichiometry_mann'
    # olsen = 'mitosis_stoichiometry_olsen'
    #proteins that have quantified cysteines in mann or olsen in LPP changing proteins (not artifactual))
    no_mann = am_diff.loc[(am_diff[mann].isnull()) ]#& (am_diff[olsen].isnull())]
    overlap = am_diff.loc[~ (am_diff[mann].isnull())]# & (am_diff[olsen].isnull()))]
    #proteins that have quantified cysteines in mann or olsen in all proteins
    all_no_mann = df.loc[(df[mann].isnull()) ]#& (df[olsen].isnull())]
    all_overlap = df.loc[~ (df[mann].isnull()) ]#& (df[olsen].isnull()))]
    #number of proteins in mann, number not in mann
    auth_no = authentic.loc[(authentic[mann].isnull())]# & (authentic[olsen].isnull())]
    auth_overlap = authentic.loc[~ (authentic[mann].isnull()) ]#& (authentic[olsen].isnull()))]
    art_no = artifact.loc[(artifact[mann].isnull()) ]#& (artifact[olsen].isnull())]
    art_overlap = artifact.loc[~ (artifact[mann].isnull()) ]#& (artifact[olsen].isnull()))]
    overlap = len(set(overlap.accession))
    no_phos = len(set(no_mann.accession))
    all_overlap = len(set(all_overlap.accession))
    all_no = len(set(all_no_mann.accession))
    auth_overlap = len(set(auth_overlap.accession))
    auth_no = len(set(auth_no.accession))
    art_overlap = len(set(art_overlap.accession))
    art_no = len(set(art_no.accession))
    phospho_all = [all_overlap, all_no]
    lpp_change = [overlap, no_phos]
    auth_change = [auth_overlap, auth_no]
    art_change = [art_overlap, art_no]
    fig4_pies = pd.DataFrame({'all': phospho_all, 'lpp_change': lpp_change,'authentic': auth_change,'artifact': art_change })
    return fig4_pies

def supp_fig4a_pie(titles, sizes, sub):
    labels = 'Proteins with high stoichiometry\nphosphorylation in mitosis','Low or unquantified stoichiometry'    
    total = sum(sizes)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.1,2)
    plt.title(titles, fontsize = 10)
    ax1.pie(sizes, autopct = (lambda p: '{:.0f}'.format(p * total / 100)) , startangle=90)
    # ax1.pie(sizes, autopct = (lambda p: '{:.0f}/476'.format(p * total / 100)) , startangle=90)
    ax1.legend(labels = labels, loc='center left', bbox_to_anchor=(-0.2, -0.15),frameon=False, handletextpad=0.3, fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/suppfig4/{}_suppfig4a_pie_mann_only{}.pdf'.format(dir, date, sub), bbox_inches='tight',transparent = True)

def supp_fig4a_generation():
    fig4_pies = supp_fig4a_prep()
    sizes = fig4_pies['all']
    titles = 'All quantified proteins'
    sub = 'all_proteins'
    supp_fig4a_pie(titles, sizes, sub)
    titles = 'Proteins with authentic decrease\nin LPP(-, -) cysteine reactivity'
    sizes = fig4_pies['authentic']
    sub = 'authentic'
    supp_fig4a_pie(titles, sizes, sub)
    titles = 'Proteins with artifactual decrease\nin LPP(-, -) cysteine reactivity'
    sizes = fig4_pies['artifact']
    sub = 'artifact'
    supp_fig4a_pie(titles, sizes, sub)

supp_fig4a_generation()

def supp_fig4cd_moving_average(row):
    pos = int(row['pos']) - 1
    seq = row['sequence_fasta']       
    sp_loc = [abs(m.start() - pos) for m in re.finditer('SP|TP', seq)] 
    if not sp_loc:
        print('No S/T-P sites')
        # count0 = count0 + 1
        total_counts = [0] * 10 #+ [count0]
    else:
        lowest_loc = [min(sp_loc)]
        print('Closest SP or TP occurs for ' + str(row['gene_res']) + ' ' + str(lowest_loc) + ' away from Cys')
        count1 = len([i for i in lowest_loc if i < 6])
        count2 = len([i for i in lowest_loc if i >= 6 and i < 11])
        count3 = len([i for i in lowest_loc if i >= 11 and i < 16])
        count4 = len([i for i in lowest_loc if i >= 16 and i < 21])
        count5 = len([i for i in lowest_loc if i >= 21 and i < 26])
        count6 = len([i for i in lowest_loc if i >= 26 and i < 31])
        count7 = len([i for i in lowest_loc if i >= 31 and i < 36])
        count8 = len([i for i in lowest_loc if i >= 36 and i < 41])
        count9 = len([i for i in lowest_loc if i >= 41 and i < 46])
        count10 = len([i for i in lowest_loc if i >= 46 and i < 51])
        # count11 = len([i for i in lowest_loc if i >= 51])
        total_counts = [count1, count2, count3, count4, count5, count6, count7, count8, count9, count10]
    return total_counts

def supp_fig4cd_stp_prep(df): 
    fasta = pd.read_csv('{}/inputs/RS_2016_fasta_20210127.csv'.format(dir))
    fasta1 = fasta[['accession', 'sequence']]
    df1 = df.merge(fasta1, on = 'accession', how = 'left', suffixes = ['', '_fasta'])
    df1 = df1.assign(pos = df1.gene_res.str.split('_', n = 1, expand = True)[1])
    counts = df1.apply(supp_fig4cd_moving_average, axis =1).to_list()
    counts = pd.DataFrame(counts)
    print(len(df1.accession_res))
    percent = counts.sum(axis = 0)/len(df.accession_res)
    labels = ['1-5','6-10','11-15','16-20','21-25','26-30', '31-35', '36-40','41-45','46-50']
    graph_data = pd.DataFrame({'data': percent, 'label': labels})
    return graph_data

def supp_fig4cd_plot(df, title, sub):
    data = supp_fig4cd_stp_prep(df)
    fig, axs = plt.subplots(figsize=(2.1, 2.2))
    # fig, axs = plt.subplots()
    g = sns.barplot(x = 'label', y = 'data', color = 'black', data = data, linewidth = 0.75)
    plt.ylabel('Percent of cysteines',family = 'ARIAL', fontsize = 10)
    plt.xlabel('Closest S/T-P distance to cysteine',family = 'ARIAL', fontsize = 10)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 45, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.title(title, fontsize = 10, family = 'ARIAL')
    plt.ylim(0,0.7)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.savefig('{}/figures/suppfig4/{}_suppfig4cd_stp_{}.pdf'.format(dir,date, sub), bbox_inches='tight',transparent = True)
    # plt.close()
    
def supp_fig4cd_stp_generate():
    df1 = df.loc[~(df[y_med].isnull())]
    cond1 = (df1[x_med] <= (1/lpp_cutoff)) |  (df1[x_med] >= lpp_cutoff)
    cond2 = (df1[y_med] <= (1/asynch_cutoff)) |  (df1[y_med] >= asynch_cutoff)
    #changing in basal only
    am_diff = df1.loc[(cond1) & (cond2)]
    artifact = am_diff.loc[( (am_diff[x_med] <= 1/lpp_cutoff) & (am_diff[adapted_med] > 1/unchanging_cutoff))]
    authentic = am_diff.loc[( (am_diff[x_med] <= 1/lpp_cutoff) & (am_diff[adapted_med] <= 1/adapted_cutoff))]
    lpp_down = am_diff.loc[am_diff[x_med] <= 1/lpp_cutoff]
    no_ppep = lpp_down.loc[~(lpp_down[adapted_med] > 1/unchanging_cutoff)]
    lpp_up = am_diff.loc[am_diff[x_med] >= lpp_cutoff]
    supp_fig4cd_plot(df1, 'All quantified cysteines', 'all')
    supp_fig4cd_plot(no_ppep, 'Decrease', 'down')
    supp_fig4cd_plot(lpp_up, 'Increase', 'up')
    supp_fig4cd_plot(authentic, 'Less reactive in LPP(-, -):  authentic', 'authentic')
    supp_fig4cd_plot(artifact, 'Less reactive in LPP(-, -):  artifact', 'artifact')
supp_fig4cd_stp_generate()

def supp_fig4g_boxplots(protein_name):
    gene1 = df.loc[df.protein == protein_name]
    sdf = gene1.copy()
    drop_cols = [col for col in sdf.columns if 'category' in col] + [col for col in sdf.columns if 'protein_' in col] + [col for col in sdf.columns if 'limited' in col]
    sdf.drop(drop_cols, axis =1, inplace = True)
    sdf.set_index(['accession', 'description', 'protein','asynch_stoichiometry_mann','mitosis_stoichiometry_mann','sequence','accession_res','gene_res'], inplace = True)
    median_cols = [col for col in sdf.columns if 'median' in col]
    sdf[median_cols] = np.log2(sdf[median_cols])
    sdf.to_excel('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig4g_{}.xlsx'.format(dir, protein_name),  float_format = '%.2f')
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    gene1 = gene1.loc[(~gene1[x_med].isnull())& (~gene1[y_med].isnull())]
    adapted_name = adapted_med.rsplit(sep = '_', maxsplit = 2)[0]
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
    fig, axs = plt.subplots( figsize=(4, 2.4))
    order = [y_name, x_name, adapted_name]
    sns.stripplot(x = 'gene_res', y = 'value', hue = 'variable', linewidth = 0.5, edgecolor = 'gray',\
                  data = gene, jitter = True, dodge = True, palette = ['b','g','#daa520'], hue_order = order, ax = axs)
    sns.boxplot(x = 'gene_res', y = 'value', hue = 'variable', data = gene, palette = ['b','g','#daa520'], \
                hue_order = order, boxprops=dict(alpha=.7), linewidth = 0.75, showfliers=False, ax = axs) 
    plt.ylabel('log2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    handles, labels = axs.get_legend_handles_labels() 
    l = plt.legend(handles[0:3], ['Mitosis/Asynch','LPP(-, -)/LPP(+, -)','LPP(-, +)/LPP(+, +)'], bbox_to_anchor=(1, 0.71), loc=2, fontsize = 10, frameon = False)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 45, ha = 'right')
    plt.title(protein_name, family = 'ARIAL', fontsize = 10)
    plt.ylim(-4.5,4.5)
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.savefig('{}/figures/suppfig4/{}_suppfig4g_bar_{}.pdf'.format(dir, date, protein_name), bbox_inches='tight',transparent = True)
    plt.close()
supp_fig4g_boxplots('KLC1')
        
def supp_fig4h_bar(df, gene_res):
    gene1 = df.loc[df.gene_res == gene_res]
    x_med1 = 'cysteine_' + x_med
    sdf = gene1.copy()
    drop_cols = [col for col in sdf.columns if 'category' in col] + [col for col in sdf.columns if 'protein_' in col]
    sdf.drop(drop_cols, axis =1, inplace = True)
    sdf.columns = sdf.columns.str.replace('cysteine_','peptide_')
    median_cols = [col for col in sdf.columns if 'median' in col]
    sdf[median_cols] = 1/sdf[median_cols]
    colors = ['m','m','m']
    sdf.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig4h_phospho_bar_{}.csv'.format(dir, gene_res), index = False, float_format = '%.2f' )
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
    order = [ y_name, 'Mitosis LPP(-)', x_name]
    fig, axs = plt.subplots(figsize = (1.5,2.4)) #figsize=(2.8, 2.8)
    # sns.stripplot(x= 'variable1', y = 'value1', data = gene, jitter = jit_num, linewidth = 0.5, order = order, dodge=True, ax = axs, palette = colors, edgecolor = ecol)
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
    plt.savefig('{}/figures/suppfig2/{}_suppfig4h_bar_{}.pdf'.format(dir, date,  gene_res1), bbox_inches='tight',transparent = True)
supp_fig4h_bar(lpp_df, 'FXR2_450')

def supp_fig4i_boxplots(gene_cys):
    gene1 = df.loc[df.gene_res == gene_cys]
    sdf = gene1.copy()
    drop_cols = [col for col in sdf.columns if 'limited' in col] + [col for col in sdf.columns if 'category' in col] + [col for col in sdf.columns if 'protein_' in col]
    sdf = sdf.drop(drop_cols, axis = 1)
    sdf.set_index(['accession', 'description', 'protein','asynch_stoichiometry_mann','mitosis_stoichiometry_mann','sequence','gene_res'], inplace = True)
    median_cols = [col for col in sdf.columns if 'median' in col]
    sdf[median_cols] = np.log2(sdf[median_cols])
    sdf.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fsuppfig4i_box_{}.csv'.format(dir, gene_cys), float_format = '%.2f')
    gene_res = gene_cys.split('_')[0] + ' C' + gene_cys.split('_')[1]
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    adapted_name = adapted_med.rsplit(sep = '_', maxsplit = 2)[0]
    adapted_cols = [col for col in gene1.columns if adapted_name + '_median_TMT' in col] 
    colors = ['#404040', 'blue', 'blue']
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
    plt.savefig('{}/figures/suppfig4/{}_suppfig4i_bar_{}.pdf'.format(dir,date, gene_cys), bbox_inches='tight',transparent = True)
    plt.close()
supp_fig4i_boxplots('FXR2_270')



