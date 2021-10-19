#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 07:28:02 2021

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

x_med ='Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
adapted_med ='Mitosis_LPP(-,+)/(+,+)_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
prots_df = pd.read_excel('{}/supplemental_tables/{}_fig2_proteins.xlsx'.format(dir,date), sheet_name = 'all prots with changes')
y_df = pd.read_excel('{}/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date), sheet_name = 'filtered')
y_exp = pd.read_excel('{}/supplemental_tables/{}_fig1_proteins.xlsx'.format(dir,date), sheet_name = 'all prots with changes')

lpp_df = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')

print('Finished reading input data')

def fig2b_heatmap(df, x_med, y_med):
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
    plt.figure(figsize = (1,3.4))
    # cmap = ['#990000','#ff3333','#ffcccc','white', '#e6e6ff', '#4d4dff', '#000099']
    # cmap = sns.light_palette((240, 100, 30), n_colors = 10, as_cmap=True, input="husl")
    cmap = sns.light_palette('#000099', as_cmap = True)
    ax = sns.heatmap(df3, cmap = cmap, vmin = 0,xticklabels=False)
    colorbar = ax.collections[0].colorbar
    colorbar.remove()
    ax.figure.axes[-1].yaxis.label.set_size(10)
    ax.figure.axes[-1].tick_params(labelsize=10)
    ax.set_ylabel("Phosphopeptides", fontsize = 10, family = 'ARIAL')
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 90, ha = 'right')
    plt.yticks(np.arange(0,11000,1000).tolist(),np.arange(0,11000,1000).tolist() , fontsize = 10, family = 'ARIAL')
    plt.title('Phosphoproteomics', family = 'ARIAL', fontsize = 10)
    # ax.vlines([0.97,2], ymin = 0, ymax = 13797, color = 'black', lw = 1)
    plt.savefig('{}/figures/fig2/{}_fig2b_phosphoenrichment_heatmap.pdf'.format(dir,date), transparent = True, bbox_inches='tight')
    plt.close()
fig2b_heatmap(lpp_df,'cysteine_' + x_med , y_med)

def fig2d_waterfall(df, label):
    #source data main dataset
    print('Running Fig2 waterfall plot')
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
    plt.scatter(x[up_cond], y[up_cond], lw = 0, c = 'blue', marker = 'o', s = 8, alpha = 0.6)
    plt.scatter(x[down_cond], y[down_cond], lw = 0, c = 'red', marker = 'o', s = 8, alpha = 0.6)
    plt.xlabel('Peptides', fontsize =10, family = 'ARIAL')
    plt.ylabel(label + 'log2 MS3 reporter ion intensity',fontsize =10, family = 'ARIAL')
    plt.xticks(fontsize =10, family = 'ARIAL')
    plt.yticks(fontsize =10, family = 'ARIAL', ticks= [-4,-2,0,2,4])
    plt.axhline(y = 1, c = 'black', linestyle = 'dotted')
    plt.axhline(y = -1, c = 'black', linestyle = 'dotted')
    plt.ylim(-5,5) 
    title = x_med.split('/')[0]
    plt.savefig('{}/figures/fig2/{}_fig2d_waterfall.pdf'.format(dir,date), bbox_inches='tight')
    plt.close()
lpp_label = 'Mitosis LPP(-)/LPP(+)\n'
fig2d_waterfall(df, lpp_label)

def fig2e_correlation_plot():
    #source data main dataset
    print('Running Fig2 correlation plot')
    plt.subplots(figsize=(2.8,2.8)) 
    df1= df.loc[~df[y_med].isnull()]
    df1 = df1.loc[~df1[x_med].isnull()]
    df1 = df1[['gene_res', 'accession', x_med,y_med]]
    df1[[x_med,y_med]] = np.log2(df1[[x_med,y_med]] )
    x = df1[x_med]
    y = df1[y_med]
    cond1 = (df1[y_med] >= np.log2(protein_cutoff)) & (df1[x_med] >= np.log2(lpp_cutoff))
    cond2 = (df1[y_med] <= np.log2(1/protein_cutoff)) & (df1[x_med] <= np.log2(1/lpp_cutoff))
    cond3 = (df1[y_med] <= np.log2(1/protein_cutoff)) & (df1[x_med] >= np.log2(lpp_cutoff))
    cond4 = (df1[y_med] >= np.log2(protein_cutoff)) & (df1[x_med] <= np.log2(1/lpp_cutoff))
    up_mask = (cond1) | (cond3) 
    down_mask = (cond2) | (cond4)
    both = (up_mask) | (down_mask)
    plt.scatter(x[~both], y[~both], lw = 0, alpha= 0.5, c = '#404040', marker = '.')
    plt.scatter(x[up_mask], y[up_mask], lw = 0, alpha = 0.5,c = 'blue', marker = '.')
    plt.scatter(x[down_mask], y[down_mask], lw = 0, alpha = 0.5,c = 'red', marker = '.')
    plt.ylabel('Mitosis/Asynch\nlog2 MS3 reporter ion intensity', fontsize =10, family = 'ARIAL')
    plt.xlabel('Mitosis LPP(-)/LPP(+)\nlog2 MS3 reporter ion intensity',fontsize =10, family = 'ARIAL')
    plt.axhline(y = np.log2(protein_cutoff), xmin = 0.611, xmax = 1, c = 'black', linestyle = 'dotted')
    plt.axhline(y = -np.log2(protein_cutoff), xmin = 0, xmax = 0.388, c = 'black', linestyle = 'dotted')
    plt.axhline(y = np.log2(protein_cutoff), xmin = 0, xmax = 0.388, c = 'black', linestyle = 'dotted')
    plt.axhline(y = -np.log2(protein_cutoff), xmin = 0.611, xmax = 1, c = 'black', linestyle = 'dotted')
    plt.axvline(x = 1, ymin = 0, ymax = 0.464, c = 'black', linestyle = 'dotted')
    plt.axvline(x = -1, ymin =0.536, ymax = 1, c = 'black', linestyle = 'dotted')
    plt.axvline(x = 1, ymin =0.536, ymax = 1, c = 'black', linestyle = 'dotted')
    plt.axvline(x = -1, ymin = 0, ymax = 0.434, c = 'black', linestyle = 'dotted')
    plt.ylim(-4.5,4.5)
    plt.xlim(-4.5,4.5) 
    plt.xticks(fontsize =10, family = 'ARIAL',ticks= [-4,-2,0,2,4])
    plt.yticks(fontsize =10, family = 'ARIAL',ticks= [-4,-2,0,2,4])
    plt.savefig('{}/figures/fig2/{}_fig2e_correlation.pdf'.format(dir,date), bbox_inches='tight')
    plt.close()
fig2e_correlation_plot()

def fig2f_venn_prep():
    #source data generated here
    x_name = x_med.rsplit(sep = '_', maxsplit = 1)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 1)[0]
    x_no = x_name + '_no'
    y_no = y_name + '_no'
    x_df = df[[ 'accession_res', 'gene_res', x_med, x_no]]
    y_df1 = y_df[[ 'accession_res', 'gene_res',y_med, y_no]]
    y_exp1 = y_exp[['accession_res','gene_res', 'exp']]
    y_df1 = y_df1.merge(y_exp1, how = 'left', on = [ 'accession_res','gene_res'])
    all_df = x_df.merge(y_df1, how = 'outer', on = ['accession_res','gene_res'])
    all_df = all_df[(all_df[x_no] >= 2) & (all_df[y_no] >=2)]
    all_df1 = all_df.loc[~(all_df.exp == 'EXPRESSION')]
    x_diff = all_df1.loc[(all_df1[x_med] >= 2) | (all_df1[x_med] <= 0.5)]
    y_diff = all_df1.loc[(all_df1[y_med] >= asynch_cutoff) | (all_df1[y_med] <= (1/asynch_cutoff))]
    both =  all_df1.loc[(all_df1.accession_res.isin(list(set(y_diff.accession_res)))) & (all_df1.accession_res.isin(list(set(x_diff.accession_res))))]
    only_x = list(set(set(x_diff.accession_res) - set(both.accession_res)))
    only_y = list(set(y_diff.accession_res) - set(both.accession_res))
    total = list(set(both.accession_res))
    all_df.loc[all_df.accession_res.isin(only_x), 'Changing LPP(-)?LPP(+) or Mitosis/Asynch?'] = 'LPP'
    all_df.loc[all_df.accession_res.isin(only_y), 'Changing LPP(-)?LPP(+) or Mitosis/Asynch?'] = 'Asynch'
    all_df.loc[all_df.accession_res.isin(total), 'Changing LPP(-)?LPP(+) or Mitosis/Asynch?'] = 'Both'
    all_df.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig2f.csv'.format(dir), index = False)
    foo = all_df.loc[all_df.duplicated()]
    sizes = [ len(only_y), len(only_x),len(both)]
    plt.figure(figsize=(2.8, 2.8))
    out = venn2(subsets = sizes, set_labels = ('Mitosis/Asynch', 'LPP(-)/LPP(+)'), set_colors = ('b','g'), alpha = 0.9)
    for text in out.set_labels:
        text.set_fontsize(10)
    for text in out.subset_labels:
        text.set_fontsize(10)
        plt.savefig('{}/figures/fig2/{}_fig2f_venn.pdf'.format(dir,date), bbox_inches= 'tight')
        # plt.close()
fig2f_venn_prep()

def fig2g_pie_source_data(df):
    #get nan in the index from excel
    df1 = df.copy()
    df1[['accession', 'protein']] = df[['accession','protein']].fillna(method = 'ffill')
    df1['react'] = df1['react'].replace('Unassigned', np.nan)
    df1['react'] = df1['react'].replace('Unclear', np.nan)
    df1['react'] = df1['gene_res'] + ': ' + df1['react']
    df1['exp'] = df1['gene_res'] + ': ' + df1['exp']
    df1['cys_category'] = df1[['react','exp']].apply(lambda x: x.str.cat(sep=''), axis=1)
    df2 = df1[['accession','protein','cys_category']]
    df3 = df2.loc[~df2.cys_category.str.contains('Unassigned')]
    df3 = df3.loc[~(df3['cys_category'] == '')]
    df4 = pd.DataFrame(df3.groupby(['accession', 'protein'])['cys_category'].sum()).reset_index()
    df4.loc[df4.cys_category.str.contains('Unclear', regex = True),'SITE SPECIFIC?'] = 'Unclear'
    df4.loc[df4.cys_category.str.contains('EXPRESSION'),'SITE SPECIFIC?'] = 'Not site specific'
    df4.loc[df4.cys_category.str.contains('REACTIVITY'),'SITE SPECIFIC?'] = 'Site specific'
    df4.drop('cys_category', axis = 1, inplace = True)
    #save for source data
    df4.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig2g.csv'.format(dir), index = False)
    return df4

def fig2g_pie(df):
    df4 = fig2g_pie_source_data(df)
    react = df4.loc[(df4['SITE SPECIFIC?'].str.contains('Site specific'))]
    exp = df4.loc[(df4['SITE SPECIFIC?'].str.contains('Not site specific'))]
    unclear = df4.loc[(df4['SITE SPECIFIC?'].str.contains('Unclear'))]
    
    labels = 'Site-specific','Not site-specific','Unclear'
    sizes = [len(set(react.accession)),len(set(exp.accession)), len(set(unclear.accession)) ]

    total = sum(sizes)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.4,2.4)
    ax1.pie(sizes, autopct = (lambda p: '{:.0f}'.format(p * total / 100)) , colors = \
        [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),\
            (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),\
                (0.8666666666666667, 0.5176470588235295, 0.3215686274509804)], textprops={'fontsize': 10},startangle=0)
    ax1.legend(labels = labels, loc='center left', bbox_to_anchor=(0.1, -0.35),frameon=False, handletextpad=0.3, fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/fig2/{}_fig2g_pie.pdf'.format(dir, date), bbox_inches='tight',transparent = True)
    plt.close()
fig2g_pie(prots_df)

def final_fig2hij_boxplots(protein_name, width):
    gene1 = df.loc[df.protein == protein_name]
    median_cols = [col for col in gene1.columns if 'median' in col]
    cat_cols = [col for col in gene1.columns if 'category' in col]
    sdf = gene1.copy()
    sdf[median_cols] = np.log2(sdf[median_cols])
    sdf.drop(cat_cols, axis =1 , inplace = True)
    sdf = sdf.sort_values('protein', ascending = True)
    sdf.set_index(['accession', 'description', 'protein','sequence','accession_res','gene_res'], inplace = True)
    sdf.to_excel('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig2hij_{}.xlsx'.format(dir, protein_name), float_format = '%.2f')
    gene_name = gene1.iloc[0,:].gene_res.split('_')[0]
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
    if len(set(gene1.gene_res)) > 9:
        x_no = x_name + '_combined_no'
        y_no = y_name + '_combined_no'
        gene1 = gene1.loc[(gene1[x_no]>=8) & (gene1[y_no] >= 8)]
    else:
        pass
    gene1 = gene1.loc[(~gene1[x_med].isnull())& (~gene1[y_med].isnull())]
    rep_cols = [col for col in gene1.columns if x_name + '_median_TMT' in col] + [col for col in gene1.columns if y_name + '_median_TMT' in col]
    rep_cols = [col for col in rep_cols if 'protein_' not in col]
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
    fig, axs = plt.subplots(figsize = (width, 2.8)) #figsize=(2.8, 2.8)
    sns.boxplot(x = 'gene_res', y = 'value',  hue = 'variable', data = gene, boxprops=dict(alpha=.7), showfliers=False, hue_order = order, linewidth = 0.75, palette = ['b','g'], ax = axs)
    sns.stripplot(x= 'gene_res', y = 'value', hue = 'variable', linewidth = 0.5, edgecolor = 'gray', hue_order = order, data = gene, jitter = 0.2, dodge=True, palette = ['b','g'], ax = axs)
    plt.ylim(-4.5,4.5)
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
    plt.savefig('{}/figures/fig2/{}_fig2hij_bar_{}.pdf'.format(dir, date, gene_name), bbox_inches='tight',transparent = True)
    plt.close()

    
final_fig2hij_boxplots('FLNB', 5) #flnb
final_fig2hij_boxplots('PIAS1', 0.8) #pias
final_fig2hij_boxplots('MAP2K4', 1.8) #map2k4

