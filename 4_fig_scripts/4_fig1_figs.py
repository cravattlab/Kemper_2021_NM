#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 08:06:35 2021

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


dir = os.getcwd()
date = '20210830'
x_med = 'cysteine_Mitosis/Asynch_combined_median'
y_med = 'protein_Mitosis/Asynch_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date), sheet_name = 'filtered')
prot_df = pd.read_excel('{}/supplemental_tables/{}_fig1_proteins.xlsx'.format(dir,date), sheet_name = 'all prots with changes')


def fig1a_waterfall_plot(df): 
    #source data is main dataset
    plt.figure(figsize = (3.3,3))
    df.reset_index()
    y = np.log2(df[x_med])
    x = np.arange(0,len(df))
    up = (y > np.log2(target_cutoff))
    down = (y < np.log2(1/target_cutoff))
    both = up|down
    plt.scatter(x[~both], y[~both], lw = 0.5, c = '#404040', marker = 'o', s = 9, alpha = 0.7)
    plt.scatter(x[up], y[up], lw = 0.5, c = 'blue', marker = 'o', s = 9, alpha = 0.7)
    plt.scatter(x[down], y[down], lw = 0.5, c = 'red', marker = 'o', s = 9, alpha = 0.7)
    plt.xlabel('Peptides', fontsize = 10, family = 'ARIAL')
    plt.ylabel('Mitosis/Asynch\nlog2 MS3 reporter ion intensity',fontsize = 10, family = 'ARIAL')
    plt.xticks(fontsize = 10, family = 'ARIAL')
    plt.yticks(ticks = [-6,-4,-2,0,2,4,6], fontsize = 10, family = 'ARIAL')
    plt.axhline(y = 1, c = 'black', linestyle = 'dotted')
    plt.axhline(y = -1, c = 'black', linestyle = 'dotted')
    plt.ylim(-5,5) 
    plt.title('Cysteine reactivity', fontsize = 10, family = 'ARIAL')
    print('Making waterfall plot')
    plt.savefig('{}/figures/fig1/{}_fig1a_waterfall_prot.pdf'.format(dir, date), bbox_inches='tight')
    plt.close()  
fig1a_waterfall_plot(df)

def fig1b_comparison_plot(df):
    #source data is main dataset
    print('making comparison plot')
    df2 = df.loc[~df[x_med].isnull()]
    df2 = df.loc[~df[y_med].isnull()]
    df2 = df2.sort_values(x_med, ascending = True)
    tc = np.log2(target_cutoff)
    plt.figure(figsize = (2.8,2.8))
    df.reset_index()
    y = np.log2(df[y_med])
    x = np.log2(df[x_med])
    down_cond = (x < -tc) & (y > np.log2(1/1.6))
    up_cond = (x > tc) & (y < np.log2(1.6))
    both = up_cond|down_cond
    plt.scatter(x[~both], y[~both], lw = 0, c = '#404040', marker = 'o',s = 8, alpha = 0.6)
    plt.scatter(x[up_cond], y[up_cond], lw = 0, c = 'blue', marker = 'o', s = 8, alpha = 0.6)
    plt.scatter(x[down_cond], y[down_cond], lw = 0, c = 'red', marker = 'o', s = 8, alpha = 0.6)
    plt.ylabel('Protein expression\nlog2 MS3 reporter ion intensity', fontsize = 10, family = 'ARIAL')
    plt.xlabel('Cysteine reactivity\nlog2 MS3 reporter ion intensity',fontsize = 10, family = 'ARIAL')
    plt.title('Mitosis/Asynch ratio',fontsize = 10, family = 'ARIAL')
    plt.xticks(ticks = [-4,-2,0,2,4],fontsize = 10, family = 'ARIAL')
    plt.yticks(ticks = [-4,-2,0,2,4], fontsize = 10, family = 'ARIAL')
    plt.axhline(y = np.log2(1.6), xmin = 0.6125, xmax = 1, c = 'black', linestyle = 'dotted')
    plt.axhline(y = -np.log2(1.6), xmin = 0, xmax = 0.3875, c = 'black', linestyle = 'dotted')
    plt.axvline(x = 1, ymin = 0, ymax = 0.5855, c = 'black', linestyle = 'dotted')
    plt.axvline(x = -1, ymin = 0.4155, ymax = 1, c = 'black', linestyle = 'dotted')
    plt.ylim(-4.5,4.5)
    plt.xlim(-4.5,4.5) 
    plt.savefig('{}/figures/fig1/{}_fig1b_correlation.pdf'.format(dir, date), bbox_inches='tight')
    plt.close()
fig1b_comparison_plot(df)


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
    
def fig1c_heatmap_plot(df):
    #source data is main dataset
    print('making heatmap')
    sns.set(style="white")
    df = df.loc[(~df[x_med].isnull()) &( ~df[y_med].isnull())]
    df = df.assign(nothing = 1)
    df2 = df[['gene_res',x_med,'nothing',y_med]]
    df2 = df2.dropna(axis = 0, how = 'any')
    df2.set_index('gene_res', inplace = True)
    df2 = np.log2(df2)
    df2.columns = ['Cysteine\nreactivity', '','Protein\nexpression']
    df3 = df2.applymap(white)
    df3 = df3.sort_values('Protein\nexpression', ascending = False)
    print(df3)
    df3 = df3.reset_index().drop('gene_res', axis = 1)
    plt.figure(figsize = (1.8,5))
    cmap = ['#990000','#ff3333','#ffcccc','white', '#e6e6ff', '#4d4dff', '#000099']
    ax = sns.heatmap(df3, yticklabels = 200,cmap = cmap)#, cbar_kws= {'label': 'log2 Mitosis/Asynch MS3 reporter ion intensity'})
    colorbar = ax.collections[0].colorbar
    colorbar.remove()
    ax.figure.axes[-1].yaxis.label.set_size(10)
    ax.figure.axes[-1].tick_params(labelsize=10)
    ax.set_ylabel("Peptides", fontsize = 10, family = 'ARIAL')
    plt.xticks(fontsize = 10, family = 'ARIAL')
    plt.yticks(np.arange(0,15000,1000).tolist(),np.arange(0,15000,1000).tolist() , fontsize = 10, family = 'ARIAL')
    plt.title('Mitosis/Asynch', family = 'ARIAL', fontsize = 10)
    ax.vlines([0.99,2], ymin = 0, ymax = 14297, color = 'black', lw = 1)
    plt.savefig('{}/figures/fig1/{}_fig1c_basal_whole-proteome_heatmap.pdf'.format(dir,date), transparent = True, bbox_inches='tight')
    plt.close()
fig1c_heatmap_plot(df)


def fig1d_pie_source_data(df):
    #get nan in the index from excel
    df1 = df.copy()
    df1[['accession', 'protein']] = df[['accession','protein']].fillna(method = 'ffill')
    # df1['react'] = df1['react'].replace('Unassigned', np.nan)
    df1['react'] = df1['gene_res'] + ': ' + df1['react']
    df1['exp'] = df1['gene_res'] + ': ' + df1['exp']
    df1['cys_category'] = df1[['react','exp']].apply(lambda x: x.str.cat(sep=''), axis=1)
    df2 = df1[['accession','protein','cys_category']]
    # df3 = df2.loc[~df2.cys_category.str.contains('Unassigned')]
    df3 = df2.loc[~(df2['cys_category'] == '')]
    df4 = pd.DataFrame(df3.groupby(['accession', 'protein'])['cys_category'].sum()).reset_index()
    df4.loc[df4.cys_category.str.contains('Unassigned', regex = True),'Reactivity or Expression?'] = 'Unassigned'
    df4.loc[df4.cys_category.str.contains('EXPRESSION'),'Reactivity or Expression?'] = 'Expression'
    df4.loc[df4.cys_category.str.contains('REACTIVITY'),'Reactivity or Expression?'] = 'Reactivity'
    df4.loc[((df4.cys_category.str.contains('REACTIVITY')) & \
             (df4.cys_category.str.contains('EXPRESSION'))),'Reactivity or Expression?'] = 'Both'
    df4.drop('cys_category', axis = 1,inplace = True)
    #save for source data
    df4.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig1d.csv'.format(dir), index = False, float_format = '%.2f')
    return df4

def fig1d_pie(df):
    #for venn diagram, if a protein has both expression and reactivity changes we will count them for both
    df4 = fig1d_pie_source_data(df)
    reactivity = df4['Reactivity or Expression?'].str.contains('Reactivity')
    expression = df4['Reactivity or Expression?'].str.contains('Expression')
    both = df4['Reactivity or Expression?'].str.contains('Both')
    react = df4.loc[(reactivity) | (both)]
    exp = df4.loc[(expression) | (both)]
    labels = 'Reactivity changes','Expression changes'
    sizes = [len(set(react.accession)),len(set(exp.accession)) ]
    total = sum(sizes)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.4,2.4)
    ax1.pie(sizes, autopct = (lambda p: '{:.0f}'.format(p * total / 100)) , colors = \
        [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),\
            (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),\
                (0.8666666666666667, 0.5176470588235295, 0.3215686274509804)], startangle=90)
    ax1.legend(labels = labels, loc='center left', bbox_to_anchor=(0.02, -0.3),frameon=False, handletextpad=0.3, fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/fig1/{}_fig1d_pie.pdf'.format(dir, date), bbox_inches='tight',transparent = True)
    plt.close()
fig1d_pie(prot_df)

def set_hue(row):
    if row['variable'] == x_med:
        if row['value'] >= 0.9:
            return 'Cysteine reactivity lower\nin mitosis'
        elif row['value'] <= -0.9:
            return 'Cysteine reactivity higher\nin mitosis'
        else:
            return 'Unchanging'
    elif row['variable'] == y_med:
        return 'Protein expression'

def fig1e_cat_scatter(df, input, sub, colors):
    sns.set(style="ticks")
    print('making final scatter plot figures')
    gene = df.loc[df.protein.isin(input)]
    median_cols = [col for col in gene.columns if 'median' in col]
    cat_cols = [col for col in gene.columns if 'category' in col]
    sdf = gene.copy()
    sdf[median_cols] = np.log2(sdf[median_cols])
    sdf.drop(cat_cols, axis =1 , inplace = True)
    sdf = sdf.sort_values('protein', ascending = True)
    sdf.set_index(['accession', 'description', 'protein','cysteine_sequence','accession_res','gene_res'], inplace = True)
    sdf.to_excel('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig1e_{}.xlsx'.format(dir, sub),float_format = '%.2f')
    gene = gene.sort_values('gene_res', ascending = True)
    gene = gene[['protein','gene_res',x_med,y_med]]
    gene = gene.melt(id_vars = ['protein','gene_res'])
    gene = gene.loc[~gene.value.isnull()]
    gene = gene.drop_duplicates(subset = ['protein','variable','value'])
    gene['value'] = np.log2(gene['value'])
    gene['hue_cat'] = gene.apply(set_hue, axis = 1)
    markers = {x_med: "o", y_med: "^"}
    g= sns.scatterplot(x="protein", y='value',style = 'variable',  alpha = 0.9, s = 50, legend = False, linewidth=0.5,  hue = 'hue_cat' ,data=gene, markers = markers, palette = colors) 
    fig = plt.gcf()
    fig.set_size_inches(2.8,2.8)
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.ylabel('Mitosis/Asynch\nlog2 MS3 reporter ion intensity',fontsize = 10, family = 'ARIAL')
    plt.xlabel('')
    plt.title(sub, fontsize = 10, family = 'ARIAL')
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    # plt.legend(title = False, loc='center left', bbox_to_anchor=(1, 0.5),frameon=False, handletextpad=0)
    plt.ylim(-4,4)
    plt.savefig('{}/figures/fig1/{}_fig1e_scatter_{}.pdf'.format(dir,date, sub), bbox_inches='tight',transparent = True)
    plt.close()

# def fig1e_cat_scatter(df, input, sub, colors):
#     sns.set(style="ticks")
#     print('making final scatter plot figures')
#     gene = df.loc[df.protein.isin(input)]
#     gene.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig1e_{}.csv'.format(dir, sub), index = False)
#     gene = gene.sort_values('gene_res', ascending = True)
#     gene = gene[['protein','gene_res',x_med,y_med]]
#     gene = gene.melt(id_vars = ['protein','gene_res'])
#     gene = gene.loc[~gene.value.isnull()]
#     gene = gene.drop_duplicates(subset = ['protein','variable','value'])
#     gene['value'] = np.log2(gene['value'])
#     gene['hue_cat'] = gene.apply(set_hue, axis = 1)
#     markers = {x_med: "o", y_med: "^"}
    
#     def jitter(row):
#         x = row['type_id']
#         prot = row['hue_cat']
#         if prot == 'Protein expression':
#             new_x = x
#         else:
#             new_x = x + random.uniform(0, .5) -.25
#         return new_x
    
#     set_np = list(range(len(input)))
#     type_ids = dict(zip(input, set_np))
#     gene['type_id'] = gene['protein'].apply(lambda x: type_ids[x])
#     x =  gene.apply(jitter, axis =1)
#     print(x)
#     g= sns.scatterplot(x=x, y='value',style = 'variable',  alpha = 0.9, s = 50, legend = False, linewidth=0.5,  hue = 'hue_cat' ,data=gene, markers = markers, palette = colors) 
#     fig = plt.gcf()
#     fig.set_size_inches(2.8,2.8)
#     plt.axhline(y=1, color = 'black', linestyle = 'dotted')
#     plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
#     plt.ylabel('Mitosis/Asynch\nlog2 MS3 reporter ion intensity',fontsize = 10, family = 'ARIAL')
#     plt.xlabel('')
#     plt.title(sub, fontsize = 10, family = 'ARIAL')
#     plt.xticks(ticks = set_np, labels = input, rotation = 40, ha = 'right', fontsize = 10, family = 'ARIAL')
#     # plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
#     plt.yticks(fontsize = 10, family = 'ARIAL')
#     # plt.legend(title = False, loc='center left', bbox_to_anchor=(1, 0.5),frameon=False, handletextpad=0)
#     plt.ylim(-4,4)
#     plt.savefig('{}/figures/fig1/{}_fig1e_scatter_{}.pdf'.format(dir,date, sub), bbox_inches='tight',transparent = True)
#     plt.close()
    
def run_fig1e_scatter(df):
    #expression
    sub = 'Expression'
    input = ['AURKA', 'CCNB1','CDKN2A','MCM6','UBE2S','CCNA2']
    colors = ['blue', 'red', 'green']
    # input = ['O14965','P14635','Q16763','Q8N726','Q14566','P30291']
    fig1e_cat_scatter(df,input, sub, colors)
    #reactivity
    input = ['ACIN1','FLNB', 'HAUS6', 'PPP2R2D','STAT3','TRIM33'] #
    sub = 'Reactivity'
    colors = ['#404040', 'red', 'blue','green']
    # input = ['O75369','Q9Y5M8','Q9UPN9','Q9NQW6','Q9Y570','O75874']
    fig1e_cat_scatter(df,input, sub, colors)   
    
run_fig1e_scatter(df)

def fig1f_boxplots(protein_name, width):
    gene1 = df.loc[df.protein == protein_name]
    median_cols = [col for col in gene1.columns if 'median' in col]
    cat_cols = [col for col in gene1.columns if 'category' in col]
    sdf = gene1.copy()
    sdf[median_cols] = np.log2(sdf[median_cols])
    sdf.drop(cat_cols, axis =1 , inplace = True)
    sdf = sdf.sort_values('protein', ascending = True)
    sdf.set_index(['accession', 'description', 'protein','cysteine_sequence','accession_res','gene_res'], inplace = True)
    sdf.to_excel('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig1f_suppfig1f_{}.xlsx'.format(dir, protein_name), float_format = '%.2f')
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
    sns.stripplot(x= 'gene_res', y = 'value', palette = hue_order, order = order, \
                  linewidth = 0.5, edgecolor= 'black', data = gene, jitter = True, dodge=True,  ax = axs)
    plt.ylim(-4.5,4.5)
    plt.title(gene_name, fontsize = 10)
    handles, labels = axs.get_legend_handles_labels() 
    plt.axhline(y=1, color = 'black', linestyle = 'dotted')
    plt.axhline(y=-1, color = 'black', linestyle = 'dotted')
    plt.ylabel('Mitosis/Asynch cysteine reactivity\nlog2 MS3 reporter ion intensity',family = 'ARIAL', fontsize = 10)
    plt.xlabel('')
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 40, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.title(gene_name, family = 'ARIAL', fontsize = 10)
    plt.savefig('{}/figures/fig1/{}_fig1f_bar_{}.pdf'.format(dir, date, gene_name), bbox_inches='tight',transparent = True)
    plt.close()

fig1f_boxplots('RRP15', 0.5) #sext
fig1f_boxplots('NOL8', 2) #arhgap39 P46013
fig1f_boxplots('DTWD1', 1.2) #dtwd


