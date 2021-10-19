#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 15:23:44 2021

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
import re

plt.rcParams['pdf.fonttype'] = 42
font = {'family' : 'ARIAL',
        'weight' : 'normal',
        'size'   : 10}         
plt.rc('font', **font) #set the font style created
plt.rcParams.update({'text.color' : "black"})

lpp_cutoff = 2
asynch_cutoff = 1.6
unchanging_cutoff = 1.6
adapted_cutoff = 1.8

protein_cols = ['accession','accession_res','description', 'protein', 'sequence','gene_res']
phospho_col = ['asynch_stoichiometry_mann','mitosis_stoichiometry_mann', 'asynch_stoichiometry_olsen','mitosis_stoichiometry_olsen']
mann = 'mitosis_stoichiometry_mann'

dir = os.getcwd()
date = '20210830'

x_med ='Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
adapted_med ='Mitosis_LPP(-,+)/(+,+)_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
prots_df = pd.read_excel('{}/supplemental_tables/{}_fig2_proteins.xlsx'.format(dir,date), sheet_name = 'all prots with changes')
y_df = pd.read_excel('{}/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date), sheet_name = 'filtered')
y_exp = pd.read_excel('{}/supplemental_tables/{}_fig1_proteins.xlsx'.format(dir,date), sheet_name = 'all prots with changes')
b_is = pd.read_excel('{}/other_datasets/backus.xlsx'.format(dir), sheet_name = 'Probe Targets in situ')
b_iv = pd.read_excel('{}/other_datasets/backus.xlsx'.format(dir), sheet_name = 'Probe Targets In vitro')
bpk = pd.read_csv('{}/other_datasets/kemper.csv'.format(dir))
vino = pd.read_csv('{}/other_datasets/vinogradova.csv'.format(dir))
lpp_df = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')

def fig4a_source_data():
    #for source table
    df1 = df[protein_cols + [mann, y_med, x_med, adapted_med]]
    lpp_sensitive = lpp_df.loc[lpp_df['cysteine_'+x_med] >= 2]
    total = []
    for id in list(set(df1.accession)):
        id_df = df1.loc[df1.accession == id].reset_index(drop = True)
        am_diff = id_df.loc[(id_df[y_med] >= asynch_cutoff) | (id_df[y_med] <= 1/asynch_cutoff)]
        changing = am_diff.loc[(id_df[x_med] <= 1/lpp_cutoff ) | (id_df[x_med] >= lpp_cutoff) ]
        down = am_diff.loc[id_df[x_med] <= 1/lpp_cutoff]
        authentic = down.loc[down[adapted_med] <= 1/adapted_cutoff]
        artifact = down.loc[down[adapted_med] > 1/unchanging_cutoff]  
        real_changing = changing.loc[~changing.accession_res.isin(list(set(artifact.accession_res)))]
        lpp = lpp_sensitive.loc[lpp_sensitive.accession.str.contains(id)]
        if real_changing.empty:
            change = False
            auth = False
            art = False
        else:
            change = 'TRUE, ' + str(real_changing.gene_res.to_list())
        if authentic.empty:
            auth = False
        else: 
            auth = 'TRUE, ' + str(authentic.gene_res.to_list())
        if artifact.empty:
            art = False
        else: 
            art = 'TRUE, ' + str(artifact.gene_res.to_list())
        if id_df[mann].dropna().empty:
            phos = False
        else:
            phos = 'True, ' + id_df[mann][0]
        if lpp.empty:
            lpp_id = False
        else:
            lpp_id = 'True, ' + str(lpp.gene_res.to_list())
        prot_info = [id, id_df.protein[0], id_df.description[0], phos, lpp_id, change, art, auth]
        total.append(prot_info)
    source_data = pd.DataFrame(total, columns = ['accession', 'protein','description',\
                                               'High stoichiometry phosphorylation? (Sharma et al)',\
                                                   'LPP-sensitive phosphorylation? (This paper)',\
                                          'Phosphorylation-dependent cysteines (not including artifactual)?','Artifactual decrease in LPP(-)?', \
                                              'Authentic decrease in LPP(-)?' ])
    source_data.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig4a_suppfig4ab_pie.csv'.format(dir), index = False)
fig4a_source_data()    
    
def fig4a_data_generation():
    #proteins with high stoich phosphorylation, fig5
    no_ppep = ~(df[adapted_med] > 1/unchanging_cutoff)
    down_cond = (df[x_med] <= 1/lpp_cutoff) &  (no_ppep)
    lpp_cond = (df[x_med] >= lpp_cutoff) | (down_cond)
    cond2 = (df[y_med] <= (1/asynch_cutoff)) |  (df[y_med] >= asynch_cutoff)
    #peptides that are changing with LPP and asynch but not artifactual changes
    am_diff = df.loc[(lpp_cond) & (cond2)]
    #no phospho
    mann = 'mitosis_stoichiometry_mann'
    #proteins that have quantified cysteines in mann or olsen in LPP changing proteins (not artifactual))
    no_mann = am_diff.loc[am_diff[mann].isnull()]
    overlap = am_diff.loc[~ am_diff[mann].isnull()]
    #proteins that have quantified cysteines in mann or olsen in all proteins
    all_no_mann = df.loc[df[mann].isnull()]
    all_overlap = df.loc[~ df[mann].isnull()]
    #number of proteins in mann, number not in mann
    df1 = df.loc[(df[y_med] >= asynch_cutoff) | (df[y_med] <= 1/asynch_cutoff)]
    down = df1.loc[df1[x_med] <= 1/lpp_cutoff]
    authentic = down.loc[down[adapted_med] <= 1/adapted_cutoff]
    artifact = down.loc[down[adapted_med] > 1/unchanging_cutoff]    
    auth_no = authentic.loc[authentic[mann].isnull()]
    auth_overlap = authentic.loc[~ authentic[mann].isnull()]
    art_no = artifact.loc[artifact[mann].isnull()]
    art_overlap = artifact.loc[~ artifact[mann].isnull()]
    
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
    print('Writing Fig4 pie chart to csv')
    return fig4_pies

def fig4a_pie(titles, sizes, sub):
    labels = 'Proteins with high stoichiometry\nphosphorylation in mitosis','Low or unquantified stoichiometry'    
    total = sum(sizes)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.1,2)
    plt.title(titles, fontsize = 10)
    ax1.pie(sizes, autopct = (lambda p: '{:.0f}'.format(p * total / 100)) , startangle=90)
    ax1.legend(labels = labels, loc='center left', bbox_to_anchor=(-0.2, -0.15),frameon=False, handletextpad=0.3, fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/fig4/{}_fig4a_pie_{}.pdf'.format(dir, date, sub), bbox_inches='tight',transparent = True)
    plt.close()

def fig4a_generation():
    fig4_pies = fig4a_data_generation()
    sizes = fig4_pies['all']
    titles = 'All quantified proteins'
    sub = 'all_proteins'
    fig4a_pie(titles, sizes, sub)
    titles = 'Proteins with phosphorylation-dependent\ncysteine reactivity changes'
    sizes = fig4_pies['lpp_change']
    sub = 'lpp_proteins'
    fig4a_pie(titles, sizes, sub)
fig4a_generation()
    
def fig4b_kegg_data():
    #this will be used to help color the map for fig 4b
    kegg = pd.read_csv('{}/other_datasets/KEGG_04110_cellcycle.csv'.format(dir))
    protein_cols = ['accession',  'protein',
       'accession_res', 'gene_res']
    y_df1 = y_df[protein_cols + [y_med]]
    y_df1 = y_df1.assign(asynch_change = 'FALSE')
    y_df1.loc[(y_df1[y_med] >= 2) | (y_df1[y_med] <= 0.5),'asynch_change'] = 'TRUE'
    x_df1 = df[protein_cols + [x_med]]
    x_df1 = x_df1.assign(lpp_change = 'FALSE')
    lpp_cond = (df[x_med] >= lpp_cutoff) | (df[x_med] <= 1/lpp_cutoff)
    am_cond = (df[y_med] >= asynch_cutoff) | (df[y_med] <= 1/asynch_cutoff)
    lpp_all = df.loc[(lpp_cond) ]
    lpp_changing = df.loc[(lpp_cond) & (am_cond)]
    x_df1.loc[x_df1.accession_res.isin(list(set(lpp_all.accession_res))),'lpp_change'] = 'TRUE'
    x_df1.loc[x_df1.accession_res.isin(list(set(lpp_changing.accession_res))),'lpp_change'] = 'TRUE1'
    df2 = y_df1.merge(x_df1, how = 'outer', on = protein_cols)
    df2 = df2.assign(key_cat = 'QUANTIFIED')
    df2.loc[(df2.lpp_change == 'TRUE'), 'key_cat'] = 'LPP_CHANGE'
    df2.loc[(df2.asynch_change == 'TRUE'), 'key_cat'] = 'ASYNCH_CHANGE'
    df2.loc[(df2.lpp_change == 'TRUE1'), 'key_cat'] = 'BOTH'
    df2.drop(['lpp_change','asynch_change'], axis =1, inplace = True)
    df2 = df2.assign(gene = df2.gene_res.str.split('_', expand = True)[0])
    df3 = kegg.merge(df2, how = 'left', on = 'gene')
    df3['gene_cat'] = df3.gene_res + ': ' + df3.key_cat
    df4 = pd.DataFrame(df3.groupby(['gene'])['gene_cat'].apply(lambda row: ', '.join(row.astype(str)))).reset_index()
    df4['accession'] = df3.groupby(['gene'])['accession'].first().to_list()
    df4['kegg_id'] = df3.groupby(['gene'])['kegg_id'].first().to_list()
    df4['KEGG MAP'] = 'NOT_QUANTIFIED'
    df4.loc[(df4['gene_cat'].str.contains('QUANTIFIED')), 'KEGG MAP'] = 'QUANTIFIED'
    df4.loc[(df4['gene_cat'].str.contains('LPP_CHANGE')), 'KEGG MAP'] = 'LPP(-)/LPP(+)'
    df4.loc[(df4['gene_cat'].str.contains('ASYNCH_CHANGE')), 'KEGG MAP'] = 'MITOSIS/ASYNCH'
    df4.loc[(df4['gene_cat'].str.contains('BOTH')), 'KEGG MAP'] = 'BOTH'
    print(df4)
    df4.to_excel('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig4b_KEGG.xlsx'.format(dir), index = False)
fig4b_kegg_data()

def fig4c_GO_lists():
    background = list(set(df.loc[~df[y_med].isnull(), 'accession']))
    #changing in LPP/CTRL, changing in Mitosis/Asynch
    cond1 = (df[x_med] <= (1/lpp_cutoff)) |  (df[x_med] >= lpp_cutoff)
    cond2 = (df[y_med] <= (1/asynch_cutoff)) |  (df[y_med] >= asynch_cutoff)
    #changing in basal only
    am_diff = df.loc[(cond1) & (cond2)]
    no_ppep = am_diff.loc[~( (am_diff[x_med] <= 1/lpp_cutoff) & (am_diff[adapted_med] > 1/unchanging_cutoff))]
    background = pd.DataFrame(background)
    change_no_ppep = pd.DataFrame({'all': list(set(no_ppep.accession))})
    change_no_ppep.to_csv('{}/for_figures/fig4/GO_analysis/inputs/{}_all_changing_noppep.txt'.format(dir, date), index = False, header = False)
    background.to_csv('{}/for_figures/fig4/GO_analysis/inputs/{}_LPP_background.txt'.format(dir, date), index = False, header = False)
fig4c_GO_lists()

def fig4c_GO_chart(sub, type):
    #input data generated in 3_fig1_proteins
    #then put through webgestalt, cutoff 0.05
    #then put through revigo to reduce redundancy - this is source data
    df = pd.read_csv('{}/for_figures/fig4/GO_analysis/fig4_revigo_{}_{}_simrei_05.csv'.format(dir, sub, type))
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
    plt.savefig('{}/figures/fig4/{}_fig4c_{}_{}.pdf'.format(dir, date, sub, type), bbox_inches='tight')
    # plt.close()
fig4c_GO_chart('noppep','cellular')

def fig4d_disorder_parse(df):
    #to merge all the IUPRED results together
    df.columns = ['foo']
    df = df.foo.str.split(pat = '|', n = 1, expand = True)
    df.columns = ['data','id']
    df['id'] = df['id'].ffill()
    headers = df.loc[3, 'data'].split(sep = '\t')
    df = df.loc[~df.data.str.contains('>')]
    df = df.loc[~df.data.str.contains('#')]
    df[headers] = df.data.str.split(pat = '\t', expand = True)
    df[['accession','dele']] = df.id.str.split(pat = '|', n = 1, expand = True)
    # df[['protein']] = df.dele.str.split(pat = ' ', n =2, expand = True)[1]
    df[['accession_res']] = df['accession'] + '_' + df['# POS'].astype(str)
    df = df.drop(['dele', 'id','data'], axis = 1)
    df['IUPRED SCORE'] = df['IUPRED SCORE'].astype(float)
    df['ANCHOR SCORE'] = df['ANCHOR SCORE'].astype(float)
    return df

def fig4d_disorder_map():
    ##IUPRED2A tool only can handle about 1500-2000 entries at once
    #I think I manually just chose ~1500-2000 entries from the entire fasta database and listed them as 20210316_all_fasta_{}_.fasta (a-k)
    #really only need to do this once 
    alphabet = list(map(chr, range(97, 108)))
    dis_map = pd.DataFrame()
    for input in alphabet:
        print(input)
        df = pd.read_csv('{}/other_datasets/IUPRED/iupred_results/20210303_all_{}_fasta.result'.format(dir, input), index_col = False)
        df1 = fig4d_disorder_parse(df)
        dis_map = dis_map.append(df1)
    dis_map.to_csv('{}/other_datasets/IUPRED/20210303_all_fasta.csv'.format(dir), index = False)
    
def fig4e_lig_map():
    #only need to do this initially
    bpk_drop = [col for col in bpk.columns if 'stdev' in col]
    bpk_drop = bpk_drop + [col for col in bpk.columns if 'Number of detections' in col]
    bpk1 = bpk.drop(bpk_drop+['description','symbol','sequence'], axis = 1)
    vino1 = vino.drop(['accession_sequence','description'], axis =1)
    b_is1 = b_is.drop('protein_description_sequence1_sequence', axis = 1)
    b_iv1 = b_iv.drop('protein_description_sequence1_sequence', axis = 1)
    bpk1 = bpk1.assign(res = bpk1.res.str.split(',')).explode('res')
    bpk1 = bpk1.assign(accession_res = (bpk1.accession + '_' + bpk1.res)).drop(['res','accession'], axis = 1)
    vino1 = vino1.assign(res = vino1.res.str.split(',')).explode('res')
    vino1 = vino1.assign(accession_res = (vino1.accession + '_' + vino1.res)).drop(['res','accession'], axis = 1)
    df_list = [bpk1, vino1, b_is1, b_iv1]
    df_list = [df.replace('--', np.nan) for df in df_list]
    df_list = [df[~df.accession_res.str.contains("\+")] for df in df_list]
    bpk1, vino1, b_is1, b_iv1 = df_list
    bpk1_cols = [col for col in bpk1.columns if '500uM' in col]
    vino1_cols = [col for col in vino1.columns if 'accession_res' not in col]
    b_is1_cols = [col for col in b_is1.columns if 'insitu' in col]
    b_iv1_cols = [col for col in b_iv1.columns if 'invitro' in col]
    bpk1[bpk1_cols] = bpk1[bpk1_cols].astype(float)
    vino1[vino1_cols] = vino1[vino1_cols].astype(float)
    b_is1[b_is1_cols] = b_is1[b_is1_cols].astype(float)
    b_iv1[b_iv1_cols] = b_iv1[b_iv1_cols].astype(float)
    bpk1 = bpk1.assign(bpk1_bool = (bpk1[bpk1_cols] >= 4).any(axis = 1))
    vino1 = vino1.assign(vino1_bool = (vino1[vino1_cols] >= 4).any(axis = 1))
    b_is1 = b_is1.assign(b_is1_bool = (b_is1[b_is1_cols] >= 4).any(axis = 1))
    b_iv1 = b_iv1.assign(b_iv1_bool = (b_iv1[b_iv1_cols] >= 4).any(axis = 1))
    ligand = bpk1.merge(vino1, on = 'accession_res', how = 'outer')
    ligand = ligand.merge(b_is1, on = 'accession_res', how = 'outer')
    ligand = ligand.merge(b_iv1, on = 'accession_res', how = 'outer')
    bool_cols = [col for col in ligand.columns if '_bool' in col]
    lig_map = ligand[['accession_res'] + bool_cols].drop_duplicates()
    lig_map = lig_map.assign(liganded_in_any = lig_map[bool_cols].any(axis = 1))
    lig_map = lig_map.sort_values('liganded_in_any', ascending = False)
    foo = lig_map.loc[lig_map.accession_res.duplicated(keep = 'first')].sort_values('accession_res')
    bar = lig_map.loc[lig_map.accession_res.duplicated(keep = 'last')].sort_values('accession_res')
    lig_map = lig_map.loc[~lig_map.index.isin(foo.index)]
    lig_map.to_csv('{}/other_datasets/ligandability_map.csv'.format(dir), index = False)
    return lig_map


def fig4d_fig4e_prep():
    # dis_map = fig4d_disorder_map()
    # lig_map = fig4d_lig_map()
    lig_map = pd.read_csv('{}/other_datasets/ligandability_map.csv'.format(dir))
    dis_map = pd.read_csv('{}/other_datasets/IUPRED/20210303_all_fasta.csv'.format(dir))
    dis_map = dis_map.assign(accession_res = (dis_map.accession + '_C' + dis_map.accession_res.str.rsplit('_',expand = True)[1]))
    dis_map = dis_map[['accession_res', 'IUPRED SCORE']]
    df1 = df.loc[~df[y_med].isnull()]
    lpp_cond = (df1[x_med] >= lpp_cutoff) | (df1[x_med] <= 1/lpp_cutoff)
    am_cond = (df1[y_med] >= asynch_cutoff) | (df1[y_med] <= 1/asynch_cutoff)
    am_diff = df1.loc[(lpp_cond) & (am_cond)]
    drop_cols = [col for col in am_diff.columns if 'median_TMT' in col]#+ [col for col in am_diff.columns if '_CV' in col]+ [col for col in am_diff.columns if '_no' in col]
    df1 = df1.drop(drop_cols, axis = 1)
    df1 = df1.assign(accession_res = (df1.accession + '_C' + df1.gene_res.str.split('_',expand = True)[1]))
    am_diff = am_diff.drop(drop_cols, axis = 1)
    am_diff = am_diff.assign(accession_res = (am_diff.accession + '_C' + am_diff.gene_res.str.split('_',expand = True)[1])) 
    bool_cols = [col for col in lig_map.columns if '_bool' in col]
    am_diff1 = am_diff.merge(lig_map, how = 'left', on = 'accession_res')
    am_diff1 = am_diff1.assign(liganded_in_any = am_diff1[bool_cols].any(axis =1)).drop_duplicates()
    am_diff2 = am_diff1.merge(dis_map, how = 'left', on = 'accession_res')
    df2 = df1.merge(lig_map, how = 'left', on = 'accession_res')
    df2 = df2.assign(liganded_in_any = df2[bool_cols].any(axis = 1)).drop_duplicates()
    df3 = df2.merge(dis_map, how = 'left', on = 'accession_res')
    down =  am_diff2.loc[(am_diff2[x_med] <= 1/lpp_cutoff)]
    down =  down[~(down[adapted_med] > 1/unchanging_cutoff)]
    up = am_diff2.loc[(am_diff2[x_med] >= lpp_cutoff)]
    # df3.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig4de.csv'.format(dir), index = False)
    return up, down,  df3

def fig4d_disorder_prep():
    up, down, df3 = fig4d_fig4e_prep()
    labels = ['Cysteines from\nall proteins','Increase', 'Decrease']
    dis_number = [len(set(df3.loc[df3['IUPRED SCORE'] >= 0.5, 'accession_res'])),\
                 len(set(up.loc[up['IUPRED SCORE'] >= 0.5, 'accession_res'])),\
                 len(set(down.loc[down['IUPRED SCORE'] >= 0.5, 'accession_res']))]
    total_number = [len(set(df3.accession_res)),len(set(up.accession_res)),len(set(down.accession_res))]
    print(labels, dis_number, total_number)
    return labels, dis_number, total_number
    
    
def fig4e_lig_prep():
    up, down, df3 = fig4d_fig4e_prep()
    a_l = df3.loc[df3['liganded_in_any']]
    down_chl = down.loc[down['liganded_in_any']]
    up_chl = up.loc[(up[x_med] >= lpp_cutoff) & (up['liganded_in_any']) ] 
    labels = ['Cysteines from\nall proteins','Increase','Decrease']
    dis_number = [len(set(a_l.accession_res)),len(set(up_chl.accession_res)),len(set(down_chl.accession_res)) ]
    total_number = [len(set(df3.accession_res)), len(set(up.accession_res)),len(set(down.accession_res))]
    return labels, dis_number, total_number

def fig4de_source_data():
    df = fig4d_fig4e_prep()[2]
    keep_cols = [col for col in df.columns if 'combined_median' in col]
    keep_cols = [col for col in keep_cols if 'protein' not in col]
    am_cond = (df[y_med] >= asynch_cutoff) | (df[y_med] <= 1/asynch_cutoff)
    up = df.loc[(df[x_med] >= lpp_cutoff) & (am_cond)]
    down = df.loc[(df[x_med] <= 1/lpp_cutoff) & (am_cond)]
    authentic = down.loc[down[adapted_med] <= 1/adapted_cutoff]
    artifact = down.loc[down[adapted_med] > 1/unchanging_cutoff]
    unassigned = down.loc[(down[adapted_med] <= 1/unchanging_cutoff) & (down[adapted_med] > 1/adapted_cutoff)]
    df1 = df[protein_cols + keep_cols + ['IUPRED SCORE', 'liganded_in_any']]
    df1['Cysteine reactivity in LPP(-)'] =  'Unchanged'
    df1.loc[df1.accession_res.isin(set(down.accession_res)), 'Cysteine reactivity in LPP(-)'] = 'Decreased'
    df1.loc[df1.accession_res.isin(list(set(up.accession_res))), 'Cysteine reactivity in LPP(-)'] = 'Increased'
    df1.loc[df1.accession_res.isin(list(set(authentic.accession_res))),\
            'Authentic or artifactual cysteine reactivity decrease in LPP(-)'] = 'Authentic'
    df1.loc[df1.accession_res.isin(list(set(artifact.accession_res))),\
            'Authentic or artifactual cysteine reactivity decrease in LPP(-)'] = 'Artifact'
    df1.loc[df1.accession_res.isin(list(set(unassigned.accession_res))),\
            'Authentic or artifactual cysteine reactivity decrease in LPP(-)'] = 'Unassigned'
    df1['Liganded in previous studies?'] =  df1['liganded_in_any'].copy()
    df1['Reside in disordered region IUPRED > 0.5?'] = False
    df1.loc[(df1['IUPRED SCORE']> 0.5), 'Reside in disordered region IUPRED > 0.5?'] = True
    df1.drop('liganded_in_any', axis =1 ,inplace = True)
    df1.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig4de.csv'.format(dir), index = False, float_format = '%.2f')
fig4de_source_data()

def fig4d_plot_disorder():   
    labels, dis_number, total_number = fig4d_disorder_prep()
    y_value =  [i / j for i, j in zip(dis_number, total_number)] 
    percent_labels = [str(i)+ '/\n' + str(j) for i,j in zip(dis_number, total_number)]    
    plt.subplots(figsize=(2,2.2))
    ax = sns.barplot(x = labels, y = y_value, palette = ('gray', 'blue', 'red'))
    plt.ylabel('Percent disordered', fontsize =10, family = 'ARIAL')
    plt.xticks(fontsize =10, family = 'ARIAL', rotation = 45, ha = 'right')
    plt.yticks(fontsize =10, family = 'ARIAL')
    plt.title('Predicted intrinsic disorder', fontsize =10, family = 'ARIAL')
    plt.ylim(0,1)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    for i in range(len(percent_labels)):
        plt.annotate(str(percent_labels[i]), xy=(i,y_value[i]), ha='center', va='bottom', fontsize = 10, family = 'ARIAL')
    plt.savefig('{}/figures/fig4/{}_fig4d_disorder.pdf'.format(dir, date), bbox_inches='tight')
fig4d_plot_disorder()


def fig4e_plot_ligandability():
    labels, dis_number, total_number = fig4e_lig_prep()
    y_value =  [i / j for i, j in zip(dis_number, total_number)] 
    percent_labels = [str(i)+ '/\n' + str(j) for i,j in zip(dis_number, total_number)]    
    plt.subplots(figsize=(2,2.2))
    ax = sns.barplot(x = labels, y = y_value, palette = ('gray', 'blue', 'red'))
    plt.ylabel('Percent ligandable', fontsize =10, family = 'ARIAL')
    plt.xticks(fontsize =10, family = 'ARIAL', rotation = 45, ha = 'right')
    plt.yticks(fontsize =10, family = 'ARIAL')
    plt.title('Ligandability', fontsize =10, family = 'ARIAL')
    plt.ylim(0,1)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    for i in range(len(percent_labels)):
        plt.annotate(str(percent_labels[i]), xy=(i,y_value[i]), ha='center', va='bottom', fontsize = 10, family = 'ARIAL')
    plt.savefig('{}/figures/fig4/{}_fig4d_ligandability.pdf'.format(dir, date), bbox_inches='tight')
fig4e_plot_ligandability()


def fig4f_moving_average(row):
    pos = int(row['pos']) - 1
    seq = row['sequence_fasta']       
    sp_loc = [abs(m.start() - pos) for m in re.finditer('SP|TP', seq)] 
    if not sp_loc:
        # print('No S/T-P sites')
        # count0 = count0 + 1
        total_counts = [0] * 10 #+ [count0]
        dist = 'No S/T-P sites'
    else:
        lowest_loc = [min(sp_loc)]
        # print('Closest SP or TP occurs for ' + str(row['gene_res']) + ' ' + str(lowest_loc) + ' away from Cys')
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
        dist = lowest_loc[0]
    return total_counts, dist

def fig4f_stp_prep(df): 
    fasta = pd.read_csv('{}/inputs/RS_2016_fasta_20210127.csv'.format(dir))
    fasta1 = fasta[['accession', 'sequence']]
    df1 = df.merge(fasta1, on = 'accession', how = 'left', suffixes = ['', '_fasta'])
    df1 = df1.assign(pos = df1.gene_res.str.split('_', n = 1, expand = True)[1])
    results = df1.apply(fig4f_moving_average, axis =1).to_list()
    counts = pd.DataFrame([item[0] for item in results]).astype(int)
    df1['S/T-P distance from modified cysteine'] = pd.Series([item[1] for item in results])
    #source data
    percent = counts.sum(axis = 0)/len(df1.accession_res)
    labels = ['1-5','6-10','11-15','16-20','21-25','26-30', '31-35', '36-40','41-45','46-50']
    graph_data = pd.DataFrame({'data': percent, 'label': labels})
    return graph_data, df1

def fig4f_stp_plot(df, title, sub):
    data, lig_df = fig4f_stp_prep(df)
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
    plt.savefig('{}/figures/fig4/{}_fig4f_stp_{}.pdf'.format(dir,date, sub), bbox_inches='tight',transparent = True)
    plt.close()

def fig4f_stp_generate():
    # lig_df = pd.read_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig4de.csv'.format(dir))
    df1 = fig4f_stp_prep(df)[1]
    lig_df = fig4d_fig4e_prep()[2]

    #THIS IS NOT WORKING
    slim_lig = lig_df[['accession_res', 'liganded_in_any']]
    keep_cols = [col for col in df1.columns if 'combined_' in col]
    df2 = df1[protein_cols + keep_cols + ['S/T-P distance from modified cysteine']]
    df2['accession_res'] = df2['accession'] + '_C' + df2.gene_res.str.split('_',expand = True)[1] 
    df2 = df2.merge(slim_lig, on = 'accession_res', how = 'left')
    am_cond = (df2[y_med] >= asynch_cutoff) | (df2[y_med] <= 1/asynch_cutoff)
    down = df2.loc[(df2[x_med] <= 1/lpp_cutoff) & am_cond]
    up = df2.loc[(df2[x_med] >= lpp_cutoff) & am_cond]
    df2['Reactivity in LPP(-)'] =  'Unchanged'
    df2.loc[df2.accession_res.isin(list(set(down.accession_res))), 'Reactivity in LPP(-)'] = 'Decreased'
    df2.loc[df2.accession_res.isin(list(set(up.accession_res))), 'Reactivity in LPP(-)'] = 'Increased'
    df2['Liganded in previous studies?'] =  df2['liganded_in_any'].copy()
    artifact = down.loc[(down[adapted_med]> (1/unchanging_cutoff)), 'accession_res']
    authentic = down.loc[down[adapted_med] <= 1/adapted_cutoff, 'accession_res']
    unassigned = down.loc[((down[adapted_med]> (1/adapted_cutoff)) & \
                          (down[adapted_med]<= (1/unchanging_cutoff))),'accession_res']
    df2.loc[df2.accession_res.isin(list(set(unassigned))),'Authentic or Artifactual decrease in LPP(-)?'] = 'Unassigned'
    df2.loc[df2.accession_res.isin(list(set(authentic))),'Authentic or Artifactual decrease in LPP(-)?'] = 'Authentic'
    df2.loc[df2.accession_res.isin(list(set(artifact))),'Authentic or Artifactual decrease in LPP(-)?'] = 'Artifact'
    df2.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig4f.csv'.format(dir), index = False, float_format = '%.2f')
    
    only_lig = lig_df.loc[lig_df.liganded_in_any == 1]
    lig_change = only_lig.loc[(only_lig[y_med] >= asynch_cutoff) | (only_lig[y_med] <= 1/asynch_cutoff)]
    lig_down = lig_change.loc[lig_change[x_med] <= 1/lpp_cutoff]
    lig_no_ppep = lig_down.loc[~(lig_down[adapted_med] > 1/unchanging_cutoff)]
    lig_up = lig_change.loc[lig_change[x_med] >= lpp_cutoff]
    fig4f_stp_plot(only_lig, 'All liganded cysteines', 'liganded')
    fig4f_stp_plot(lig_no_ppep, 'Decrease', 'liganded_noppep')
    fig4f_stp_plot(lig_up, 'Increase', 'liganded_up')
fig4f_stp_generate()

def fig4h_boxplots(protein_name):
    gene1 = df.loc[df.protein == protein_name]
    sdf = gene1.copy()
    drop_cols = [col for col in sdf.columns if 'category' in col] + [col for col in sdf.columns if 'protein_' in col]
    sdf.drop(drop_cols, axis =1, inplace = True)
    sdf.set_index(['accession', 'description', 'protein','asynch_stoichiometry_mann','mitosis_stoichiometry_mann','sequence','accession_res','gene_res'], inplace = True)
    median_cols = [col for col in sdf.columns if 'median' in col]
    sdf[median_cols] = np.log2(sdf[median_cols])
    sdf.to_excel('{}/supplemental_tables/source_data/Kemper_etal_2021_source_fig4h_{}.xlsx'.format(dir, protein_name), float_format = '%.2f')
    x_name = x_med.rsplit(sep = '_', maxsplit = 2)[0]
    y_name = y_med.rsplit(sep = '_', maxsplit = 2)[0]
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
    fig, axs = plt.subplots( figsize=(3.2, 2.4))
    order = [y_name, x_name, adapted_name]
    o = sns.color_palette()[1]
    # sns.stripplot(x= 'gene_res', y = 'value', hue = 'variable', alpha = 0.7, linewidth = 0.5, edgecolor = 'gray', hue_order = order, data = gene, jitter = True, dodge=True, palette = ['b','g','y'], ax = axs)
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
    plt.savefig('{}/figures/fig4/{}_fig4h_bar_{}.pdf'.format(dir, date, protein_name), bbox_inches='tight',transparent = True)    # plt.close()
    
fig4h_boxplots('KLC2')


