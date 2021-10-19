#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 13:28:58 2021

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
denature_med = 'Mitosis_Native/Denatured_combined_median'

denature_df = pd.read_excel('{}/supplemental_tables/{}_fig2_denature_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
dis_map = pd.read_csv('{}/other_datasets/IUPRED/20210303_all_fasta.csv'.format(dir))
dis_map = dis_map.assign(accession_res = (dis_map.accession + '_C' + dis_map.accession_res.str.rsplit('_',expand = True)[1]))
dis_map = dis_map[['accession_res', 'IUPRED SCORE']]


def moving_ave(df):
    new_df = df.loc[df['IUPRED SCORE'] >= 0.5]
    counts = len(new_df.accession_res)
    return counts


def fig4f_stp_prep(df): 
    # df[debaty_med] = np.log2(df[y_med] )
    df = df.assign(accession_res = (df.accession + '_C' + df.gene_res.str.split('_',expand = True)[1]))
    df = df.loc[df[denature_med].notnull()]
    df1 = df.merge(dis_map, how = 'left', on = 'accession_res').reset_index()
    df1[denature_med] = np.log2(df1[denature_med])
    denom = df1.loc[df1['IUPRED SCORE'] >= 0.5]
    c1 = np.log2(2.5)
    c2 = np.log2(1.8)
    c3 = np.log2(1.4)
    x = denom[denature_med]
    g = denom.loc[x >=c1]
    f = denom.loc[(x < c1) & (x >= c2)]
    e = denom.loc[(x < c2) & (x >= c3)]
    d = denom.loc[(x < c3) & (x >= -c3)]
    c = denom.loc[(x < -c3) & (x >= -c2)]
    b = denom.loc[(x < -c2) & (x >= -c1)]
    a = denom.loc[x <= -c1]
    df_lst = [a, b, c, d, e, f,g]
    comb = pd.concat(df_lst)
    total_counts = [moving_ave(x) for x in df_lst]
    percent = [x/len(set(denom.accession_res)) for x in total_counts]
    labels = ['< -1.3', '> -1.3 and < -0.85', '> -0.85 and < -0.48','> -0.48 and < 0.48', '> 0.48 and < 0.84', '> 0.84 and < 1.32','> 1.3']
    graph_data = pd.DataFrame({'data': percent, 'label': labels})
    print(sum(total_counts), len(set(denom.accession_res)))
    fig, axs = plt.subplots(figsize=(2.1, 2.2))
    # fig, axs = plt.subplots()
    g = sns.barplot(x = 'label', y = 'data', color = 'black', data = graph_data, linewidth = 0.75)
    plt.ylabel('% of cysteines ',family = 'ARIAL', fontsize = 10)
    plt.xlabel('Native/Denatured cysteine reactivity\nlog2 MS3 reporter ion instensity',family = 'ARIAL', fontsize = 10)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    plt.xticks(fontsize = 10, family = 'ARIAL', rotation = 45, ha = 'right')
    plt.yticks(fontsize = 10, family = 'ARIAL')
    plt.title('Cysteines in predicted disordered regions', fontsize = 10, family = 'ARIAL' )
    # plt.ylim(0,0.7)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.savefig('{}/figures/suppfig2/potential/{}_denature_moving_ave_disordered.pdf'.format(dir,date), bbox_inches='tight',transparent = True)
    # plt.savefig('{}/figures/suppfig2/potential/{}_denature_moving_ave_structured.pdf'.format(dir,date), bbox_inches='tight',transparent = True)
    # plt.close()

fig4f_stp_prep(denature_df)
import matplotlib.patches as mpatches


# def scatter_plot(df, x_med, y_med):
#     # dis_map = fig4d_disorder_map()
#     # lig_map = fig4d_lig_map()
#     df1 = df.loc[~df[y_med].isnull()]
#     df1[y_med] = np.log2(df1[y_med] )
#     df1 = df1.assign(accession_res = (df1.accession + '_C' + df1.gene_res.str.split('_',expand = True)[1]))
#     df3 = df1.merge(dis_map, how = 'left', on = 'accession_res')
#     df3 = df3.sort_values(y_med).reset_index()
#     # df3.loc[df3['IUPRED SCORE'] >= 0.5, 'hue_color'] 
#     # df3.loc[df3['IUPRED SCORE'] < 0.5, 'hue_color'] = 'blue'
#     cond1 = df3['IUPRED SCORE'] >= 0.5
#     x = df3[y_med]
#     y = df3['IUPRED SCORE']
#     print(x, y)
#     plt.figure(figsize = (5,5))
#     plt.scatter(x[~cond1], y[~cond1], lw = 0, alpha= 0.5, c = '#404040', marker = '.')
#     plt.scatter(x[cond1], y[cond1], lw = 0, alpha = 0.5,c = 'blue', marker = '.')
#     left, bottom, width, height = (1, 0.49, 3.4, 0.5)
#     # rect=mpatches.Rectangle((left,bottom),width,height, 
#     #                         #fill=False,
#     #                         alpha=0.3,
#     #                        facecolor="red")
#     # plt.gca().add_patch(rect)
#     # plt.axvline(x = 1,  c = 'red')
#     plt.ylabel('IUPRED SCORE', fontsize = 10, family = 'ARIAL')
#     plt.xlabel('Native/Denatured', fontsize = 10, family = 'ARIAL')
#     plt.xticks(fontsize = 10, family = 'ARIAL')
#     plt.yticks(fontsize = 10, family = 'ARIAL')
#     plt.savefig('{}/figures/suppfig2/potential/{}_denature_IUPRED.pdf'.format(dir,date), transparent = True, bbox_inches='tight')
    
    # plt.scatter(x, y_med,  data = df3,lw = 0, alpha= 0.5,  marker = '.',)
    
scatter_plot(denature_df, x_med, denature_med)

# def supp_fig2b_correlation_plot(df, x_med, y_med):
#     #source data in main datase
#     print('Running Supp Fig2 correlation plot')
#     plt.subplots(figsize=(2.8,2.8)) 
#     df1= df.loc[~(df[y_med] == 'NQ')]
#     df1 = df1.loc[~(df1[x_med] == 'NQ')]
#     df1 = df1[['gene_res', 'accession', x_med,y_med,'phospho_bool']]
#     df1[[x_med,y_med]] = np.log2(df1[[x_med,y_med]].astype('float') )
#     x = df1[x_med]
#     y = df1[y_med]
#     # cond1 = (df1[x_med] >= np.log2(lpp_cutoff))
#     # cond2 = (df1[x_med] <= np.log2(1/lpp_cutoff))
#     # up_mask = (cond1) # | (cond3) 
#     # down_mask = (cond2)# | (cond4)
#     # both = (up_mask) | (down_mask)
#     # plt.scatter( y[~both],x[~both], lw = 0, alpha= 0.5, c = '#404040', marker = '.')
#     # plt.scatter( y[up_mask], x[up_mask],lw = 0, alpha = 0.5,c = 'blue', marker = '.')
#     # plt.scatter( y[down_mask], x[down_mask], lw = 0, alpha = 0.5,c = 'red', marker = '.')
#     both = df1['phospho_bool']
#     plt.scatter( y[both],x[both], lw = 0, alpha= 0.3, c = 'green', marker = '.')
#     plt.scatter( y[~both],x[~both], lw = 0, alpha= 0.3, c = '#404040', marker = '.')

#     plt.xlabel('Mitosis Native/Denatured\nlog2 MS3 reporter ion intensity', fontsize =10, family = 'ARIAL')
#     plt.ylabel('Mitosis LPP(-)/LPP(+)\nlog2 MS3 reporter ion intensity',fontsize =10, family = 'ARIAL')
#     # plt.axhline(y = 1, c = 'black', linestyle = 'dotted')
#     # plt.axhline(y = -1, c = 'black', linestyle = 'dotted')
#     # plt.axvline(x = 1, c = 'black', linestyle = 'dotted')
#     # plt.axvline(x = -1, c = 'black', linestyle = 'dotted')
#     plt.ylim(-5,5) 
#     plt.ylim(-5, 5)
#     plt.xticks(fontsize =10, family = 'ARIAL',ticks= [-4,-2,0,2,4])
#     plt.yticks(fontsize =10, family = 'ARIAL',ticks= [-4,-2,0,2,4])
#     plt.savefig('{}/figures/suppfig2/potential/{}_suppfig2b_rot_correlation4.pdf'.format(dir,date), bbox_inches='tight')
# supp_fig2b_correlation_plot(denature_df, x_med, denature_med)

# def fig4a_data_generation():
#     #proteins with high stoich phosphorylation, fig5
#     phospho_cols = [col for col in denature_df.columns if 'mitosis_stoichiometry' in col]
    
#     denature_df['phospho_bool'] = denature_df[phospho_cols].notnull().any(axis = 1)
    
#     t = denature_df.loc[(denature_df[x_med].notnull()) & (denature_df[denature_med].notnull())]
#     p = t.loc[t['phospho_bool']]
    
#     anti_corr = t.loc[(t[denature_med] >= 2) & (t[x_med] <= 0.5)]
#     phospho_anti_corr = anti_corr.loc[anti_corr['phospho_bool']]
    
#     anti_corr_down = t.loc[(t[denature_med] <= 0.5) & (t[x_med] >= 2)]
#     phospho_anti_corr_down = anti_corr_down.loc[anti_corr_down['phospho_bool']]
    
#     d_up = t.loc[(t[denature_med] >= 2) ]
#     phospho_d_up = d_up.loc[d_up['phospho_bool']]
    
#     d_down = t.loc[(t[denature_med] <= 0.5) ]
#     phospho_d_down = d_down.loc[d_down['phospho_bool']]
    
#     x_down = t.loc[t[x_med] <= 0.5]
#     p_x_down = x_down.loc[x_down['phospho_bool']]
       
#     x_up = t.loc[t[x_med] >=2]
#     p_x_up = x_up.loc[x_up['phospho_bool']]
    

    
#     # t = denature_df.copy()
#     # p = denature_df.loc[denature_df['phospho_bool']]
    
#     t0 = len(set(t.accession_res))
#     p0 = len(set(p.accession_res))
#     t1 = len(set(anti_corr.accession_res))
#     p1 = len(set(phospho_anti_corr.accession_res))
#     t2 = len(set(anti_corr_down.accession_res))
#     p2 = len(set(phospho_anti_corr_down.accession_res))
#     t3 = len(set(d_up.accession_res))
#     p3 = len(set(phospho_d_up.accession_res))
#     t4 = len(set(d_down.accession_res))
#     p4 = len(set(phospho_d_down.accession_res))
#     t5 = len(set(d_up.accession_res))
#     p5 = len(set(phospho_d_up.accession_res))
#     t6 = len(set(x_down.accession_res))
#     p6 = len(set(p_x_down.accession_res))
#     t7 = len(set(x_up.accession_res))
#     p7 = len(set(p_x_up.accession_res))
    
#     labels = ['0','1','2', '3','4', '5','6', '7']
#     p_number = [p0, p1, p2, p3, p4, p5, p6, p7]
#     t_number = [t0, t1, t2, t3, t4, t5, t6, t7]
    
#     labels = ['All quantified cys','Native/Denatured > 2\nAND LPP(-)/LPP(+) < 0.5','Native/Denatured > 2', 'LPP(-)/LPP(+) < 0.5']
#     p_number = [p0, p1, p3, p6]
#     t_number = [t0, t1, t3, t6]
    
#     y_value =  [i / j for i, j in zip(p_number, t_number)] 
#     percent_labels = [str(i)+ '/\n' + str(j) for i,j in zip(p_number, t_number)]    
#     plt.subplots(figsize=(3,2.2))
#     ax = sns.barplot(x = labels, y = y_value)
#     plt.ylabel('Percent of cysteines on proteins w/\nhigh stoichiometry phosphorylation', fontsize =10, family = 'ARIAL')
#     plt.xticks(fontsize =10, family = 'ARIAL', rotation = 45, ha = 'right')
#     plt.yticks(fontsize =10, family = 'ARIAL')
#     # plt.title('Predicted intrinsic disorder', fontsize =10, family = 'ARIAL')
#     plt.ylim(0,1)
#     plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
#     for i in range(len(percent_labels)):
#         plt.annotate(str(percent_labels[i]), xy=(i,y_value[i]), ha='center', va='bottom', fontsize = 10, family = 'ARIAL')
#     plt.savefig('{}/figures/suppfig2/potential/{}_denature_percent_phosphorylated1.pdf'.format(dir, date), bbox_inches='tight')
# fig4d_plot_disorder()
    

