#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 11:28:16 2021
This creates the categories for analysis and protein pages
@author: estherkemper
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt, seaborn as sns
from matplotlib import rcParams
plt.rcParams['pdf.fonttype'] = 42
font = {'family' : 'ARIAL',
        'weight' : 'normal',
        'size'   : 10}         
plt.rc('font', **font) #set the font style created
sns.set_style("white" )
plt.rcParams.update({'text.color' : "black"})

dir = os.getcwd()
#set a target cutoff defined in table
date = '20210830'

lpp_cutoff = 2
asynch_cutoff = 1.6 
unchanging_cutoff = 1.6
protein_cutoff = 1.6

protein_cols = ['accession','description', 'protein', 'accession_res', 'gene_res']
phospho_col = [ 'asynch_stoichiometry_mann','mitosis_stoichiometry_mann']

x_med ='Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir,date), sheet_name = 'filtered')
lpp = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')

print('Finished reading input data')

def peptide_main(df):
    #input df3
    #changing in LPP/CTRL, changing in Asynch/Mitosis
    cond1 = (df[x_med] <= (1/lpp_cutoff)) |  (df[x_med] >= lpp_cutoff)
    cond2 = (df[y_med] <= (1/asynch_cutoff)) |  (df[y_med] >= asynch_cutoff)
    #changing in basal only
    am_diff = df.loc[(cond1) & (cond2)]
    return am_diff
    
def protein_category_prep(df):
    #protein pages
    #input is df3
    df_pep_count = df.groupby(['accession'])[x_med].count()
    med = df.groupby(['accession', 'protein', 'description'])[x_med].agg([np.median, np.max, np.min, 'count']).\
    rename(columns={'median': x_med, 'amax':'max', 'amin':'min'})
    med['max/min'] = med['max']/med['min']
    df = df.merge(med, on = 'accession', how = 'left', suffixes = ('','_AGG'))
    am_diff = peptide_main(df)
    prot_ids = list(set(am_diff['accession']))
    prots = df.loc[df.accession.isin(prot_ids)].copy()
    prots.set_index(['accession','description','protein']+phospho_col+['sequence','gene_res'], inplace = True)
    prots.reset_index(inplace = True)
    prots.loc[:,'x/median'] = prots[x_med]/prots[x_med + '_AGG']
    # y_df = df[['accession', x_med]].drop_duplicates()
    # pep_count = y_df.merge(med, on = 'accession', how = 'right', suffixes = ('','_AGG'))
    return prots

def reactivity(row):
    med_agg = x_med + '_AGG'
    if row['count'] > 4:
        #low changes: AND min/median AND max unchanging
        if (row[x_med] <= (1/lpp_cutoff)) and (row['max'] > (1/unchanging_cutoff)) and (row[x_med]/row[med_agg] <= (1/lpp_cutoff)):
            return row['accession_res']
        #high changes: AND max/median
        elif (row[x_med] >= lpp_cutoff) and (row['min'] < unchanging_cutoff) and (row[x_med]/row[med_agg] >= lpp_cutoff):
            return row['accession_res']
        #if not changing but rest are
        elif (row[x_med] >= 1/unchanging_cutoff) and (row[x_med] <= unchanging_cutoff): 
            if (row[med_agg] < 1/unchanging_cutoff) or (row[med_agg] > unchanging_cutoff): 
                if (row[x_med]/row[med_agg] >= lpp_cutoff) or ((row[x_med]/row[med_agg] <= lpp_cutoff)):
            #not changing but rest are
                    return row['accession_res']
                # if row['max']>(1/unchanging_cutoff)  and row[x_med] <= (1/lpp_cutoff) and row['max']/row[x_med] >= lpp_cutoff:
                #     return row['accession_res']
                # elif row['min']<(unchanging_cutoff)  and row[x_med] >= lpp_cutoff and row['min']/row[x_med] <= (1/lpp_cutoff):
                #     return row['accession_res']
    elif (row['count']<=4) & (row['count']>= 3) & (row['max']/row['min'] >=2):
        #low changes: max/min >=2 AND 1 unchanging AND (min/median OR or min/protein) 
        if (row[x_med] <= (1/lpp_cutoff)) and (row['max'] > (1/unchanging_cutoff)):
            if (row[x_med]/row[med_agg] <= (1/lpp_cutoff)):
                return row['accession_res']
        #high changes: AND (max/median )
        elif (row[x_med] >= lpp_cutoff) and (row['min'] < unchanging_cutoff):
            if (row[x_med]/row[med_agg] >= lpp_cutoff):
                return row['accession_res']
    if row['count'] ==2:
        #low changes: max/min AND 1 unchanging 
        if (row[x_med] <= (1/lpp_cutoff)) and (row['max'] > (1/unchanging_cutoff))  and (row['max']/row['min'] >=2):
            return row['accession_res']
        #high changes: max/min AND 1 unchanging AND ()
        elif (row[x_med] >= lpp_cutoff) and (row['min'] < unchanging_cutoff) and (row['max']/row['min'] >=2):
            return row['accession_res']

def expression(row):
    med_agg = x_med + '_AGG'
    if row['count'] > 4:
        # low changes: AND max 1.5 changing
        if (row['min'] <= (1/lpp_cutoff)) and (row['max'] < (1/unchanging_cutoff)):
            return row['accession_res']
        #high changes: AND min 1.5 changing
        elif (row['max'] >= lpp_cutoff) and (row['min'] > unchanging_cutoff):
            return row['accession_res']
        # if (row[x_med] >= lpp_cutoff) or (row[x_med] <= 1/lpp_cutoff):
        #     if  (row[med_agg] >= lpp_cutoff) or (row[med_agg] <= 1/lpp_cutoff): 
        #         return row['accession_res']
    elif (row['count']<=4) and (row['count']>= 3):
        #low changes: AND max 1.5 changing OR protein 1.5 changing
        if row[x_med] <= (1/lpp_cutoff):
            if row['max'] < (1/unchanging_cutoff):
                return row['accession_res']
        #high changes: AND min 1.5 changing OR protein is 1.5 changing
        elif row[x_med] >= lpp_cutoff:
            if row['min'] > unchanging_cutoff:
                return row['accession_res']
    if row['count'] == 2:
        #low changes: both changing 1.5 
        if (row[x_med] <= (1/lpp_cutoff)) and (row['max'] <= (1/unchanging_cutoff)):
            return row['accession_res']
        #high changes: max/min AND 1 unchanging AND ()
        elif (row[x_med] >= lpp_cutoff) and (row['min'] >= unchanging_cutoff):
            return row['accession_res']

def run_protein_categories(df):
    prots = protein_category_prep(df)
    react_id = prots.apply(reactivity, axis =1)
    react_id_df = pd.DataFrame({'accession_res':react_id} )
    react_id_df = react_id_df.assign(react = 'Site specific')
    prots = prots.merge(react_id_df, on = 'accession_res', how = 'left')
    react_id = react_id_df.accession_res.str.split('_', n = 2, expand = True)[0]
    all_react_df = prots.loc[prots.accession.isin(react_id)]
    print(str(len(set(all_react_df.accession))) + ' proteins with reactivity changes')
    exp_id = prots.apply(expression, axis =1)
    exp_id_df = pd.DataFrame({'accession_res':exp_id} )
    exp_id_df = exp_id_df.assign(exp = 'Not site specific')
    prots = prots.merge(exp_id_df, on = 'accession_res', how = 'left')
    changing = ((prots[x_med] >= 2) | (prots[x_med] <= 0.5))
    prots.loc[((prots['exp'].isnull()) & (prots['react'].isnull())) & changing, ['react','exp']] = 'Unassigned'
    exp_id = exp_id_df.accession_res.str.split('_', n = 2, expand = True)[0]
    all_exp_df = prots.loc[prots.accession.isin(exp_id)]
    print(str(len(set(all_exp_df.accession))) + ' proteins with cys that change in same direction')
    one_pep = prots.groupby('accession').count()
    one_pep = prots.groupby('accession', as_index = False)['gene_res'].count()
    one_pep_list = one_pep.loc[one_pep.gene_res == 1, 'accession']
    prots.loc[(prots.accession.isin(list(set(one_pep_list)))),['react','exp']] = 'Unclear'
    unclear = prots.loc[prots.accession.isin(one_pep_list)].sort_values(['protein'], ascending = True)
    print(str(len(set(unclear.accession))) + ' proteins with 1 cys')
    print(str(len(set(prots.accession))) + ' total')
    prots.loc[prots.protein == 'Uncharacterized','protein'] = 'Uncharacterized_' + prots.accession
    prots_df = prots.sort_values(['protein'], ascending = True)
    return prots_df

def lpp_merge(df):
    prots_df = run_protein_categories(df)
    p_med = 'cysteine_' + x_med
    lpp1 = lpp[['accession', 'gene_res', p_med]]
    lpp2 = pd.DataFrame(lpp1.groupby('accession')['gene_res'].apply(lambda row: ", ".join(row))).reset_index()
    lpp2.rename({'gene_res':'LPP-sensitive phosphosite'}, axis =1, inplace = True)
    prots_df1 = prots_df.merge(lpp2, on = 'accession', how = 'left')
    # prots_df1['protein'] = str(prots_df1.protein + ';' + prots_df1['LPP-sensitive phosphosite']
    prots_df1['LPP-sensitive phosphosite'] = prots_df1['LPP-sensitive phosphosite'].astype(str)
    prots_df1['LPP-sensitive phosphosite'] = prots_df1['LPP-sensitive phosphosite'].replace('nan','')
    prots_df1 = prots_df1.assign(accession_res = prots_df1.accession_res.str.rsplit('_', n=1, expand = True)[0])
    prots_df1.set_index(['accession','description','protein', 'LPP-sensitive phosphosite', 'asynch_stoichiometry_mann','mitosis_stoichiometry_mann','accession_res','sequence','gene_res'], inplace = True)
    
    return prots_df1

def protein_main(df):
    am_diff = peptide_main(df)
    prots_df = lpp_merge(df)
    combined_cols = [col for col in prots_df.columns if 'combined_median' in col]
    TMT_reps = [col for col in prots_df.columns if 'median_TMT' in col]
    CV = [col for col in prots_df.columns if '_CV' in col]
    no_cols = [col for col in prots_df.columns if 'combined_no' in col]
    sdf = prots_df[combined_cols + CV + no_cols+ TMT_reps] #'LPP-sensitive phosphosite'] + phospho_col + 
    sdf.rename({'react':'Reactivity-based change?', 'exp': 'Expression-based change?'}, axis = 1, inplace = True)
    sdf[TMT_reps + combined_cols] = np.log2(sdf[TMT_reps + combined_cols])
    sdf[combined_cols] = sdf[combined_cols].fillna('NQ')
    print(sdf[combined_cols])

    print('Writing protien data to excel')
    with pd.ExcelWriter('{}/supplemental_tables/{}_fig2_proteins.xlsx'.format(dir, date)) as writer:  # doctest: +SKIP
        # am_diff.to_excel(writer, sheet_name = 'LPP_diff-pep', index = False)
        prots_df.to_excel(writer, sheet_name = 'all prots with changes', float_format = '%.2f')
        sdf.to_excel(writer, sheet_name = 'supplemental', float_format = '%.2f')
    
if __name__ == '__main__':    
    protein_main(df)
    
    

    
   

    
    


