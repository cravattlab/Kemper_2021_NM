#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 12:46:09 2021

@author: estherkemper
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 08:44:06 2021

@author: estherkemper
"""
import numpy as np
import pandas as pd
import os

dir = os.getcwd()

#set a target cutoff defined in table

date = '20210830'

# protein_cols = ['accession','description', 'protein', 'gene_res']
protein_cols = ['accession','accession_res','description', 'protein', 'gene_res']
# phospho_col = ['pos_seq_mann', 'asynch_stoichiometry_mann','mitosis_stoichiometry_mann', 'stoichiometry','pos_seq_olsen', 'asynch_stoichiometry_olsen','mitosis_stoichiometry_olsen']
x_med = 'cysteine_Mitosis/Asynch_combined_median'
y_med = 'protein_Mitosis/Asynch_combined_median'
x_no = 'cysteine_Mitosis/Asynch_combined_no'
basal_reps = ['Exp476', 'Exp483', 'Exp489', 'Exp494']
cat_col ='cysteine_Mitosis/Asynch_category'

expt_table = pd.read_csv('{}/inputs/protein_inputs.csv'.format(dir), index_col = None)

target_cutoff = 2
lpp_cutoff = 2
asynch_cutoff = 1.6
unchanging_cutoff = 1.5
protein_cutoff = 1.6

df = pd.read_excel('{}/supplemental_tables/{}_fig1_combined.xlsx'.format(dir,date), sheet_name = 'filtered')


def peptide_changes(df):
    #input is df
    df['cysteine/protein'] = df[x_med]/df[y_med]
    #standard, high or low cysteines
    std_cond = (df[x_med] <= (1/target_cutoff)) |  (df[x_med] >= target_cutoff)
    #protein or all median changing, but one cysteine not and also 2x diff from the rest, cystesine not changing
    #can also use median value fo all quantified cysteines if over 5 quantified perprotein
    prot_cond = (df[y_med] <= (1/target_cutoff)) |  (df[y_med] >= target_cutoff)
    prot_cond1 = ((df[x_med+'_AGG'] <= (1/target_cutoff)) |  (df[x_med+'_AGG'] >= target_cutoff)) & (df['count']>= 5)
    cys_cond = (df[x_med] > (1/unchanging_cutoff)) &  (df[x_med] < unchanging_cutoff)
    cys_cond1 = (df['cysteine/protein'] >= 2) | (df['cysteine/protein'] <= 0.5)
    # cys_cond2 =  (df['count']>= 5) & ((df[x_med]/df[x_med+'_AGG'] >= 2) | ((df[x_med]/df[x_med+'_AGG'] <= 0.5))) 
    cys_cond3 = (df[x_med]/df['max'] <= 0.5) | (df[x_med]/df['min'] >= 2)
    unchanging_cond = (cys_cond) & (cys_cond3) & ((prot_cond) | (prot_cond1))
    # am_diff = df.loc[(unchanging_cond)|std_cond]
    am_diff = df.loc[std_cond]
    return am_diff
    
def protein_category_prep(df):
    #input is df
    #input is df
    #prepping for reactivity/expression stuff
    df_pep_count = df.groupby(['accession'])[x_med].count()
    med = df.groupby(['accession', 'protein', 'description'])[x_med].agg([np.median, np.max, np.min, 'count']).\
    rename(columns={'median': x_med, 'amax':'max', 'amin':'min'})
    med['max/min'] = med['max']/med['min']
    df = df.merge(med, on = 'accession', how = 'left', suffixes = ('','_AGG'))
    am_diff = peptide_changes(df)

    prot_ids = list(set(am_diff['accession']))
    prots = df.loc[df.accession.isin(prot_ids)].copy()
    # prots = prots.merge(med, on = 'accession', how = 'left', suffixes = ('','_AGG'))
    prots.set_index(['accession','description','protein','cysteine_sequence','gene_res'], inplace = True)
    prots.reset_index(inplace = True)
    prots.loc[:,'x/median'] = prots[x_med]/prots[x_med + '_AGG']
    # y_df = df[['accession', y_med]].drop_duplicates()
    # pep_count = y_df.merge(med, on = 'accession', how = 'right', suffixes = ('','_AGG'))
    # am_diff = am_diff.merge(med, on = 'accession', how = 'left', suffixes = ('', '_AGG'))
    return df, am_diff, prots

def reactivity(row):
    """assigning reactivity change for each individual cysteine"""
    #x_med is cysteine, y_med is protein, _agg is the aggregate of all cysteine values
    med_agg = x_med + '_AGG'
    if row['count'] > 4:
        #low changes: AND min/median AND max unchanging
        if (row[x_med] <= (1/lpp_cutoff)) and (row['max'] > (1/unchanging_cutoff)) and (row[x_med]/row[med_agg] <= (1/lpp_cutoff)):
            if row[y_med] != row[y_med]:
                return row['accession_res']
            elif (row[x_med]/row[y_med] <= (1/lpp_cutoff)):
                return row['accession_res']
        #high changes: AND max/median
        elif (row[x_med] >= lpp_cutoff) and (row['min'] < unchanging_cutoff) and (row[x_med]/row[med_agg] >= lpp_cutoff):
            if row[y_med] != row[y_med]:
                return row['accession_res']
            elif (row[x_med]/row[y_med] >= lpp_cutoff):
                return row['accession_res']
        #if not changing but rest are
        elif (row[x_med] >= 1/unchanging_cutoff) and (row[x_med] <= unchanging_cutoff): 
            if (row[med_agg] < 1/unchanging_cutoff) or (row[med_agg] > unchanging_cutoff): 
                if (row[x_med]/row[med_agg] >= lpp_cutoff) or ((row[x_med]/row[med_agg] <= lpp_cutoff)):
                    if row[y_med] != row[y_med]:
                        return row['accession_res']
                    elif (row[x_med]/row[y_med] <= (1/lpp_cutoff)) or (row[x_med]/row[y_med] >= lpp_cutoff):
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
            elif (row[x_med]/row[y_med] <= (1/lpp_cutoff)):
                return row['accession_res']
        #high changes: AND (max/median OR max/protein)
        elif (row[x_med] >= lpp_cutoff) and (row['min'] < unchanging_cutoff):
            if (row[x_med]/row[med_agg] >= lpp_cutoff):
                return row['accession_res']
            elif (row[x_med]/row[y_med] >= lpp_cutoff):
                return row['accession_res']
    if row['count'] ==2:
        #low changes: max/min AND 1 unchanging OR max/protein 
        if (row[x_med] <= (1/lpp_cutoff)) and (row['max'] > (1/unchanging_cutoff))  and (row['max']/row['min'] >=2):
            return row['accession_res']
        elif (row[x_med] <= (1/lpp_cutoff)) and (row[x_med]/row[y_med] <= (1/lpp_cutoff)) and (row[y_med] > (1/unchanging_cutoff)):
            return row['accession_res']
        #high changes: max/min AND 1 unchanging AND ()
        elif (row[x_med] >= lpp_cutoff) and (row['min'] < unchanging_cutoff) and (row['max']/row['min'] >=2):
            return row['accession_res']
        elif (row[x_med] >= lpp_cutoff) and (row[x_med]/row[y_med] >= lpp_cutoff) and (row[y_med] < unchanging_cutoff):
            return row['accession_res']
    elif row['count'] ==1:
        #low changes: max/min AND 1 unchanging OR max/protein 
        if (row[x_med] <= (1/lpp_cutoff))and (row['min']/row[y_med] <= (1/lpp_cutoff)) and (row[y_med] > (1/unchanging_cutoff)):
            return row['accession_res']
        #high changes: max/min AND 1 unchanging AND ()
        elif (row[x_med] >= lpp_cutoff) and (row['max']/row[y_med] >= lpp_cutoff) and (row[y_med] < unchanging_cutoff):
            return row['accession_res']

# def expression(row):
#     """assigning expression change to individual cysteine"""
#     med_agg = x_med + '_AGG'
#     if row['count'] > 4:
#         # # low changes: AND max 1.5 changing
#         # if (row[x_med] <= (1/lpp_cutoff)) and (row['max'] < (1/unchanging_cutoff)):
#         #     return row['accession_res']
#         # #high changes: AND min 1.5 changing
#         # elif (row[x_med] >= lpp_cutoff) and (row['min'] > unchanging_cutoff):
#         #     return row['accession_res']
#         #low changes: med OR protein expression changing (if available)
#         if (row[x_med] >= lpp_cutoff) and (row[x_med]/row[med_agg] < (lpp_cutoff)):
#             if  (row[med_agg] >= lpp_cutoff): 
#                 return row['accession_res']
#             elif  (row[y_med] >= lpp_cutoff): 
#                 return row['accession_res']
#         #high changes: med OR protein expression changing (if available)
#         elif (row[x_med] <= 1/lpp_cutoff) and (row[x_med]/row[med_agg] > (1/lpp_cutoff)):
#             if (row[med_agg] <= 1/lpp_cutoff):
#                 return row['accession_res']
#             elif (row[y_med] <= 1/lpp_cutoff):
#                 return row['accession_res']
#     elif (row['count']<=4) and (row['count']>= 3):
#         #low changes: if aggg of peptides OR protein 1.5 changing
#         if row[x_med] <= (1/lpp_cutoff) and (row[med_agg] > (1/lpp_cutoff)):
#             # if row['max'] < (1/unchanging_cutoff):
#             #     return row['accession_res']
#             if row[y_med] <= (1/protein_cutoff) or row[med_agg] <= (1/protein_cutoff):
#                 return row['accession_res']
#         #high changes: AND min 1.5 changing OR protein is 1.5 changing
#         elif row[x_med] >= lpp_cutoff and (row[x_med]/row[med_agg] < (lpp_cutoff)):
#             # if row['min'] > unchanging_cutoff:
#             #     return row['accession_res']
#             if row[y_med] >= protein_cutoff or row[med_agg] >= protein_cutoff and (row[x_med]/row[y_med] < (lpp_cutoff)):
#                   return row['accession_res']
#     if row['count'] <= 2:
#         #low changes: both changing 1.5 AND protein is changing
#         if (row[x_med] <= (1/lpp_cutoff)) and (row[y_med] <=(1/protein_cutoff)) and (row['max'] <= (1/unchanging_cutoff)) and (row[x_med]/row[med_agg] > (1/lpp_cutoff)):
#             return row['accession_res']
#         #high changes: max/min AND 1 unchanging AND ()
#         elif (row[x_med] >= lpp_cutoff) and (row[y_med] >= protein_cutoff) and (row['min'] >= unchanging_cutoff) and (row[x_med]/row[med_agg] < (lpp_cutoff)):
#             return row['accession_res']

def expression(row):
    """assigning expression change to individual cysteine"""
    med_agg = x_med + '_AGG'
    if row['count'] > 4:
        # # low changes: AND max 1.5 changing
        # if (row[x_med] <= (1/lpp_cutoff)) and (row['max'] < (1/unchanging_cutoff)):
        #     return row['accession_res']
        # #high changes: AND min 1.5 changing
        # elif (row[x_med] >= lpp_cutoff) and (row['min'] > unchanging_cutoff):
        #     return row['accession_res']
        #low changes: med OR protein expression changing (if available)
        if (row[x_med] >= lpp_cutoff) and (row[x_med]/row[med_agg] < (lpp_cutoff)):
            if  (row[med_agg] >= lpp_cutoff): 
                return row['accession_res']
            elif  (row[y_med] >= lpp_cutoff): 
                return row['accession_res']
        #high changes: med OR protein expression changing (if available)
        elif (row[x_med] <= 1/lpp_cutoff) and (row[x_med]/row[med_agg] > (1/lpp_cutoff)):
            if (row[med_agg] <= 1/lpp_cutoff):
                return row['accession_res']
            elif (row[y_med] <= 1/lpp_cutoff):
                return row['accession_res']
    elif (row['count']<=4) and (row['count']>= 3):
        #low changes: AND max 1.5 changing OR protein 1.5 changing
        if row[x_med] <= (1/lpp_cutoff)and (row[x_med]/row[med_agg] > (1/lpp_cutoff)):
            if row['max'] < (1/unchanging_cutoff):
                return row['accession_res']
            elif row[y_med] <= (1/protein_cutoff)and (row[x_med]/row[y_med] > (1/lpp_cutoff)):
                return row['accession_res']
        #high changes: AND min 1.5 changing OR protein is 1.5 changing
        elif row[x_med] >= lpp_cutoff and (row[x_med]/row[med_agg] < (lpp_cutoff)):
            if row['min'] > unchanging_cutoff:
                return row['accession_res']
            elif row[y_med] > protein_cutoff and (row[x_med]/row[y_med] < (lpp_cutoff)):
                  return row['accession_res']
    if row['count'] <= 2:
        #low changes: both changing 1.5 AND protein is changing
        if (row[x_med] <= (1/lpp_cutoff)) and (row[y_med] <=(1/protein_cutoff)) and (row['max'] <= (1/unchanging_cutoff)) and (row[x_med]/row[med_agg] > (1/lpp_cutoff)):
            return row['accession_res']
        #high changes: max/min AND 1 unchanging AND ()
        elif (row[x_med] >= lpp_cutoff) and (row[y_med] >= protein_cutoff) and (row['min'] >= unchanging_cutoff) and (row[x_med]/row[med_agg] < (lpp_cutoff)):
            return row['accession_res']

def run_protein_categories(df):
    #get dataframes for reactivity, expression and uncharacterized of proteins with changing cysteines
    df, am_diff, prots = protein_category_prep(df)
    react_id = prots.apply(reactivity, axis =1)
    react_id_df = pd.DataFrame({'accession_res':react_id} )
    react_id_df = react_id_df.assign(react = 'REACTIVITY')
    prots = prots.merge(react_id_df, on = 'accession_res', how = 'left')
    react_id = react_id_df.accession_res.str.rsplit('_', n = 2, expand = True)[0]
    all_react_df = prots.loc[prots.accession.isin(react_id)]
    print(str(len(set(all_react_df.accession))) + ' proteins with reactivity changes')
    exp_id = prots.apply(expression, axis =1)
    exp_id_df = pd.DataFrame({'accession_res':exp_id} )
    exp_id_df = exp_id_df.assign(exp = 'EXPRESSION')
    prots = prots.merge(exp_id_df, on = 'accession_res', how = 'left')
    changing = ((prots[x_med] >= 2) | (prots[x_med] <= 0.5))
    prots.loc[((prots['exp'].isnull()) & (prots['react'].isnull())) & changing, ['exp', 'react']] = 'Unassigned'
    exp_id = exp_id_df.accession_res.str.rsplit('_', n = 2, expand = True)[0]
    all_exp_df = prots.loc[prots.accession.isin(exp_id)]
    print(str(len(set(all_exp_df.accession))) + ' proteins with cys that change in same direction')
    uncharacterized = prots.loc[~(prots.accession.isin(exp_id)) & ~(prots.accession.isin(react_id))]
    exp_list = pd.DataFrame({'accession': list(set(exp_id))})
    react_list = pd.DataFrame({'accession': list(set(react_id))})
    print(str(len(set(prots.accession))) + ' total')
    print(str(len(set(uncharacterized.accession))) + ' uncharacterized')
    all_exp_df.set_index(['accession','description', 'protein','cysteine_sequence','gene_res'], inplace = True)
    all_exp_df = all_exp_df.sort_values(['protein'], ascending = True)
    all_react_df.set_index(['accession','description', 'protein','cysteine_sequence','gene_res'], inplace = True)
    all_react_df = all_react_df.sort_values(['protein'], ascending = True)
    prots.loc[prots.protein == 'Uncharacterized', 'protein'] = 'Uncharacterized_' + prots.accession
    prots_df = prots.sort_values(['protein'], ascending = True)
    prots_df = prots_df.assign(accession_res = prots_df.accession_res.str.rsplit('_', n = 1, expand = True)[0])
    prots_df.set_index(['accession','description', 'protein','accession_res','cysteine_sequence','gene_res'], inplace = True)
    pep = [col for col in prots_df.columns if 'no_unique' in col]
    pep = pep + [col for col in prots_df.columns if 'all_peptides' in col]
    go_df = prots_df.copy()
    go_df[['react', 'exp']] = go_df[['react', 'exp']].fillna('')
    # go_df[['react', 'exp']] = go_df[['react', 'exp']].replace('Unassigned', '')
    go_df = go_df.loc[~((go_df.react == '') & (go_df.exp == ''))]
    go_df = go_df.assign(new = go_df['react']+ go_df['exp'] )
    go_df1 = pd.DataFrame(go_df.groupby('accession')['new'].apply(sum)).reset_index()
    #go_no_exp contains what we will use for GO analysis,  is what is taken out
    go_react_only = go_df1.loc[go_df1.new.str.contains('REACTIVITY')]
    go_reactivity_changing = pd.DataFrame({'accession': list(set(go_react_only.accession))})
    go_reactivity_changing.to_csv('{}/for_figures/suppfig1/GO_analysis/inputs/suppfig1_changing_REACTONLY.txt'.format(dir), index = False)
    background_list = pd.DataFrame({'accession':list(set(df.accession))})
    background_list.to_csv('{}/for_figures/suppfig1/GO_analysis/inputs/suppfig1_background.txt'.format(dir), index = False)

    return prots_df
    
def protein_main(df):
    prots_df = run_protein_categories(df)
    combined_cols = [col for col in prots_df.columns if 'combined_median' in col]
    category = ['react', 'exp']
    TMT_reps = [col for col in prots_df.columns if 'median_TMT' in col]
    CV = [col for col in prots_df.columns if '_CV' in col]
    no_cols = [col for col in prots_df.columns if 'combined_no' in col]
    sdf = prots_df[category + combined_cols + TMT_reps + CV + no_cols]
    sdf.rename({'react':'Reactivity-based change?', 'exp': 'Expression-based change?'}, axis = 1, inplace = True)
    sdf[TMT_reps + combined_cols] = np.log2(sdf[TMT_reps + combined_cols])
    sdf[[x_med, y_med]] = sdf[[x_med, y_med]].fillna('NQ')
    # misc = prot_change_cys_constant(df)
    print('Writing protien data to excel')
    with pd.ExcelWriter('{}/supplemental_tables/{}_fig1_proteins.xlsx'.format(dir,date)) as writer:  # doctest: +SKIP
        # am_diff.to_excel(writer, sheet_name = 'diff_pep', float_format='%.2f', index = False)
        # all_react_df.to_excel(writer, sheet_name = 'reactivity', float_format='%.2f')
        # all_exp_df.to_excel(writer, sheet_name = 'expression', float_format='%.2f')
        # uncharacterized.to_excel(writer, sheet_name = 'uncharacterized', float_format='%.2f')
        prots_df.to_excel(writer, sheet_name = 'all prots with changes', float_format = '%.2f')
        sdf.to_excel(writer, sheet_name = 'supplemental', float_format='%.2f')
#     return all_exp_df, all_react_df, prots_df

# def protein_figures(df):
#     all_exp_df, all_react_df, prots_df = protein_main(df)
#     prots_df.reset_index(inplace = True)

protein_main(df)
