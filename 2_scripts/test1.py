#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 10:58:16 2021

@author: ekk
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:03:56 2019
##This is for individual reps of mitosis +/- lpp (no asynch, no LPP after)
##LPP1 is the experimental group, just comparing CTRL1 and LPP1 
This combines everything together, then filters based on LPP/CTRL
@author: estherkemper
"""
import numpy as np
import pandas as pd
import os

dir = os.getcwd()

#set a target cutoff defined in table
target_cutoff = 2
date = '20210830'

phospho_location = '{}/final_phospho_mann_olsen.csv'.format(dir)
expt_table = pd.read_csv('{}/inputs/peptide_inputs.csv'.format(dir), index_col = None)
expt_table['ratio_name'] = expt_table.num + '/' + expt_table.denom + '_ratio_combined_median'

# protein_cols = ['accession','description', 'protein', 'gene_res']
protein_cols = ['unique','accession','description', 'protein', 'accession_res','gene_res']
lim_med = list(expt_table.loc[expt_table.group == 'limited_original', 'ratio_name'])[0]
before_med = list(expt_table.loc[expt_table.group == 'all_original', 'ratio_name'])[0]
before_no = before_med.rsplit('_', maxsplit = 1)[0] + '_no'
before_cat = before_med.rsplit('_', maxsplit = 2)[0] + '_category'
after_med = list(expt_table.loc[expt_table.group == 'adapted', 'ratio_name'])[0]
after_cat = after_med.rsplit('_', maxsplit = 2)[0] + '_category'

def combine():
    print('Reading input datasets')
    original_df = pd.read_csv('{}/1_scripts/combined/{}_combined_all_original.csv'.format(dir, date))
    limited_original_df = pd.read_csv('{}/1_scripts/combined/{}_combined_limited_original.csv'.format(dir, date))
    adapted_df = pd.read_csv('{}/1_scripts/combined/{}_combined_adapted.csv'.format(dir, date))
    basal_df = pd.read_csv('{}/1_scripts/combined/{}_combined_basal.csv'.format(dir, date))
    protein = pd.read_csv('{}/1_scripts/combined/{}_combined_protein_whole_proteome.csv'.format(dir, date))
    
    
    limited_original_df.set_index(protein_cols, inplace = True)
bar = [str(x) for x in list(range(0,len(foo)))]
    original
    
    
    lim_ratio_name = lim_med.rsplit('_', maxsplit = 2)[0]
    limited_original1_df.columns = limited_original1_df.columns.str.replace(lim_ratio_name,'Mitosis_CTRL1/Mitosis_LPP1_limited_before_ratio_combined_median')
    limited_original1_df.columns = limited_original1_df.columns.str.replace('LPP1_','Mitosis_LPP/CTRL_limited_before_')
    original1_df.columns = original1_df.columns.str.replace('LPP1_','Mitosis_LPP/CTRL_before_')
    adapted1_df.columns = adapted1_df.columns.str.replace('LPP3_','Mitosis_LPP/CTRL_after_')
    basal1_df.columns = basal1_df.columns.str.replace('Asynch_CTRL2','cysteine_Asynch/Mitosis')
    #ctrl1_df.columns = ctrl1_df.columns.str.replace('CTRL1_','Before/After_')
    protein_df.columns = protein_df.columns.str.replace('Asynch_CTRL2','protein_Asynch/Mitosis')
    limited_original1_df.reset_index(inplace = True)
    print('Merging data')
    df = original1_df.merge(adapted1_df, how = 'outer', on = protein_cols, suffixes = [' original', ' adapted'])
    df = df.merge(limited_original1_df, how = 'outer', on = protein_cols, suffixes = ['', ' original_LIMITED'])
    #df = df.merge(ctrl1_df, how = 'outer', on = protein_cols, suffixes = ['', ' original/adapted'])
    df = df.merge(basal1_df, how = 'outer', on = protein_cols, suffixes = ['', ' basal'])
    df = df.merge(protein_df, how = 'left', on = ['accession','description','protein'], suffixes = ['', ' protein'])
    return df

def sequence_combine(row):
     reps = [col for col in row.index if 'sequence' in col]
     sequences = row[reps].dropna()
     shortest_sequence =  min(sequences, key=len)
     return pd.Series(shortest_sequence)
 
def phospho_merge():
    df = combine()
    df.rename(columns = {'sequence':'sequence ctrl'}, inplace = True)
    #get rid of empty columns
    df = df.dropna(axis = 1, how = 'all')
    sequence_cols = [col for col in df.columns if 'sequence' in col]
    df['sequence'] = df[sequence_cols].apply(sequence_combine, axis = 1)
    df.drop(sequence_cols, axis = 1 , inplace = True)
    df['gene'] = df.gene_res.str.split(pat = '_', expand = True)[0]
    phospho_df = pd.read_csv('{}'.format(phospho_location))
    phospho_col = [col for col in phospho_df.columns if 'accession' not in col]
    df1 = df.merge(phospho_df, on = 'accession', how = 'left')
    print('Finished merging')
    return df1, phospho_col

def qc_checks(row):
    #if after value exists but corresponding value in before lPP/ctrl not there, do not interpret
    after_reps = [col for col in row.index if 'CTRL_after_ratio_Exp' in col]
    before_reps = [col for col in row.index if 'Mitosis_LPP/CTRL_before_ratio_Exp' in col]
    before = pd.to_numeric(row[before_med], errors = 'coerce')
    limited = pd.to_numeric(row[lim_med], errors = 'coerce')
    after = pd.to_numeric(row[after_med], errors = 'coerce')
    #for changing after, bu tnot changing before, with before not matching or nulll
    if ((after >=2) and (before <= 1.6)) or \
        ((after <= 0.5) and (before >= (1/1.6))):  
        if (limited != limited) or (type(limited) == str):
            return 'Delete'
        elif row[after_cat] == '0 or 1 rep':
            after_max = row[after_reps].notnull()
            exp_name = row[after_reps][after_max].index.to_list()[0]
            exp_name = exp_name.rsplit(sep = '_', maxsplit = 1)[1]
            before_value = [col for col in before_reps if exp_name in col]
            if row[before_value].item() != row[before_value].item():
                return 'Delete'
            elif (row[before_value].item()/after < 2) and (before/after >2 ):
                return 'Delete'
            else:
                return row[after_cat]
        else:
            return row[after_cat]
    #for unchanging after but changing before
    elif ((after <= 1.6) and  (before >= 2)) or\
        ((after>= (1/1.6)) and (before <= 0.5)):
        if (limited != limited) or (type(limited) == str):
            return 'Delete1'
        elif row[after_cat] == '0 or 1 rep':
            after_max = row[after_reps].notnull()
            exp_name = row[after_reps][after_max].index.to_list()[0]
            exp_name = exp_name.rsplit(sep = '_', maxsplit = 1)[1]
            before_value = [col for col in before_reps if exp_name in col]
            if row[before_value].item() != row[before_value].item():
                return 'Delete1'
            elif (row[before_value].item()/before >= 2) or (row[before_value].item()/before <= 0.5):
                return 'Check me, before value that corresponds to after different from med'
            else:
                return row[after_cat]
        else: 
            return row[after_cat]

    else:
        return row[after_cat]

def qc_prep(df):
    print(df[before_med])
    before = pd.to_numeric(df[before_med], errors = 'coerce')
    limited = pd.to_numeric(df[lim_med], errors = 'coerce')
    after = pd.to_numeric(df[after_med], errors = 'coerce')
    df['all/limited'] = before/limited
    all_after = before/after
    after_after = limited/after
    #can't be 2x different and one changing but the other not
    #ratio of before LPP/CTRL to after LPP/CTRL can't change 2x in one but not other
    cond0 = (~df[after_cat].str.contains('CV too high', na = False))
    cond1 = (df['all/limited'] <= 0.5) & (limited >=2) & (before<2) & cond0
    cond2 = (df['all/limited'] <= 0.5) & (before <= 0.5) & (limited >0.5) & cond0
    cond3 = (df['all/limited'] >= 2) & (before >=2) & (limited<2)& cond0
    cond4 = (df['all/limited'] >= 2) & (limited <= 0.5) & (before >0.5)& cond0
    cond5 = ((all_after >=2) | (all_after <= 0.5)) & cond0
    cond6 = (df['Mitosis_LPP/CTRL_limited_before_ratio_category'].isnull()) & cond5 
    #second 
    cond10 = (df['all/limited'] <= 0.5) | (df['all/limited'] >= 2) 
    cond11 =  (df['all/limited'] <= 0.5) & (all_after <2) & (after_after >2) & cond0
    cond12 = (df['all/limited'] >= 2 ) & (all_after >0.5) & (after_after <0.5)& cond0
    cond13 =  (df['all/limited'] <= 0.5) & (all_after <0.5) & (after_after >0.5) & cond0
    cond14 =  (df['all/limited'] >= 2) & (all_after >2) & (after_after <2)& cond0
    cond15 = ~((before >2) & (limited >2))
    # cond9 = limited.isnull()
    #don't look at after conditions if too variable
    qc_labels = df.apply(qc_checks, axis = 1)
    df.loc[:, after_cat] = qc_labels
    # df.loc[:, after_cat[1]] = qc_labels
    df.loc[(cond10) & cond0 & (cond15) , after_cat] = 'two methods different, check me!'
    df.loc[(cond11)|(cond12)|(cond13)|(cond14), after_cat] = 'Two methods diff'
    df.loc[(cond1)|(cond2)|(cond3)|(cond4)|(cond6), after_cat] = 'Delete, two methods too different (not CV too high)'
    return df

def qc_stuff():
    df1, phospho_col = phospho_merge()
    comb_median_cols = sorted([col for col in df1.columns if 'ratio_combined_median' in col])
    exp_median_cols = sorted([col for col in df1.columns if 'ratio_Exp' in col])
    CV_cols = sorted([col for col in df1.columns if 'ratio_combined_CV' in col] )
    no_cols = sorted([col for col in df1.columns if 'ratio_combined_no' in col])
    category_cols = sorted([col for col in df1.columns if 'category' in col])
    #median for exp_groups
    all_reps = sorted([col for col in df1.columns if 'rep' in col])
    keep_protein_cols = ['unique','accession', 'description', 'protein', 'sequence', 'accession_res','gene_res']
    #reorder
    df2 = df1[keep_protein_cols + phospho_col + category_cols + comb_median_cols +exp_median_cols + CV_cols + no_cols  ].copy()

    sdf = df1[keep_protein_cols + phospho_col + category_cols + comb_median_cols +exp_median_cols + CV_cols + no_cols].copy()
    sdf.loc[sdf[before_no].isnull(), before_med] = 'NQ'
    sdf.loc[(sdf[before_no] < 2), before_med] = 'NQ'
    sdf.loc[(sdf[before_cat].str.contains('CV too high')) & (~sdf[before_cat].isnull()), before_med] = 'NQ'
    
    cv_med_cols = [col for col in comb_median_cols if 'protein' not in col]
    sdf_cols = [col for col in category_cols if before_cat not in col]
    for col in sdf_cols:
        print(col)
        name = col.rsplit(sep = '_', maxsplit = 2)[0]
        med_col = [col for col in cv_med_cols if name in col][0]
        cat_col = [col for col in sdf_cols if name in col][0]
        no_col = [col for col in no_cols if name in col][0]
        sdf.loc[(sdf[cat_col].str.contains('CV too high', na = True)), med_col] = 'NQ'
        sdf.loc[(sdf[cat_col].str.contains('CV too high', na = True)), med_col] = 'NQ'
        sdf.loc[(sdf[no_col]<2) & (sdf[before_med] == 'NQ'), med_col] = 'NQ'
    df3 = df2.loc[df2[before_no] >=2]
    df3 = df3.loc[~df3['Mitosis_LPP/CTRL_before_ratio_category'].str.contains('CV too high')]
    #get rid of CV too high
    print(df3['cysteine_Asynch/Mitosis_ratio_category'])
    # df3['cysteine_Asynch/Mitosis_ratio_category'] = df3.apply(basal_cv, axis = 1)
    for col in category_cols:
        print(col)
        name = col.rsplit(sep = '_', maxsplit = 2)[0]
        med_col = [col for col in cv_med_cols if name in col][0]
        cat_col = [col for col in category_cols if name in col][0]
        df3.loc[(df3[cat_col].str.contains('CV too high', na = True)), med_col] = np.nan
        df3.loc[~ (df3[cat_col].str.contains('CV too high', na = True)), med_col] = df3[med_col]
    df3 = qc_prep(df3)
    sdf = qc_prep(sdf)
    df3.loc[(df3[after_cat].str.contains('Delete', na = False)), after_med] = np.nan
    sdf.loc[(sdf[after_cat].str.contains('Delete', na = False)), after_med] = 'NQ'
    return sdf, df3

def main():
    sdf, df3 = qc_stuff()
    # df2.loc[(cond5)|(cond6)|(cond7)|(cond8), adapted_cat_cols] = 'maybe_delete?'
    # df_dups = pd.concat(g for _, g in df2.groupby("sequence") if len(g) > 1)
    # df_dups1 = pd.concat(g for _, g in df2.groupby("gene_res") if len(g) > 1)
    # df_dups2 = pd.concat([df2, df_dups, df_dups1]).drop_duplicates(keep = False)
    # df_dups2 = df_dups2.loc[df_dups2.unique != 'U']
    #other filter: if LPP_after is changing but LPP_before isn't, check if the corresponding LPP_before rep is very different from median values
    # df3.to_csv('{}/total_combined/{}_total_combined.csv'.format(dir,date), index = False)
    print('Writing to excel1')
    with pd.ExcelWriter('{}/total_combined/{}_fig2_fig4_combined1.xlsx'.format(dir,output_date), engine='xlsxwriter', options={'strings_to_numbers': True}) as writer:  # doctest: +SKIP
        df3.to_excel(writer, sheet_name = 'filtered', index = False, float_format='%.2f')
        sdf.to_excel(writer, sheet_name = 'supplemental', index = False, float_format='%.2f')
        # df4.to_excel(writer, sheet_name = 'lost from QC8', index = False, float_format='%.2f')
        # df_dups.to_excel(writer, sheet_name = 'all_nonunique', index = False, float_format='%.2f')
        # df_dups1.to_excel(writer, sheet_name = 'gene_nonunique', index = False, float_format='%.2f')
        # df_dups2.to_excel(writer, sheet_name = 'other nonunique', index = False, float_format='%.2f')
         
main()


