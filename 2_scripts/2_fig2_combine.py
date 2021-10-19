#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 15:03:08 2021

@author: estherkemper
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

protein_cols = ['accession','description', 'protein', 'accession_res','gene_res']
lim_med = 'Mitosis_LPP(-,-)/(+,-)_limited_combined_median'
lim_cat = 'Mitosis_LPP(-,-)/(+,-)_limited_category'
original_med = 'Mitosis_LPP(-,-)/(+,-)_combined_median'
original_no = 'Mitosis_LPP(-,-)/(+,-)_combined_no'
adapted_no = 'Mitosis_LPP(-,+)/(+,+)_combined_no'
original_cat = 'Mitosis_LPP(-,-)/(+,-)_category'
adapted_med = 'Mitosis_LPP(-,+)/(+,+)_combined_median'
adapted_cat = 'Mitosis_LPP(-,+)/(+,+)_category'

expt_table = pd.read_csv('{}/inputs/protein_inputs.csv'.format(dir), index_col = None)
lpp = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')

numerator = 'Mitosis_CTRL1'
adapted_num = 'Mitosis_CTRL3'
basal_denom = 'Asynch_CTRL2'
original_denom = 'Mitosis_LPP1'
adapted_denom = 'Mitosis_LPP3'

def combine():
    print('Reading input datasets')
    original_df = pd.read_csv('{}/1_scripts/combined/{}_combined_all_original.csv'.format(dir,date))
    lim_original_df = pd.read_csv('{}/1_scripts/combined/{}_combined_limited_original.csv'.format(dir,date))
    adapted_df = pd.read_csv('{}/1_scripts/combined/{}_combined_adapted.csv'.format(dir,date))
    basal_df = pd.read_csv('{}/1_scripts/combined/{}_combined_basal.csv'.format(dir,date))
    protein = pd.read_csv('{}/1_scripts/combined/{}_combined_protein_whole_proteome.csv'.format(dir,date))
    
    protein.columns = ['protein_' + x if x not in protein_cols else x for x in protein.columns]
    basal_df.columns = ['cysteine_' + x if x not in protein_cols else x for x in basal_df.columns]
    basal_df.columns = basal_df.columns.str.replace(numerator+'/'+basal_denom, 'Mitosis/Asynch')
    protein.columns = protein.columns.str.replace(numerator+'/'+basal_denom, 'Mitosis/Asynch')
    lim_keep_cols = [col for col in lim_original_df.columns if 'combined_median' in col]
    lim_keep_cols = lim_keep_cols + [col for col in lim_original_df.columns if 'category' in col]
    lim_original_df = lim_original_df[protein_cols + ['sequence'] + lim_keep_cols]
    lim_original_df.columns = lim_original_df.columns.str.replace(numerator + '/' + original_denom,'Mitosis_LPP(-,-)/(+,-)_limited')
    original_df.columns = original_df.columns.str.replace(numerator + '/' + original_denom,'Mitosis_LPP(-,-)/(+,-)')
    adapted_df.columns = adapted_df.columns.str.replace(adapted_num + '/' + adapted_denom,'Mitosis_LPP(-,+)/(+,+)')

    print('Merging data')
    df = original_df.merge(basal_df, how = 'outer', on = protein_cols, suffixes = [' original', ' basal'])
    df = df.merge(adapted_df, how = 'outer', on = protein_cols, suffixes = ['', ' adapted'])
    df = df.merge(lim_original_df, how = 'outer', on = protein_cols, suffixes = ['', ' original_LIMITED'])
    df = df.merge(protein, how = 'left', on = ['accession','description','protein'], suffixes = ['', ' protein'])
    rep_cols = [col for col in df.columns if '_channel_' in col]
    df.drop(rep_cols, axis = 1, inplace = True)
    df.columns = df.columns.str.replace('_ratio_','_')
    
    return df

def sequence_combine(row):
     reps = [col for col in row.index if 'sequence' in col]
     sequences = row[reps].dropna()
     shortest_sequence =  min(sequences, key=len)
     return pd.Series(shortest_sequence)
 
def phospho_merge():
    df = combine()
    df.rename(columns = {'sequence':'sequence original'}, inplace = True)
    #get rid of empty columns
    df = df.dropna(axis = 1, how = 'all')
    sequence_cols = [col for col in df.columns if 'sequence' in col]
    df['sequence'] = df[sequence_cols].apply(sequence_combine, axis = 1)
    df.drop(sequence_cols, axis = 1 , inplace = True)
    df['gene'] = df.gene_res.str.split(pat = '_', expand = True)[0]
    phospho_df = pd.read_csv('{}/other_datasets/final_phospho_mann_olsen.csv'.format(dir))
    keep_cols = [col for col in phospho_df.columns if '_olsen' not in col]
    phospho_df = phospho_df[keep_cols]
    phospho_col = [col for col in phospho_df.columns if 'accession' not in col]
    df1 = df.merge(phospho_df, on = 'accession', how = 'left')
    print('Finished merging')
    return df1, phospho_col

def qc_checks(row):
    #if after value exists but corresponding value in before lPP/ctrl not there, do not interpret
    adapted_reps = [col for col in row.index if '+)_median_TMT' in col]
    original_reps = [col for col in row.index if '-)_median_TMT' in col]
    original = pd.to_numeric(row[original_med], errors = 'coerce')
    limited = pd.to_numeric(row[lim_med], errors = 'coerce')
    after = pd.to_numeric(row[adapted_med], errors = 'coerce')
    #for changing after, bu tnot changing original, with original not matching or nulll
    if ((after >=2) and (original <= 1.6)) or \
        ((after <= 0.5) and (original >= (1/1.6))):  
        if (limited != limited) or (type(limited) == str):
            return 'Delete, corresponding is null'
        elif row[adapted_cat] == '0 or 1 rep':
            after_max = row[adapted_reps].notnull()
            exp_name = row[adapted_reps][after_max].index.to_list()[0]
            exp_name = exp_name.rsplit(sep = '_', maxsplit = 1)[1]
            original_value = [col for col in original_reps if exp_name in col][0]
            if pd.isna(row[original_value]) or row[original_value] == 'NQ':
                return 'Delete, changing after but not quant in original'
            elif (row[original_value]/after < 2) and (original/after >2 ):
                return 'Delete, corresponding not changing much but combined does'
            else:
                return row[adapted_cat]
        else:
            return row[adapted_cat]
    #for unchanging after but changing original
    elif ((after <= 1.6) and  (original >= 2)) or\
        ((after>= (1/1.6)) and (original <= 0.5)):
        if (limited != limited) or (type(limited) == str):
            return 'Delete, corresponding is null'
        elif row[adapted_cat] == '0 or 1 rep':
            after_max = row[adapted_reps].notnull()
            exp_name = row[adapted_reps][after_max].index.to_list()[0]
            exp_name = exp_name.rsplit(sep = '_', maxsplit = 1)[1]
            original_value = [col for col in original_reps if exp_name in col][0]
            if ((pd.isna(row[original_value]))) or (row[original_value] == 'NQ'):
                return 'Delete, corresponding is null'
            elif row[original_value]/original >= 2 or row[original_value]/original <= 0.5:
                return 'Check me, original value that corresponds to after different from med'
            else:
                return row[adapted_cat]
        else: 
            return row[adapted_cat]

    else:
        return row[adapted_cat]

def qc_prep(df):
    original = pd.to_numeric(df[original_med], errors = 'coerce')
    limited = pd.to_numeric(df[lim_med], errors = 'coerce')
    after = pd.to_numeric(df[adapted_med], errors = 'coerce')
    df['all/limited'] = original/limited
    all_after = original/after
    after_after = limited/after
    #can't be 2x different and one changing but the other not
    #ratio of original LPP/CTRL to after LPP/CTRL can't change 2x in one but not other
    cond0 = (~df[adapted_cat].str.contains('CV too high', na = False))
    cond1 = (df['all/limited'] <= 0.5) & (limited >=2) & (original<2) & cond0
    cond2 = (df['all/limited'] <= 0.5) & (original <= 0.5) & (limited >0.5) & cond0
    cond3 = (df['all/limited'] >= 2) & (original >=2) & (limited<2)& cond0
    cond4 = (df['all/limited'] >= 2) & (limited <= 0.5) & (original >0.5)& cond0
    cond5 = ((all_after >=2) | (all_after <= 0.5)) & cond0
    
    cond5a = ( df['all/limited'] >=2) & (after< 0.55) & (original > 1)
    cond5b =  ( df['all/limited'] <= 0.5) & (original< 0.55) & (after > 1)

    cond6 = (df[lim_cat].isnull()) & cond5 
    #second 
    cond10 = (df['all/limited'] <= 0.5) | (df['all/limited'] >= 2) 
    cond11 =  (df['all/limited'] <= 0.5) & (all_after <2) & (after_after >2) & cond0
    cond12 = (df['all/limited'] >= 2 ) & (all_after >0.5) & (after_after <0.5)& cond0
    cond13 =  (df['all/limited'] <= 0.5) & (all_after <0.5) & (after_after >0.5) & cond0
    cond14 =  (df['all/limited'] >= 2) & (all_after >2) & (after_after <2)& cond0
    cond15 = ~((original >2) & (limited >2))
    # cond9 = limited.isnull()
    #don't look at after conditions if too variable
    qc_labels = df.apply(qc_checks, axis = 1)
    df.loc[:, adapted_cat] = qc_labels
    df.loc[:, adapted_cat[1]] = qc_labels
    df.loc[(cond10) & cond0 & (cond15) , adapted_cat] = 'two methods different, check me!'
    df.loc[(cond11)|(cond12)|(cond13)|(cond14), adapted_cat] = 'Two methods diff'
    df.loc[(cond5a)|(cond5b), adapted_cat] = 'Delete, CV too high'
    df.loc[(cond1)|(cond2)|(cond3)|(cond4)|(cond6), adapted_cat] = 'Delete, two methods too different (not CV too high)'
    return df

def qc_stuff():
    df1, phospho_col = phospho_merge()
    
    comb_median_cols = sorted([col for col in df1.columns if '_combined_median' in col])
    exp_median_cols = sorted([col for col in df1.columns if 'median_TMT' in col])
    CV_cols = sorted([col for col in df1.columns if 'combined_CV' in col] )
    no_cols = sorted([col for col in df1.columns if 'combined_no' in col])
    category_cols = sorted([col for col in df1.columns if 'category' in col])
    keep_protein_cols = ['accession', 'description', 'protein', 'sequence', 'accession_res','gene_res']
    #reorder
    df2 = df1[keep_protein_cols + phospho_col + category_cols + comb_median_cols + CV_cols + no_cols + exp_median_cols].copy()

    sdf = df1[keep_protein_cols + phospho_col + category_cols + comb_median_cols + CV_cols + no_cols + exp_median_cols].copy()
    sdf[comb_median_cols + exp_median_cols] = np.log2(sdf[comb_median_cols + exp_median_cols])
    sdf = sdf.sort_values(original_med, ascending = False)


    sdf.loc[sdf[original_no].isnull(), original_med] = 'NQ'
    sdf.loc[(sdf[original_no] < 2), original_med] = 'NQ'
    sdf.loc[(sdf[original_cat].str.contains('CV too high')) & (~sdf[original_cat].isnull()), original_med] = 'NQ'
    
    cv_med_cols = [col for col in comb_median_cols if 'protein' not in col]
    sdf_cols = [col for col in category_cols if lim_cat not in col]
    for col in sdf_cols:
        name = col.rsplit(sep = '_', maxsplit = 1)[0]
        med_col = [col for col in cv_med_cols if name in col][0]
        cat_col = [col for col in sdf_cols if name in col][0]
        no_col = [col for col in no_cols if name in col][0]
        sdf.loc[(sdf[cat_col].str.contains('CV too high', na = True)), med_col] = 'NQ'
        sdf.loc[(sdf[cat_col].str.contains('CV too high', na = True)), med_col] = 'NQ'
        sdf.loc[(sdf[no_col]<2) & (sdf[original_med] == 'NQ'), med_col] = 'NQ'
    df3 = df2.loc[df2[original_no] >=2]
    df3 = df3.loc[~df3[original_cat].str.contains('CV too high')]
    for col in sdf_cols:
        name = col.rsplit(sep = '_', maxsplit = 1)[0]
        med_col = [col for col in cv_med_cols if name in col][0]
        cat_col = [col for col in sdf_cols if name in col][0]
        df3.loc[(df3[cat_col].str.contains('CV too high', na = True)), med_col] = np.nan
        df3.loc[~ (df3[cat_col].str.contains('CV too high', na = True)), med_col] = df3[med_col]
    df3 = qc_prep(df3)
    sdf = qc_prep(sdf)
    df3.loc[(df3[adapted_cat].str.contains('Delete', na = False)), adapted_med] = np.nan
    sdf.loc[(sdf[adapted_cat].str.contains('Delete', na = False)), adapted_med] = 'NQ'
    # sdf = sdf.loc[ (sdf[original_no].notnull()) & (sdf[adapted_no].notnull())]
    sdf[original_med] = sdf[original_med].replace('NQ', np.nan)
    sdf = sdf.sort_values(original_med, ascending = False)
    sdf[original_med] = sdf[original_med].replace(np.nan, 'NQ')
    print(sdf)
    return sdf, df3

def lpp_merge():
    sdf, df3 = qc_stuff()
    p_med = 'cysteine_' + original_med
    print(lpp.columns)
    lpp1 = lpp[['accession', 'gene_res', p_med]]
    lpp1 = lpp1.round(decimals = 1)
    lpp2 = lpp1.loc[lpp1[p_med] >= 2]
    lpp2 = pd.DataFrame(lpp2.groupby('accession')['gene_res'].apply(lambda row: ', '.join(row.astype(str)))).reset_index()
    lpp2.rename({'gene_res':'LPP-sensitive phosphosite'}, axis =1, inplace = True)
    sdf1 = sdf.merge(lpp2, on = 'accession', how = 'left')
    return sdf1, df3

def main():
    sdf, df3 = lpp_merge()
    drop_cols = [col for col in sdf.columns if 'category' in col]
    sdf.drop(drop_cols, axis =1, inplace = True)
    sdf = sdf.assign(accession_res = sdf.accession_res.str.rsplit('_', n = 1, expand = True))
    print('Writing to excel1')
    with pd.ExcelWriter('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir,date)) as writer:  # doctest: +SKIP
        df3.to_excel(writer, sheet_name = 'filtered', index = False)
        sdf.to_excel(writer, sheet_name = 'supplemental', index = False, float_format = '%.2f')
         
main()


