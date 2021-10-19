#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
from functools import reduce
from qc_filters_invert_combine1 import qc_filtering

#main directory, where script is run from
dir = os.getcwd()

protein_cols = ['accession','description', 'protein', 'accession_res','gene_res']
#set a target cutoff defined in table
target_cutoff = 2
date = '20210830'

def combine_reps(compile_df, df, inputs):
    """combining our individual datasets"""
    # add exp name to columns of data columns
    dataset_name = inputs.dataset.item()
    rep = inputs.rep.item()
    numerator = inputs.num.item()
    denominator = inputs.denom.item()
    df.columns = ['{}{}'.format(c, '' if c in protein_cols else '_' + dataset_name) for c in df.columns]
    med_col = [col for col in df.columns if numerator + '_median' in col][0]
    # #rename all TMT replicates with rep1, rep2, etc
    new_med_col = med_col.split('Exp')[0] + 'TMT_' + rep
    new_med_col1 = new_med_col.replace(numerator, numerator+'/'+denominator)
    df.rename(columns = {med_col : new_med_col1}, inplace = True)
    compile_df.append(df)
    combined_df = reduce(lambda  left,right: pd.merge(left,right,on=protein_cols, how='outer'), compile_df)
    return combined_df

def sequence_combine(row):
     reps = [col for col in row.index if 'sequence' in col]
     sequences = row[reps].dropna()
     shortest_sequence =  min(sequences, key=len)
     return pd.Series(shortest_sequence)

def calculate_agg_stats(df, expt_table):
    #we just need the inputs for the num and denom name, so taking first row should be fine
    inputs = expt_table.iloc[0,:]
    #get the shortest sequence from all TMT reps to represent in final data
    sequence_cols = [col for col in df.columns if 'sequence' in col]
    df['sequence'] = df[sequence_cols].apply(sequence_combine, axis = 1)
    df.drop(sequence_cols, axis = 1 , inplace = True)
    #find median, CV, count reps for the ratio of numerator/denominator 
    numerator = inputs.num
    denominator = inputs.denom
    median_cols = [col for col in df.columns if numerator+'/'+denominator + '_median' in col]
    print(median_cols)
    #store the ratio name in variable
    ratio_name = numerator+ '/'+denominator
    #make LPP3 into ratio 
    #df[mitosis_ratio+'_ratio_combined_median'] = df['Mitosis_LPP3_combined_median']/df['Mitosis_CTRL3_combined_median']
    df[ratio_name+'_ratio_combined_max'] = df[median_cols].max(axis = 1)
    df[ratio_name+'_ratio_combined_min'] = df[median_cols].min(axis = 1)
    old_df = df.copy()
    df[ratio_name+'_ratio_combined_median'] = df[median_cols].median(axis = 1)
    df[ratio_name+'_ratio_combined_no'] = df[median_cols].count(axis = 1)
    mean_value = df[median_cols].mean(axis = 1)
    sd_value = df[median_cols].std(axis = 1)
    df[ratio_name+'_ratio_combined_CV'] = sd_value/mean_value
    ratio_col = [col for col in df.columns if ratio_name in col]
    rep_col = [col for col in df.columns if numerator + '_median' in col]
    all_cols = ratio_col +rep_col
    qc_df = df[all_cols]
    mean_mr = qc_df.apply(qc_filtering, axis = 1)
    mean_mr.columns = mean_mr.columns.str.replace("qc", str(ratio_name))
    df = old_df.merge(mean_mr, how = 'outer', left_index = True, right_index = True)
    return df

def cleanup_df(compile_df, expt_table):
    all = calculate_agg_stats(compile_df, expt_table)
    protein_cols = ['accession','description', 'protein', 'accession_res','gene_res','sequence']
    aggregate = all.set_index(protein_cols)
    inputs = expt_table.iloc[0,:]
    denominator = inputs.denom
    numerator = inputs.num
    
    #median cols for combined ratio cols
    comb_ratio_cols = sorted([col for col in aggregate.columns if 'ratio_combined_median' in col])
    comb_ratio_cv_cols = sorted([col for col in aggregate.columns if 'ratio_combined_CV' in col])
    comb_ratio_no_cols = sorted([col for col in aggregate.columns if 'ratio_combined_no' in col])
    category_cols = [col for col in aggregate.columns if 'ratio_category' in col]
    #median for TMT experiments
    median_cols = [col for col in aggregate.columns if 'median_TMT_rep' in col]
    #rename all channel replicates with rep1, rep2, etc
    num_reps = [col for col in aggregate.columns if numerator + '_rep' in col]
    denom_reps = [col for col in aggregate.columns if denominator + '_rep' in col]
    num_reps_no = [str(x) for x in list(range(1,len(num_reps)+1))]
    denom_reps_no = [str(x) for x in list(range(1,len(denom_reps)+1))]
    num_names = pd.Series(num_reps).str.split('_rep_', expand = True)[0] + '_channel_rep_'
    denom_names = pd.Series(denom_reps).str.split('_rep_', expand = True)[0] + '_channel_rep_'
    new_num_reps = [i + j for i, j in zip(list(num_names), num_reps_no)]
    new_denom_reps = [i + j for i, j in zip(list(denom_names), denom_reps_no)]
    num_rep_dict = dict(zip(num_reps, new_num_reps))
    denom_rep_dict = dict(zip(denom_reps, new_denom_reps))
    aggregate.rename(columns = num_rep_dict, inplace = True)
    aggregate.rename(columns = denom_rep_dict, inplace = True)
    all_reps = new_num_reps + new_denom_reps
    categories = sorted([col for col in aggregate.columns if 'categor' in col])
    #select the columns we want in the order we want
    #change for repxx??
    all_df = aggregate[category_cols + comb_ratio_cols + median_cols + comb_ratio_cv_cols + comb_ratio_no_cols + all_reps].copy()
    all_df.reset_index(inplace = True)
    return all_df

def main():
#for a dataset: perform functions defined above

    a = input('Which group from your inputs table would you like to run? If you would like to run all of them, write \"all\": ')
    expt_table = pd.read_csv('{}/inputs/peptide_inputs.csv'.format(dir), index_col = None)
    expt_table.index = expt_table.group + '_' + expt_table.dataset
    if a == 'all':
        for exp_group in list(set(expt_table.group)):
            compile_df = []
            expt_table_subset = expt_table.loc[expt_table.group == exp_group]
            group_name = expt_table_subset['group'][0]
            print('Working on ' + group_name)
            for group_dataset in expt_table_subset.index:
                print('Combining ' + str(group_dataset))
                inputs = expt_table_subset.loc[expt_table_subset.index == group_dataset].dropna(axis = 1)
                df = pd.read_csv('{}/0_scripts/replicates/{}_peptide_{}.csv'.format(dir, date,group_dataset))
                combined_df = combine_reps(compile_df, df, inputs)

            df1 = cleanup_df(combined_df, expt_table_subset)
                
            df1.to_csv('{}/1_scripts/combined/{}_combined_{}.csv'.format(dir, date,group_name), index = False)
    elif a not in list(expt_table.group):
        print('Error, group not listed. Make sure to copy the group from inputs table exactly')
    else:
        compile_df = []
        expt_table = expt_table.loc[expt_table.group == a]
        group_name = expt_table['group'][0]
        # print('These are your inputs: ', expt_table)
        print('Working on ' + group_name)
        for group_dataset in expt_table.index:
            print('Combining ' + str(group_dataset))
            inputs = expt_table.loc[expt_table.index == group_dataset].dropna(axis = 1)
            #keep the experimental group names to remember for later
            df = pd.read_csv('{}/0_scripts/replicates/{}_peptide_{}.csv'.format(dir, date,group_dataset))
            combined_df = combine_reps(compile_df, df, inputs)
        df1 = cleanup_df(combined_df, expt_table) 
        df1.to_csv('{}/1_scripts/combined/{}_combined_{}.csv'.format(dir, date,group_name), index = False)

main()