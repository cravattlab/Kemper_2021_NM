#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from functools import reduce
import os


#main directory, where script is run from
dir = os.getcwd()

protein_cols = ['accession','description', 'protein']
#set a target cutoff defined in table
target_cutoff = 2
date = '20210830'

def combine_reps(compile_df, df, inputs):
    """combining our individual datasets"""
    # add exp name to columns of data columns
    dataset_name = inputs.dataset.item()
    df.columns = ['{}{}'.format(c, '' if c in protein_cols else '_' + dataset_name) for c in df.columns]
    compile_df.append(df)
    combined_df = reduce(lambda  left,right: pd.merge(left,right,on=protein_cols, how='outer'), compile_df)
    return combined_df

# def sequence_combine(row):
#      reps = [col for col in row.index if 'sequence' in col]
#      sequences = row[reps].dropna()
#      shortest_sequence =  min(sequences, key=len)
#      return pd.Series(shortest_sequence)

def calculate_agg_stats(df, expt_table):
    #we just need the inputs for the num and denom name, so taking first row should be fine
    inputs = expt_table.iloc[0,:]
    #get the shortest sequence from all TMT reps to represent in final data
    #find median, CV, count reps for the ratio of numerator/denominator 
    numerator = inputs.num
    denominator = inputs.denom
    median_cols = [col for col in df.columns if numerator + '_median' in col]
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
    return df

def cleanup_df(compile_df, expt_table):
    all = calculate_agg_stats(compile_df, expt_table)
    aggregate = all.set_index(protein_cols)
    # cols for TMT reps
    rep_median_cols = sorted([col for col in aggregate.columns if 'median_Exp' in col])
    inputs = expt_table.iloc[0,:]
    denominator = inputs.denom
    numerator = inputs.num
    # cols for TMT reps, only keep numerator because denom will just be 1
    rep_median_cols = sorted([col for col in aggregate.columns if numerator + '_median_Exp' in col])
    #rename all TMT replicates with rep1, rep2, etc
    rep_median_names = pd.Series(rep_median_cols).str.split('Exp', expand = True)[0] + 'TMT_rep_'
    rep_median_no = [str(x) for x in list(range(1,len(rep_median_cols)+1))]
    #replace just numerator in column names as ratio names
    new_rep_median_cols = [i + j for i, j in zip(list(rep_median_names), rep_median_no)]
    new_rep_median_cols = [x.replace(numerator, numerator+'/'+denominator) for x in new_rep_median_cols]
    rename_median_dict = dict(zip(rep_median_cols, new_rep_median_cols))
    aggregate.rename(columns = rename_median_dict, inplace = True)

    #median cols for combined ratio cols
    comb_ratio_cols = sorted([col for col in aggregate.columns if 'ratio_combined_median' in col])
    comb_ratio_cv_cols = sorted([col for col in aggregate.columns if 'ratio_combined_CV' in col])
    comb_ratio_no_cols = sorted([col for col in aggregate.columns if 'ratio_combined_no' in col])
    category_cols = [col for col in aggregate.columns if 'ratio_category' in col]
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
    all_df = aggregate[category_cols + comb_ratio_cols + new_rep_median_cols + comb_ratio_cv_cols + comb_ratio_no_cols + all_reps].copy()
    all_df.reset_index(inplace = True)
    return all_df

def main():
#for a dataset: perform functions defined above

    a = input('Which group from your inputs table would you like to run? If you would like to run all of them, write \"all\": ')
    expt_table = pd.read_csv('{}/inputs/protein_inputs.csv'.format(dir), index_col = None)
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
                df = pd.read_csv('{}/0_scripts/replicates/{}_protein_{}.csv'.format(dir, date,group_dataset))
                combined_df = combine_reps(compile_df, df, inputs)
            df1 = cleanup_df(combined_df, expt_table_subset)
            df1.to_csv('{}/1_scripts/combined/{}_combined_protein_{}.csv'.format(dir, date,group_name), index = False)
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
            df = pd.read_csv('{}/0_scripts/replicates/{}_protein_{}.csv'.format(dir, date,group_dataset))
            combined_df = combine_reps(compile_df, df, inputs)
        df1 = cleanup_df(combined_df, expt_table) 
        df1.to_csv('{}/1_scripts/combined/{}_combined_protein_{}.csv'.format(dir, date,group_name), index = False)

main()