#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os
from protein_tmt_census_cleanup import read_census

#main directory, where script is run from
dir = os.getcwd()

protein_cols = ['accession','description', 'protein', 'accession_res','sequence','gene_res','unique']
#set a target cutoff defined in table
target_cutoff = 2
date = '20210830'

def calculate_stats(group_dataset, inputs):
    """Calculate stats, filter. Function variable is the individual exp (i.e. Exp456)"""
    dataset_name = group_dataset.split('_')[-1]
    df = read_census(dataset_name, inputs)
    #going to rename channel columns to include dataset and group name i.e. c1 : Mitosis_CTRL1_rep_Exp420_c1
    channels = [col for col in inputs.columns if col.startswith('c')]
    groups = inputs[channels].values[0].tolist()
    new_list = [j + '_rep_' + dataset_name + '_' + i for i, j in zip(channels, groups)]
    df.rename(columns=dict(zip(channels, new_list)), inplace = True)
    #get rid of channels that aren't necessary
    drop_channels = [col for col in df.columns if col.startswith('c')]
    df.drop(drop_channels, axis = 1, inplace = True)
    group_names = [inputs.num.item(), inputs.denom.item()]
    #find the sum/sd/mean/CV of ctrl channels
    #filter by CTRL 1+2 channel
    rep_cols = [col for col in df.columns if'_rep' in col]
    df['total_sum'] = df[rep_cols].sum(axis =1)
    #require total sum of 6,000 for 6plex and 10,000 for 10plex
    sum_control = 1000 * len(channels)
    df = df.loc[df['total_sum'] >= sum_control]
    #prep for QC filtering on each experimental group, perform QC filtering
    for group in group_names:
        print('Calculating stats for ' + group)
        col = [col for col in df.columns if group + '_rep' in col]
        sum_value = df[col].sum(axis = 1)
        sd_value = df[col].std(axis = 1)
        mean_value = df[col].mean(axis = 1)
        df[group+'_median'] = df[col].median(axis = 1)
        df[group+'_CV'] = sd_value/mean_value
        df[group+'_bool'] =  (df[group+'_CV'] <= 0.5)
    df_old = df.copy()

    bool_cols = [col for col in df.columns if '_bool' in col]
    df = df.loc[df[bool_cols].any(bool_only = True, axis = 1)]
    return df

def aggregate_stats(group_dataset, inputs):
    df = calculate_stats(group_dataset, inputs)
    # See what you've lost!
    # df_removed = pd.concat([df,df_old]).drop_duplicates(keep=False)
    #get all columns with rep or median in title, then all the LPP1 cols within that subset
    rep_cols = [col for col in df.columns if '_rep' in col]
    all_medians = [col for col in df.columns if '_median' in col]
    rep_median_cols = rep_cols + all_medians 
    denom_mean = [col for col in all_medians if inputs.denom.item() in col][0]
    #subset data if denom_data value is greater than 0, then separate when it is 0 so we don't divide by 0
    exp_df = df.loc[df[denom_mean] > 0].copy()
    no_exp_data = pd.concat([df,exp_df]).drop_duplicates(keep=False)
    #divide by denominator (experimental group) unless the experimental group is 0, then just divide by 0.01. We will cap ratios later. Concat these two groups back together
    exp_df[rep_median_cols] = exp_df[rep_median_cols].div(exp_df[denom_mean], axis=0)
    no_exp_data[rep_median_cols] = no_exp_data[rep_median_cols].div(0.01).copy()
    df1 = pd.concat([exp_df, no_exp_data], sort = True)
    #Cap ratios at 20 or 0.05
    df1[rep_median_cols] = df1[rep_median_cols].apply(lambda x: [y if y <= 20 else 20 for y in x])
    df1[rep_median_cols] = df1[rep_median_cols].apply(lambda x: [y if y >= 0.05 else 0.05 for y in x])
    #make unique identifiers combining unique vs redundant + accession (this can be different for missed cleavage, etc)
    dataset_name = group_dataset.split('_')[-1]
    df1.reset_index(inplace = True)
    df1['all_pep'] = df1['tryptic_position'].astype(str) + '-' + df1['tryp_end'].astype(str)+ ':' + df1.sequence
    #Take first value for protein, description, get median and std for Asynch/Mitosis_ratio for a given protein, get the number of peptides that start at different positions for a given protein, might include overlapping peptides
    #Use only with at least 2 unique tryptic positions
    seq_list = df1.groupby('accession')['all_pep'].apply(set).apply(list)
    seq_list.replace("\'", "", inplace = True)
    seq_list.name = 'all_peptides'
    nonunique_df = df1.copy()
    df1 = df1.loc[df1['unique'] == 'U']
    identifiers = df1.groupby('accession').first()
    identifiers = identifiers[['protein','description']]
    data = df1.groupby('accession').median()
    sum_data = df1.groupby(['accession'])['tryptic_position'].nunique()
    sum_data.name = 'no_unique'
    agg_df = identifiers.merge(data, left_index = True, right_index = True)
    agg_df = agg_df.merge(sum_data, left_index = True, right_index = True)
    agg_df = agg_df.merge(seq_list, left_index = True, right_index = True)
    agg_df = agg_df.drop('tryptic_position', axis = 1)
    agg_df.reset_index(inplace = True)
    agg_df1 = agg_df.loc[agg_df.no_unique >=2]
    return agg_df1

def main():
    a = input('Which group from your inputs table would you like to run? If you would like to run all of them, write \"all\": ')
    expt_table = pd.read_csv('{}/inputs/protein_inputs.csv'.format(dir), index_col = None)
    expt_table.index = expt_table.group + '_' + expt_table.dataset
    if a == 'all':
        for group_dataset in expt_table.index:
            print('Working on ' + str(group_dataset))
            inputs = expt_table.loc[expt_table.index == group_dataset].dropna(axis = 1)
            #keep the experimental group names to remember for later
            short_df = aggregate_stats(group_dataset, inputs)
            short_df.to_csv('{}/0_scripts/replicates/{}_protein_{}.csv'.format(dir, date,group_dataset), index = False)
            # with pd.ExcelWriter('{}/basal/replicates/{}_{}.xlsx'.format(dir,date,name)) as writer:  # doctest: +SKIP
            #     short_df.to_excel(writer, sheet_name = 'filtered', index = False, float_format='%.2f')
            #     df.to_excel(writer, sheet_name = 'all', index = False, float_format = '%.2f')
    elif a not in list(expt_table.group):
        print('Error, group not listed. Make sure to copy the group from inputs table exactly')
    else:
        expt_table = expt_table.loc[expt_table.group == a]
        print('These are your inputs: ', expt_table)
        for group_dataset in expt_table.index:
            print('Working on ' + str(group_dataset))
            inputs = expt_table.loc[expt_table.index == group_dataset].dropna(axis = 1)
            #keep the experimental group names to remember for later
            short_df = aggregate_stats(group_dataset, inputs)
            short_df.to_csv('{}/0_scripts/replicates/{}_protein_{}.csv'.format(dir, date,group_dataset), index = False)

main()