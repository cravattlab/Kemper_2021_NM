#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np 
import pandas as pd
import os
from peptide_tmt_census_cleanup import read_census
# from qc_cysABPP_invert_rep import qc_filtering

#main directory, where script is run from
dir = os.getcwd()

protein_cols = ['accession','description', 'protein', 'accession_res','sequence','gene_res','unique']
#set a target cutoff defined in table
target_cutoff = 2
date = '20210830'


def calculate_stats(group_dataset, no_missed_sites, inputs, phosphopep):
    """Calculate stats, filter. Function variable is the individual exp (i.e. Exp456)"""
    dataset_name = group_dataset.split('_')[-1]
    df = read_census(dataset_name, no_missed_sites, inputs, phosphopep)
    #going to rename channel columns to include dataset and group name i.e. c1 : Mitosis_CTRL1_rep_Exp420_c1
    channels = [col for col in inputs.columns if col.startswith('c')]
    groups = inputs[channels].values[0].tolist()
    new_list = [j + '_rep_' + dataset_name + '_' + i for i, j in zip(channels, groups)]
    df.rename(columns=dict(zip(channels, new_list)), inplace = True)
    #get rid of channels that aren't necessary
    drop_channels = [col for col in df.columns if col.startswith('c')]
    df.drop(drop_channels, axis = 1, inplace = True)
    group_names = [inputs.num.item(), inputs.denom.item()]
    #for each group in a given dataset, extract group, sum, median, sd, CV across replicates
    df.set_index(protein_cols, inplace = True)
    for group in group_names:
        print('Calculating stats for ' + group)
        col = [col for col in df.columns if group + '_rep' in col]
        sum_value = df[col].sum(axis = 1)
        sd_value = df[col].std(axis = 1)
        df[group+'_mean'] = df[col].mean(axis = 1)
        df[group+'_median'] = df[col].median(axis = 1)
        df[group+'_CV'] = sd_value/df[group+'_mean']
        rep_no = len([col for col in df.columns if group + '_rep' in col])
        #if an group has a ctrl_sum <= 5,000 x number of channels and CV <= 0.5
        df[group+'_bool'] =  (sum_value >= (rep_no * 5000)) & (df[group+'_CV'] <= 0.5)
    df_old = df.copy()
    #at least one of the groups must have both sum over 5,000 and a CV of less than 0.5
    bool_cols = [col for col in df.columns if '_bool' in col]
    df = df.loc[df[bool_cols].any(bool_only = True, axis = 1)]
    print('Kept ' + str(len(df)) + ' peptide entries out of an original ' + str(len(df_old)))
    return df

def aggregate_stats(group_dataset, no_missed_sites, inputs, phosphopep):
    df = calculate_stats(group_dataset, no_missed_sites, inputs, phosphopep)
    #get all columns with rep or median in title, then all the LPP1 cols within that subset
    rep_cols = [col for col in df.columns if '_rep' in col]
    all_medians = [col for col in df.columns if '_median' in col]
    all_means = [col for col in df.columns if '_mean' in col]
    rep_median_cols = rep_cols + all_medians + all_means
    denom_median = [col for col in all_medians if inputs.denom.item() in col][0]
    num_median = [col for col in all_medians if inputs.num.item() in col][0]
    #subset data if denom_data value is greater than 0, then separate when it is 0 so we don't divide by 0
    exp_df = df.loc[df[denom_median] > 0].copy()
    no_exp_data = pd.concat([df,exp_df]).drop_duplicates(keep=False)
    #divide by denominator (experimental group) unless the experimental group is 0, then just divide by 0.01. We will cap ratios later. Concat these two groups back together
    exp_df[rep_median_cols] = exp_df[rep_median_cols].div(exp_df[denom_median], axis=0)
    no_exp_data[rep_median_cols] = no_exp_data[rep_median_cols].div(0.01).copy()
    df1 = pd.concat([exp_df, no_exp_data], sort = True)
    #Cap ratios at 20 or 0.05
    df1[rep_median_cols] = df1[rep_median_cols].apply(lambda x: [y if y <= 20 else 20 for y in x])
    df1[rep_median_cols] = df1[rep_median_cols].apply(lambda x: [y if y >= 0.05 else 0.05 for y in x])
    #make unique identifiers combining unique vs redundant + accession (this can be different for missed cleavage, etc)
    df1.reset_index(inplace = True)
    df1['unique'] = df1['unique'].replace(np.nan, 'R')
    # df1['accession_res'] = df1['unique'] + '_' + df1['accession_res']
    bool_cols = [col for col in df1.columns if '_bool' in col]
    df1.drop(bool_cols, axis =1, inplace = True)
    #get the smallest sequence for each accession_res, and get only the protein_cols 
    identifiers = df1.loc[df1['sequence'].str.len().groupby(df1['accession_res']).idxmin()]
    identifiers = identifiers[protein_cols]
    #for all other columns get the median value
    data = df1.groupby('accession_res', as_index = False).median()
    df1 = identifiers.merge(data, on = 'accession_res')
    ##if the ratio is less than 1, we will flip the ratio then take the CV
    # df1['small_ratio_bool'] = (df1[num_median] < 1)
    return df1

# def qc_data(group_dataset, inputs):
#     df1 = aggregate_stats(group_dataset, inputs)
#     unfiltered_df = df1.copy()
#     # #For ratios greater than 4, ctrl1 values are small and the QC filtering fucks it up
#     rep_cols = sorted([col for col in df1.columns if '_rep' in col])
#     #keep only rep columns to do QC filtering
#     short_df = df1[protein_cols + rep_cols]
#     short_df = short_df.set_index(protein_cols)
#     #prep for QC filtering on each experimental group, perform QC filtering
#     group_names = [inputs.num.item(), inputs.denom.item()]
#     for group in sorted(set(group_names)):
#         print(group)
#         col = [col for col in df1.columns if group in col]
#         test_df = df1[protein_cols + col]
#         test_df = test_df.set_index(protein_cols)
#         mean_mr = test_df.apply(qc_filtering, args=(group,), axis = 1)
#         short_df = short_df.merge(mean_mr, how = 'outer', left_index = True, right_index = True)
#     short_df.reset_index(inplace = True)
#     combined_df = short_df.dropna(axis=1, how='all')
#     #recalculate median, not counting if CV is too high in denominator
#     numerator_median = inputs.num.item() + '_median'
#     denominator_median = inputs.denom.item() + '_median'
#     combined_df['Need_to_delete_median?'] =  (combined_df[numerator_median].isnull()) | (combined_df[denominator_median].isnull())
#     combined_df[numerator_median] = combined_df[numerator_median]/combined_df[denominator_median]

#     # short_df.drop('small_div_col', axis = 1, inplace = True)
#     return combined_df

def no_missed_tryptic(mod):
    mod = str(mod)
    if '79.96633' in mod:
        return 2
    else:
        return 1

def main():
#for a dataset: perform functions defined above
    a = input('Which group from your inputs table would you like to run? If you would like to run all of them, write \"all\": ')
    phosphopep = input('Is this to find cysteines with phosphosites on the same tryptic? If yes, write \"yes\", else write \"no\": ')
    expt_table = pd.read_csv('{}/inputs/peptide_inputs.csv'.format(dir), index_col = None)
    expt_table.index = expt_table.group + '_' + expt_table.dataset
    if a == 'all':
        for group_dataset in expt_table.index:
            print('Working on ' + str(group_dataset))
            inputs = expt_table.loc[expt_table.index == group_dataset].dropna(axis = 1)
            #keep the experimental group names to remember for later
            mod = inputs.modification.item()
            no_missed_sites = no_missed_tryptic(mod)
            print('Checking for ' + str(no_missed_sites) + ' missed cleavage sites')
            df = aggregate_stats(group_dataset, no_missed_sites, inputs, phosphopep)
            if phosphopep == 'yes':
                df.to_csv('{}/0_scripts/replicates/{}_phospho_peptide_{}.csv'.format(dir, date,group_dataset), index = False)
            else:
                df.to_csv('{}/0_scripts/replicates/{}_peptide_{}.csv'.format(dir, date,group_dataset), index = False)

    elif a not in list(expt_table.group):
        print('Error, group not listed. Make sure to copy the group from inputs table exactly')
    else:
        expt_table = expt_table.loc[expt_table.group == a]
        print('These are your inputs: ', expt_table)
        for group_dataset in expt_table.index:
            print('Working on ' + str(group_dataset))
            inputs = expt_table.loc[expt_table.index == group_dataset].dropna(axis = 1)
            mod = inputs.modification.item()
            print(mod)
            no_missed_sites = no_missed_tryptic(mod)
            print('Checking for ' + str(no_missed_sites) + ' missed cleavage sites')
            #keep the experimental group names to remember for later
            df = aggregate_stats(group_dataset, no_missed_sites, inputs, phosphopep)
            if phosphopep == 'yes':
                df.to_csv('{}/0_scripts/replicates/{}_phospho_peptide_{}.csv'.format(dir, date,group_dataset), index = False)
            else:
                df.to_csv('{}/0_scripts/replicates/{}_peptide_{}.csv'.format(dir, date,group_dataset), index = False)

main()



