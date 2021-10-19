#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 8th

@author: EKK
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
plt.rcParams['pdf.fonttype'] = 42


dir = os.getcwd()
output_dir= dir + '/1_scripts/combined_isotop'
inputs = '/inputs/isotop_peptide_inputs.csv'
target_cutoff = 2
protein_cols = ['accession','description','symbol', 'protein','sequence', 'gene_res']
date = '20210830'
expt_table = pd.read_csv('{}/{}'.format(dir, inputs))


def multicombine(group):
    df = pd.DataFrame(columns=protein_cols)
    reps = expt_table.loc[expt_table.reps.str.contains(group), 'reps']
    for rep in reps:
        rep_num = rep.rsplit('_')[-1]
        replicate = pd.read_csv('{}/0_scripts/replicates_isotop/{}_{}.csv'.format(dir, date,rep))
        data_cols = [col for col in replicate.columns if col not in protein_cols]
        data_cols = data_cols + ['sequence']
        rename_dict = {col:col+'_'+rep_num for col in data_cols}
        replicate = replicate[protein_cols + data_cols]
        replicate = replicate.rename(columns = rename_dict)
        df = df.merge(replicate, on=['accession','description', 'symbol', 'protein', 'gene_res'], how='outer')
        # df = df[df.columns.drop(list(df.filter(regex='std')))]
    return df

def sequence_combine(row):
     reps = [col for col in row.index if 'sequence' in col]
     sequences = row[reps].dropna()
     shortest_sequence =  min(sequences, key=len)
     return pd.Series(shortest_sequence)

def calculate_mean_IR(row):
    headers = ['combined_median','combined_CV', 'combined_no','category']
    IR_cols = [col for col in row.index if 'IR_median' in col]
    combined_median = round(row[IR_cols].median(), 2)
    old_cv = [col for col in row.index if 'combined_CV' in col][0]
    combined_min = [col for col in row.index if '_min' in col][0]
    combined_max = [col for col in row.index if '_max' in col][0]
    no_cols = [col for col in row.index if 'combined_no' in col][0]
    if np.isnan(row[old_cv]) == True:
        return pd.Series([combined_median, row[old_cv],row[no_cols],'0 or 1 rep'], index=headers)
    elif row[old_cv] <= 0.6:
        return pd.Series([combined_median, row[old_cv],row[no_cols],'Keep'], index=headers)
    elif row[old_cv] > 0.6:
        if row[combined_min] >= 1.6:
            return pd.Series([combined_median, row[old_cv],row[no_cols],'Keep, CV high but all values changing'], index=headers)
        elif row[combined_max] <= (1/1.6):
            return pd.Series([combined_median, row[old_cv],row[no_cols],'Keep, CV high but all values changing'], index=headers)
        elif row[combined_max] < 1.6:
            if row[combined_min] > (1/1.6):
                return pd.Series([combined_median, row[old_cv],row[no_cols],'Keep, high CV but medians unchanging'], index=headers)
            else:
                return pd.Series([combined_median, row[old_cv],row[no_cols],'CV too high, do not interpret'], index=headers)
    else:
        return pd.Series([combined_median, row[old_cv],row[no_cols],'CV too high, do not interpret'], index=headers)
        # elif row[no_cols] >=4:
        #     mean_value = row[IR_cols].mean()
        #     mxid = abs(row[IR_cols] - mean_value).idxmax()
        #     new_reps = row[IR_cols].drop(mxid)
        #     new_mean = new_reps.mean()
        #     new_cv = new_reps.std()/new_mean
        #     if row[combined_min] > 1.5 or row[combined_max] < (1/1.5):
        #         return pd.Series([combined_median, row[old_cv],row[no_cols],'Keep, high CV but all above 1.5'], index=headers)
        #     elif new_cv <= 0.6:         
        #         return pd.Series([combined_median, row[old_cv],row[no_cols],'CV okay if get rid of rep farthest from mean'], index=headers)
        #     elif combined_median >=2:
        #         mnid = row[IR_cols].idxmin()
        #         new_reps = row[IR_cols].drop(mnid)
        #         new_min = new_reps.min()
        #         if new_min >= 1.8:
        #             return pd.Series([combined_median, row[old_cv],row[no_cols],'All changing if get rid of min'], index=headers)
        #         else:
        #             return pd.Series(['NQ', row[old_cv],row[no_cols],'CV too high, do not interpret!'], index=headers)
        #     elif combined_median <= 0.5:
        #         mxid = row[IR_cols].idxmax()
        #         new_reps = row[IR_cols].drop(mxid)
        #         new_max = new_reps.max()
        #         if new_max <= (1/1.8):
        #             return pd.Series([combined_median, row[old_cv],row[no_cols],'All changing if get rid of max'], index=headers)
        #         else:
        #             return pd.Series([combined_median, row[old_cv],row[no_cols],'CV too high, do not interpret!!'], index=headers)
        #     else:
        #         return pd.Series(['NQ', row[old_cv],row[no_cols],'CV too high, do not interpret!!'], index=headers)
        # else:
        #     return pd.Series(['NQ', row[old_cv],row[no_cols],'CV too high, do not interpret!!!'], index=headers)


    # # if combined_median == 20.0:
    # #     return pd.Series([combined_median,'keep; and meets target cutoff'], index=headers)
    # elif combined_median >= target_cutoff and combined_median < 20:
    #         #SD filter
    #         #report the lower value if stdev is too large
    #         if stdev>0.6*combined_median and row[IR_cols].min()<target_cutoff:
    #             print(type(abs(row[IR_cols] - combined_median)[0]))
    #             mxid = abs(row[IR_cols] - combined_median).idxmax()
    #             new_reps = row[IR_cols].drop(mxid)
    #             new_mean = new_reps.mean()
    #             new_cv = new_reps.std()/new_mean
    #             if new_cv <= 0.6:         
    #                 return pd.Series([new_mean,'CV okay if get rid of rep farthest from mean'], index=headers)
    #             else:
    #                 return pd.Series([np.nan,'CV too high'], index=headers)
    #         else:
    #             return pd.Series([combined_median,'keep; and meets target cutoff'], index=headers)
    #     #finally, default well-behaved proteins
    # else:
    #     return pd.Series([combined_median,'keep; passed through all filters'], index=headers)

def generate_combined(exp_group):
    #iterate over replicates
    #combine and generate mean SILAC ratio
    print('Working on:', exp_group, '_______')
    df = multicombine(exp_group)
    sequence_cols = [col for col in df.columns if 'sequence' in col]
    df['sequence'] = df[sequence_cols].apply(sequence_combine, axis = 1)
    IR_cols = [col for col in df.columns if 'IR_median' in col]
    qp_cols = [col for col in df.columns if 'No_qp' in col]
    combined_mean = df[IR_cols].mean(axis = 1)
    df['combined_no'] = df[IR_cols].apply(lambda values: np.count_nonzero(~np.isnan(values)), axis=1)
    df[exp_group + '_total_qp'] = df[qp_cols].apply(lambda values: np.sum(values), axis=1)
    df['combined_CV'] = df[IR_cols].std(axis = 1)/combined_mean
    df['combined_max'] = df[IR_cols].max(axis = 1)
    df['combined_min'] = df[IR_cols].min(axis = 1)
    # df = df.loc[df.combined_no >=2]
    combined_mean = df.apply(calculate_mean_IR, axis=1)
    combined_mean.columns = [exp_group + '_' + col for col in combined_mean.columns]
    combined_cols = [col for col in df.columns if 'combined_' in col]
    link_cols = [col for col in df.columns if 'link' in col]
    rep_mean_cols = [col for col in df.columns if 'IR_mean' in col]
    rep_std_cols = [col for col in df.columns if 'IR_std' in col]
    rep_CV_cols = [col for col in df.columns if 'CV_rep' in col]
    rep_noqp_cols = [col for col in df.columns if 'No_qp' in col]
    old_seq_cols = [col for col in df.columns if 'sequence_' in col]
    drop_cols = combined_cols + link_cols + rep_mean_cols + old_seq_cols + rep_std_cols + rep_noqp_cols + rep_CV_cols
    df.drop(drop_cols, axis = 1, inplace = True)
    df1 = df.join(combined_mean)
    df1.columns = df1.columns.str.replace('IR_median',exp_group)
    data_cols = sorted([col for col in df1.columns if col not in protein_cols])
    df1 = df1[protein_cols + data_cols]
    median_col = [col for col in df1.columns if 'combined_median' in col]
    print(df1[median_col])

    return df1

# def plot(df, exp_group):
#     df = df.loc[df['combined_mean'] != 'NQ']
#     y_axis = 'MS1 area log2(' + exp_group + ')'
#     mask = df.combined_mean < 0.05
#     df.loc[mask, 'combined_mean'] = 0.05
#     y = np.log2(df.combined_mean)
#     x = np.arange(0,len(df))
#     y_diff = (df.combined_mean > target_cutoff) | (df.combined_mean < (1/target_cutoff))
#     plt.scatter(x[~y_diff], y[~y_diff], lw = 0.5, c = 'black', marker = '.')
#     plt.scatter(x[y_diff], y[y_diff], lw = 0.5, c = 'red', marker = '.')
#     plt.xlabel('Peptides', fontsize = 16, family = 'ARIAL')
#     plt.ylabel(y_axis,fontsize = 16, family = 'ARIAL')
#     plt.xticks(fontsize = 16, family = 'ARIAL')
#     plt.yticks(fontsize = 16, family = 'ARIAL')
#     plt.ylim(-5,5)
#     plt.savefig('{}/combined/combined_{}.png'.format(home_dir, exp_group), bbox_inches='tight', transparent=True)
#     plt.close()

def main():
    groups = list(set(expt_table.reps.str.rsplit( '_', n = 1, expand = True)[0]))
    for exp_group in groups:
        #combine datasets
        combined_df = generate_combined(exp_group)
        # plot(combined, exp_group)
        combined_df.to_csv('{}/1_scripts/combined_isotop/{}_{}.csv'.format(dir,date, exp_group), index = False)

if __name__=="__main__":
    main()