#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is for QC filtering channels for a given experimental group (i.e. Mitosis vs Asynch) for cys-ABPP at the replicate level
"""

import numpy as np
import pandas as pd

#to-do: for num: take num/denom (unfiltered) then do QC, if the ratio is way high or way low, keep it. Then do the opposite for denom!!

def qc_filtering(row, group):
    """QC filtering for individual rep """
    #output will be a series / dataframe with median, CV, category
    headers = [group+'_median', group+'_CV',  group+'_category']
    reps = [col for col in row.index if '_rep' in col]
    old_cv = [col for col in row.index if '_CV' in col][0]
    old_median = [col for col in row.index if '_median' in col][0]
    #if CV is np.nan, means there's only 1 rep
    #If less than 0.5 keep, if greater than 0.5, 
    #if all values are really big or really small keep
    #if there are enough reps can just take median
    #otherwise return np.nan
    if np.isnan(row[old_cv]) == True:
        return pd.Series([row[old_median], row[old_cv],'Keep, test'], index=headers)
    elif row[old_cv] <= 0.5:
        return pd.Series([row[old_median], row[old_cv],'Keep, passed through filters'], index=headers)
    elif row[old_cv] > 0.5:
        min = row[reps].min()
        max = row[reps].max()
        no_reps = row[reps].count() 
        if min >= 2:
            return pd.Series([row[old_median], row[old_cv],'Keep, CV high but all values greater than 2'], index=headers)
        elif max <= 0.5:
            return pd.Series([row[old_median], row[old_cv],'Keep, CV high but all values less than 0.5'], index=headers)
        # elif no_reps == 2:
        #     return pd.Series([np.nan,np.nan, 'CV too high, value dropped'], index=headers)
        # elif no_reps >2:
        #     mean_value = row[reps].mean()
        #     #index of value that is farthest from the mean value
        #     mean_value = row[reps].mean()
        #     mxid = abs(row[reps] - mean_value).idxmax()
        #     new_reps = row[reps].drop(mxid)
        #     new_median = new_reps.median()
        #     min = new_reps.min()
        #     max = new_reps.max()
        #     if new_median > 0:
        #         new_cv = new_reps.std()/new_reps.mean()
        #         if new_cv <= 0.5:
        #             return pd.Series([new_median, new_cv,'Dropped value farthest from mean, CV okay now'], index=headers)
        #         elif new_cv > 0.5:
        #             if min >= 2:
        #                 return pd.Series([row[old_median], row[old_cv],'Keep, CV high but all values greater than 2 after deleting farthest from mean'], index=headers)
        #             elif max <= 0.5:
        #                 return pd.Series([row[old_median], row[old_cv],'Keep, CV high but all values less than 0.5 after deleting farthest from mean'], index=headers)
        #             else:
        #                 return pd.Series([np.nan,np.nan, 'CV too high, value dropped'], index=headers)          
        else:
            return pd.Series([np.nan, np.nan,'CV too high, value dropped'], index=headers)