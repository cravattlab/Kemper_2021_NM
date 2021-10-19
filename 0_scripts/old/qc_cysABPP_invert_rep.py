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
    if row[old_median] >= 1:
        # if np.isnan(row[old_cv]) == True:
        #     return pd.Series([row[old_median], row[old_cv],'Keep, test'], index=headers)
        if row[old_cv] <= 0.5:
            return pd.Series([row[old_median], row[old_cv],'Keep, passed through filters'], index=headers)
        elif row[old_cv] > 0.5:
            min = row[reps].min()
            max = row[reps].max()
            no_reps = row[reps].count() 
            if min >= 2:
                return pd.Series([row[old_median], row[old_cv],'Keep, CV high but all values greater than 2'], index=headers)
            elif max <= 0.5:
                return pd.Series([row[old_median], row[old_cv],'Keep, CV high but all values less than 0.5'], index=headers)
            else:
                return pd.Series([np.nan, np.nan,'CV too high, value dropped'], index=headers)
    elif row[old_median] < 1:
        new_reps = list(map((lambda x: 1/x), list(row[reps])))
        new_std = np.std(new_reps)
        new_mean = np.mean(new_reps)
        new_cv = new_std/new_mean
        # if np.isnan(row[new_cv]) == True:
        #     return pd.Series([row[old_median], new_cv,'Keep, test'], index=headers)
        if new_cv <= 0.5:
            return pd.Series([row[old_median], new_cv,'Keep, passed through filters'], index=headers)
        elif new_cv > 0.5:
            new_min = np.min(new_reps)
            new_max = np.max(new_reps)
            no_reps = row[reps].count() 
            if new_min >= 2:
                return pd.Series([row[old_median], new_cv,'Keep, CV high but all values greater than 2'], index=headers)
            elif new_max <= 0.5:
                return pd.Series([row[old_median], new_cv,'Keep, CV high but all values less than 0.5'], index=headers)
            else:
                return pd.Series([np.nan, np.nan,'CV too high, value dropped'], index=headers)
    else:
        print('error')