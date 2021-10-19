#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QC filtering for combining groups for Phosphos TMT
"""
import pandas as pd
import numpy as np
from functools import reduce


def qc_filtering(row):
    """QC filtering for variability across individual experiments"""
    #name is the experimental group name
    headers = ['qc_ratio_combined_median', 'qc_ratio_combined_CV', 'qc_ratio_combined_no', 'qc_ratio_category']
    reps = [col for col in row.index if '_median_Exp' in col]
    old_cv = [col for col in row.index if '_combined_CV' in col][0]
    old_median = [col for col in row.index if '_combined_median' in col][0]
    combined_min = [col for col in row.index if '_min' in col][0]
    combined_max = [col for col in row.index if '_max' in col][0]
    # mxid = [col for col in row.index if '_mxid' in col]
    # mnid = [col for col in row.index if '_mnid' in col]
    no_cols = [col for col in row.index if '_combined_no' in col][0]
    if np.isnan(row[old_cv]) == True:
        return pd.Series([row[old_median], row[old_cv],row[no_cols],'0 or 1 rep'], index=headers)
    elif row[old_median] >= 1:
        if row[old_cv] <= 0.6:
            return pd.Series([row[old_median], row[old_cv],row[no_cols],'Keep'], index=headers)
        elif row[old_cv] > 0.6:
            if row[combined_min] >= 1.8:
                return pd.Series([row[old_median], row[old_cv],row[no_cols],'Keep, CV high but all values changing'], index=headers)
            elif row[no_cols] > 4:
                no_passing = sum(map(lambda x: x >= 1.5, row[reps]))
                total_no = row[no_cols]
                percent_passing = no_passing/total_no
                if percent_passing >= 0.8:
                    return pd.Series([row[old_median], row[old_cv],row[no_cols],'Keep, high CV but 80 percent of values 1.5'], index=headers)
                else:
                    return pd.Series([row[old_median], row[old_cv],row[no_cols],'CV too high, do not interpret!'], index=headers)
            else:
                return pd.Series([row[old_median], row[old_cv],row[no_cols],'CV too high, do not interpret!'], index=headers)
        else:
            print('Error0')
            return pd.Series([np.nan, np.nan, np.nan,'Could not characterize!'], index=headers)
    elif row[old_median] < 1:
        new_reps = list(map((lambda x: 1/x), list(row[reps])))
        new_std = np.nanstd(new_reps)
        new_mean = np.nanmean(new_reps)
        new_cv = new_std/new_mean
        if new_cv <= 0.6:
            return pd.Series([row[old_median], new_cv, row[no_cols], 'Keep'], index=headers)
        elif new_cv > 0.6:
            if row[combined_min] >= 1.8:
                return pd.Series([row[old_median], new_cv, row[no_cols], 'Keep, CV high but all values changing'], index=headers)
            elif row[no_cols] > 4:
                no_passing = sum(map(lambda x: x >= 1.5, new_reps))
                total_no = row[no_cols]
                percent_passing = no_passing/total_no
                if percent_passing >= 0.8:
                    return pd.Series([row[old_median], new_cv, row[no_cols], 'Keep, high CV but 80 percent of values 1.5'], index=headers)
                else:
                    return pd.Series([row[old_median], new_cv, row[no_cols], 'CV too high, do not interpret!'], index=headers)
            else:
                return pd.Series([row[old_median], new_cv, row[no_cols], 'CV too high, do not interpret!'], index=headers)
        else:
            # print('Error???')
            return pd.Series([np.nan, np.nan, np.nan,'Could not categorize'], index=headers)
    else:
        print('Error1')
        return pd.Series([np.nan, np.nan, np.nan,'Could not categorize'], index=headers)
               

