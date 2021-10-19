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
    elif row[old_cv] <= 0.6:
        return pd.Series([row[old_median], row[old_cv],row[no_cols],'Keep'], index=headers)
    elif row[old_cv] > 0.6:
        if row[combined_min] >= 1.8:
            return pd.Series([row[old_median], row[old_cv],row[no_cols],'Keep, CV high but all values changing'], index=headers)
        elif row[combined_max] <= (1/1.8):
            return pd.Series([row[old_median], row[old_cv],row[no_cols],'Keep, CV high but all values changing'], index=headers)
        elif row[no_cols] < 4:
            if row[combined_max] < 1.8:
                if row[combined_min] > (1/1.8):
                    return pd.Series([row[old_median], row[old_cv],row[no_cols],'Keep, high CV but medians unchanging'], index=headers)
                else:
                    return pd.Series([row[old_median], row[old_cv],row[no_cols],'CV too high, do not interpret'], index=headers)
            else:
                return pd.Series([row[old_median], row[old_cv],row[no_cols],'CV too high, do not interpret'], index=headers)
        elif row[no_cols] >=4:
            mean_value = row[reps].mean()
            mxid = abs(row[reps] - mean_value).idxmax()
            new_reps = row[reps].drop(mxid)
            new_mean = new_reps.mean()
            new_cv = new_reps.std()/new_mean
            if row[combined_min] > 1.5 or row[combined_max] < (1/1.5):
                return pd.Series([row[old_median], row[old_cv],row[no_cols],'Keep, high CV but all above 1.5'], index=headers)
            elif new_cv <= 0.6:         
                return pd.Series([row[old_median], row[old_cv],row[no_cols],'CV okay if get rid of rep farthest from mean'], index=headers)
            elif row[old_median] >=2:
                mnid = row[reps].idxmin()
                new_reps = row[reps].drop(mnid)
                new_min = new_reps.min()
                if new_min >= 1.8:
                    return pd.Series([row[old_median], row[old_cv],row[no_cols],'All changing if get rid of min'], index=headers)
                else:
                    return pd.Series([row[old_median], row[old_cv],row[no_cols],'CV too high, do not interpret!'], index=headers)
            elif row[old_median] <= 0.5:
                mxid = row[reps].idxmax()
                new_reps = row[reps].drop(mxid)
                new_max = new_reps.max()
                if new_max <= (1/1.8):
                    return pd.Series([row[old_median], row[old_cv],row[no_cols],'All changing if get rid of max'], index=headers)
                else:
                    return pd.Series([row[old_median], row[old_cv],row[no_cols],'CV too high, do not interpret!!'], index=headers)
            else:
                return pd.Series([row[old_median], row[old_cv],row[no_cols],'CV too high, do not interpret!!'], index=headers)
        else:
            return pd.Series([row[old_median], row[old_cv],row[no_cols],'CV too high, do not interpret!!!'], index=headers)