#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is for cleaning up census out for 6-lex, 10-plex, or 16-plex
"""
import numpy as np 
import pandas as pd
import os
from peptide_tmt_census_cleanup import clean_census_out_6plex
from peptide_tmt_census_cleanup import clean_census_out_10plex
# from peptide_tmt_census_cleanup import clean_census_out_16plex

dir = os.getcwd()
fasta = pd.read_csv('{}/inputs/RS_2016_fasta_20210127.csv'.format(dir))

#2016 fasta database from radu

def read_census(dataset_name, inputs):
    """Getting sequence and matching with fasta to determine position, half tryptics and missed cleavages"""
    plex = inputs.plex.item()
    if plex == 6:
        print('6plex detected')
        df = clean_census_out_6plex(dataset_name, inputs)
    elif plex == 10:
        print('10plex detected')
        df = clean_census_out_10plex(dataset_name, inputs)
    elif plex == 16:
        print('16plex detected')
        df = clean_census_out_16plex(dataset_name, inputs)
    else:
        print('Error cannot figure out what multiplex to use')
    #replace mod with *, get rid of periods (defined in inputs table)
    mod = str(inputs.modification.item())
    mod = '\(' + mod + '\)'
    df['sequence'] = df['sequence'].str.replace(mod,'*', regex = True)
    df['old_sequence'] = df['sequence']
    #only looking at sequence between the '.' of the tryptic peptide entry
    df['sequence'] = df.sequence.str.slice(2, -2)
    #merge with fasta data
    df = df.merge(fasta, how = 'left', on = ['accession'])
    #find residue number by finding location of tryptic peptide(_x) within whole protein sequence(_y)
    #make new columns with accession_res
    df['tryptic_position'] = df.apply(lambda r: str(r['sequence_y']).find(r['sequence_x']), axis = 1)
    df['tryp_end'] = df.tryptic_position + df.apply(lambda x: len(x['sequence_x']), axis = 1)
    old_df = df.copy()
    #Create condition for position to be the 1st or 2nd amino acid
    #at n terminus require either KR before cleavage or cleavage after methionine AND
    #at c-terminus require K/R before cleavage or - at the end (- means it's at the end of the sequence)
    position_condition = (df.tryptic_position <= 2)
    ntryptic_condition = (df.old_sequence.str.slice(0,1).str.contains('-|K|R', regex = True) == True)
    ctryptic_condition = (df.old_sequence.str.slice(-3,-2 ).str.contains('K|R', regex = True))|(df.old_sequence.str.slice(-1, ).str.contains('-', regex = False))  == True
    df = df.loc[position_condition|ntryptic_condition]
    df = df.loc[ctryptic_condition]
    #DF with peptides that do not end in K or R or -
    # DF that don't start with K or R or -
    #drop the stuff we don't need
    df_half_tryp = pd.concat([df,old_df]).drop_duplicates(keep=False)
    #AND less than 2 missed cleavaged sites, but don't count if it's R.RXXX
    #Make a DF with all the peptides that we've dropped (half_tryp_missed)
    df = df.loc[df.old_sequence.str.slice(3,-3).str.count('K|R') < 2]
    tryp_missed = pd.concat([df,old_df]).drop_duplicates(keep=False)
    tryp_missed = pd.concat([tryp_missed, df_half_tryp]).drop_duplicates(keep=False)
    print(str(len(df_half_tryp.index))+ ' half tryptics and ' +str(len(tryp_missed))+ ' missed_cleavage')
    #drop unnecesary columns leftover from all this (_x is from our data, _y from fasta)
    df = df.drop(['gene_x', 'gene_y', 'description_x', 'sequence_x', 'sequence_y'], axis = 1)
    df.rename(columns = {'description_y':'description','old_sequence':'sequence'}, inplace = True)

    # with pd.ExcelWriter('{}/replicates/removed/{}_improper_tryptic.xlsx'.format(dir,dataset_name)) as writer:
    #     df_half_tryp.to_excel(writer, sheet_name = 'half_tryptic', index = False, float_format='%.2f')
    #     tryp_missed.to_excel(writer, sheet_name = 'half_missed', index = False,float_format='%.2f')
    #     old_df.to_excel(writer, sheet_name = 'unfiltered', index = False,float_format='%.2f')
    print('Determined position and filtered bad tryptic peptides')
    df['sequence'] = df.sequence.str.slice(2, -2)
    return df


