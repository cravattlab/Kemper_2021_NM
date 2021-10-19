#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is for cleaning up census out for 6-lex, 10-plex, or 16-plex
"""
import numpy as np 
import pandas as pd
import os

dir = os.getcwd()

#2016 fasta database from radu

fasta = pd.read_csv('{}/inputs/RS_2016_fasta_20210127.csv'.format(dir))

def clean_census_out_6plex(dataset_name, inputs, phosphopep):
    """Get census_out from IP2 and get it into working form, for 10-plex, for 6-plex see 201908091_Exp357_TMT_processing.py"""
    #census_out.txt from IP2
    #Had to re-save .txt as tab-delimited otherwise wouldn't read... will figure this out later
    channels = map(str,range(1,7,1))
    all_channels = ['c' + x for x in channels]
    if phosphopep == 'yes':
        df = pd.read_csv('{}/inputs/phospho/{}.txt'.format(dir, dataset_name), sep = '\t', header = None, low_memory = False)
    else:
        df = pd.read_csv('{}/inputs/{}.txt'.format(dir, dataset_name), sep = '\t', header = None, low_memory = False)
    #separate headers from data
    header = df.loc[df[0] == 'H']
    header = list(header.iloc[13][2:44])
    header[40] = 'description'
    header = header + ['42'] + ['unique']
    df = df.loc[df[0] != 'H']
    df['unique'] = df[1].copy()
    #replace U and nan with NaN to fill in accession and description for all rows
    df[1] = df[1].replace('U',np.NaN)
    df[1] = df[1].replace('nan',np.NaN)
    df[1] = df[1].fillna(method = 'ffill')
    df[41] = df[41].fillna(method = 'ffill')
    #useless protein row
    #get rid of extra columns
    df = df[~df[0].str.contains("P")]
    df = df.iloc[:,1:44]
    df.columns = header
    df = df.dropna(how = 'all', axis = 1)
    #column names
    good_columns = ['UNIQUE','SEQUENCE', 'description', 'unique',
    'm/z_126.127725_int', 
    'm/z_127.12476_int', 
    'm/z_128.134433_int',
    'm/z_129.131468_int',
    'm/z_130.141141_int',
    'm/z_131.138176_int']
    #choose by only column names we wnat and rename columns
    df = df[good_columns]
    rename = ['accession','sequence', 'gene_description','unique'] + all_channels
    df.columns = (rename)
    df = df.loc[:, df.columns.notnull()]
    #get rid of channels with nan in title
    channels = [x for x in all_channels if str(x) != 'nan']
    df[channels] = df[channels].astype('int')
    # #bandaid: If we only want to analyze certain channels, make sure it has c in the column header
    # c1_channels = [col for col in channels if 'c' in col]
    # df = df[['accession','sequence', 'gene_description'] + c1_channels]
    # #change channel values to integer values, get rid of 0 values
    # df[c1_channels] = df[c1_channels].astype('int')
    #replace all 0s with NaN (we will not include 0 in mean calculations, etc)
    #replace empty space with nan
    df = df.replace(r'^\s+$', np.nan, regex=True)
    print(str(len(df.index)) +' total peptide entries for '+ dataset_name)
    #split gene_description into 2 columns
    df[['gene','description']] = df['gene_description'].str.split(' ', n = 1, expand = True)
    df = df.drop('gene_description', axis=1)
    #remove keratin and reverse sequences
    # ker_rev = df.loc[df.description.str.contains("Keratin")|(df.accession.str.contains("Reverse"))]
    # ker_rev.to_excel('{}/replicates/removed/keratin_reverse_{}.xlsx'.format(dir,dataset_name))
    df = df.loc[~df.description.str.contains("Keratin")]
    df = df.loc[~df.accession.str.contains("Reverse")]
    return df

def clean_census_out_10plex(dataset_name, inputs, phosphopep):
    """Get census_out from IP2 and get it into working form, for 10-plex, for 6-plex see 201908091_Exp357_TMT_processing.py"""
    #census_out.txt from IP2
    #Had to re-save .txt as tab-delimited otherwise wouldn't read... will figure this out later
    channels = map(str,range(1,11,1))
    all_channels = ['c' + x for x in channels]
    if phosphopep == 'yes':
        df = pd.read_csv('{}/inputs/phospho/{}.txt'.format(dir, dataset_name), sep = '\t', header = None, low_memory = False)
    else:
        df = pd.read_csv('{}/inputs/{}.txt'.format(dir, dataset_name), sep = '\t', header = None, low_memory = False)
    #separate headers from data
    header = df.loc[df[0] == 'H']
    header = list(header.iloc[13][2:64])
    header[60] = 'description'
    header = header + ['62'] + ['unique']
    df['unique'] = df[1].copy()
    df = df.loc[df[0] != 'H']
    #replace U and nan with NaN to fill in accession and description for all rows
    df[1] = df[1].replace('U',np.NaN)
    df[1] = df[1].replace('nan',np.NaN)
    df[1] = df[1].fillna(method = 'ffill')
    df[61] = df[61].fillna(method = 'ffill')
    #get rid of keratin and reverse sequences and useless protein row
    df = df[~df[0].str.contains("P")]
    #get rid of extra columns
    df = df.iloc[:,1:64]
    df.columns = header
    df = df.dropna(how = 'all', axis = 1)
    #column names
    good_columns = ['UNIQUE','SEQUENCE', 'description', 'unique',
    'm/z_126.127726_int',
    'm/z_127.124761_int',
    'm/z_127.131081_int',
    'm/z_128.128116_int',
    'm/z_128.134436_int',
    'm/z_129.131471_int',
    'm/z_129.13779_int',
    'm/z_130.134825_int',
    'm/z_130.141145_int',
    'm/z_131.13818_int']
    #choose by only column names we wnat and rename columns
    df = df[good_columns]
    rename = ['accession','sequence', 'gene_description','unique'] + all_channels
    df.columns = (rename)
    df = df.loc[:, df.columns.notnull()]
    #get rid of channels with nan in title
    channels = [x for x in all_channels if str(x) != 'nan']
    df[channels] = df[channels].astype('int')
    # #bandaid: If we only want to analyze certain channels, make sure it has c in the column header
    # c1_channels = [col for col in channels if 'c' in col]
    # df = df[['accession','sequence', 'gene_description'] + c1_channels]
    # #change channel values to integer values, get rid of 0 values
    # df[c1_channels] = df[c1_channels].astype('int')
    #replace all 0s with NaN (we will not include 0 in mean calculations, etc)
    #df= df.replace(0, np.NaN)
    #replace empty space with nan
    df = df.replace(r'^\s+$', np.nan, regex=True)
    print(str(len(df.index)) +' total peptide entries for '+ dataset_name)
    #split gene_description into 2 columns
    df[['gene','description']] = df['gene_description'].str.split(' ', n = 1, expand = True)
    df = df.drop('gene_description', axis=1)
    #remove keratin and reverse sequences
    ker_rev = df.loc[df.description.str.contains("Keratin")|(df.accession.str.contains("Reverse"))]
    # ker_rev.to_excel('{}/replicates/removed/keratin_reverse_{}.xlsx'.format(dir,dataset_name))
    df = df.loc[~df.description.str.contains("Keratin")]
    df = df.loc[~df.accession.str.contains("Reverse")]
    return df

def read_census(dataset_name, no_missed_sites, inputs, phosphopep):
    """Getting sequence and matching with fasta to determine position, half tryptics and missed cleavages"""
    plex = inputs.plex.item()
    if plex == 6:
        print('6plex detected')
        df = clean_census_out_6plex(dataset_name, inputs, phosphopep)
    elif plex == 10:
        print('10plex detected')
        df = clean_census_out_10plex(dataset_name, inputs, phosphopep)
    else:
        print('Error cannot figure out what multiplex to use')

    #replace mod with *, get rid of periods (defined in inputs table)
    #replace mod and find position
    mod = inputs.modification
    mod = '\(' + str(mod.item()) + '\)'
    phos_mod = '\(79.966331\)'
    if phosphopep == 'yes':
        df = df.loc[(df.sequence.str.contains(mod)) & (df.sequence.str.contains(phos_mod))]
    else:
        pass
    df['sequence'] = df['sequence'].str.replace(mod,'*', regex = True)
    df['sequence'] = df['sequence'].str.replace(phos_mod,'#', regex = True)
    df['old_sequence'] = df['sequence']
    #only looking at sequence between the '.' of the tryptic peptide entry
    df['sequence'] = df.sequence.str.slice(2, -2)
    #find pos of * in sequence, get rid of *
    # df['no_mods'] = df['sequence'].str.count('\*')
    # one_mod_df = df.loc[df['no_mods'] ==1]
    # two_mod_df = df.loc[df['no_mods'] ==2]
    # one_mod_df = one_mod_df.assign(pos = one_mod_df['sequence'].str.find('*', 0))
    #duplicate the rows for those peptides that have multiple modifications
    #Find position for the first mod, then the second mod, then bring everything back together (have to add for last b/c the * takes up space)
    # two_mod_df_first = two_mod_df.assign(pos = two_mod_df['sequence'].str.find('*', 0))
    # two_mod_df_last = two_mod_df.assign(pos = two_mod_df['sequence'].str.rfind('*', 0) - 1)
    # df = pd.concat([one_mod_df, two_mod_df_first, two_mod_df_last], ignore_index = True)
    temp_seq = df.sequence.str.replace('#', '', regex = True)
    df = df.assign(pos = temp_seq.str.find('*'))
    df['sequence'] = df['sequence'].str.replace('*','', regex = True)
    df = df.assign(phos_pos = df['sequence'].str.find('#', 0))
    df['sequence'] = df['sequence'].str.replace('#','')
    #merge with fasta data
    df = df.merge(fasta, how = 'left', on = ['accession'], suffixes = ['', '_fasta'])
    #find residue number by finding location of tryptic peptide(_x) within whole protein sequence(_y)
    #make new columns with accession_res
    df['tryptic_position'] = df.apply(lambda r: str(r['sequence_fasta']).find(r['sequence']), axis = 1)
    df['residue'] = df.tryptic_position + df.pos
    df['gene_res'] = df['gene_fasta'] + '_' + df['residue'].astype(str)
    df['accession_res'] = df['accession'] + '_' + df['residue'].astype(str) #+ '_' + df['no_mods'].astype(str)
    if phosphopep == 'yes':
        print('Calculating only for phosphopeptides with modified cysteines')
        df.loc[df.phos_pos != -1, 'phos_pos'] = df.tryptic_position + df.phos_pos
        df['accession_res'] =  df['accession_res'] + '_' + df.phos_pos.astype(str)
    else:
        print('Calculating for TMT-ABPP or phosphoenrichment (not mixed)')
        pass
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
    df = df.loc[df.old_sequence.str.slice(3,-3).str.count('K|R') <= no_missed_sites]
    tryp_missed = pd.concat([df,old_df]).drop_duplicates(keep=False)
    tryp_missed = pd.concat([tryp_missed, df_half_tryp]).drop_duplicates(keep=False)
    print(str(len(df_half_tryp.index))+ ' half tryptics and ' +str(len(tryp_missed))+ ' missed_cleavage')
    #drop unnecesary columns leftover from all this (_x is from our data, _y from fasta)
    df = df.drop(['pos', 'gene', 'gene_fasta', 'residue', 'tryptic_position','description',  'sequence', 'sequence_fasta'], axis = 1) #
    df.rename(columns = {'description_fasta':'description','old_sequence':'sequence'}, inplace = True)

    # with pd.ExcelWriter('{}/replicates/removed/{}_improper_tryptic.xlsx'.format(dir,dataset_name)) as writer:
    #     df_half_tryp.to_excel(writer, sheet_name = 'half_tryptic', index = False, float_format='%.2f')
    #     tryp_missed.to_excel(writer, sheet_name = 'half_missed', index = False,float_format='%.2f')
    #     old_df.to_excel(writer, sheet_name = 'unfiltered', index = False,float_format='%.2f')
    print('Determined position and filtered bad tryptic peptides')
    df['sequence'] = df.sequence.str.slice(2, -2)
    return df

# inputs = pd.read_csv('{}/inputs/peptide_inputs.csv'.format(dir), index_col = None)
# inputs.set_index('dataset', inplace=True)
# df = read_census('Exp511', inputs)
# print(df)

