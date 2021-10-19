#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 14:11:10 2021

@author: estherkemper
"""
import os
import numpy as np
import pandas as pd
from difflib import SequenceMatcher
import re
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

"""creating tryptic fasta database"""

dir = os.getcwd()
date = '20210830'

fasta = pd.read_csv('{}/inputs/RS_2016_fasta_20210127.csv'.format(dir))
mann = pd.read_csv('{}/other_datasets/mann_all.csv'.format(dir))
mann_high = pd.read_csv('{}/other_datasets/mann_high.csv'.format(dir))
lpp_df = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
print('reading inputs')
# slim_kemp = slim_kemp.assign(start = slim_kemp.pos - slim_kemp.sequence.str.find('*') + 1)
# slim_kemp = slim_kemp.assign(end = slim_kemp.start + slim_kemp.sequence.str.len() - 2)


def mann_cleanup(df):
    """cleaning up column names and getting rid of Y aa and sites without intensity in control or Noco, etc"""
    #renaming mann columns, only considering things with intensity in contorl or nocodazole, only considering S/T sites
    df_ctrl_cols = [col for col in df.columns if 'Intensity C' in col ] 
    df_mit_cols = [col for col in df.columns if 'Intensity N' in col ]
    df_dict = {'Uniprot ID':'accession', 'Protein names':'description', 'Gene names':'gene', 'Sequence window':'sequence', 'Position':'pos', 'Amino acid': 'aa'}
    df = df.rename(df_dict, axis = 1)
    df1 = df[['accession','description','gene','pos','aa','sequence']+df_ctrl_cols + df_mit_cols]
    sum_ints = np.sum(df1[df_ctrl_cols + df_mit_cols], axis =1)
    df2 = df1.loc[sum_ints > 0]
    df2 = df2.loc[~(df2.aa == 'Y')]
    #there are 7,664 proteins at this point
    df2.drop(df_ctrl_cols + df_mit_cols, axis = 1, inplace = True)
    return df2
def high_mann_cleanup(df):
    """cleaning up column names and getting rid of Y aa and sites without intensity in control or Noco, etc"""
    #renaming columns in the reported >50% stoichioemtry, 829 accessions at this point
    df_dict = {'Protein ID':'accession_high', 'ProteinName':'description', 'GeneName':'gene', 'Co-existing sites: AminoAcid_Position(Occupancy)':'res_stoichiometry'}
    df = df.rename(df_dict, axis = 1)
    df1 = df[['accession_high','res_stoichiometry']].dropna(how = 'all', axis = 0)
    #there are multiple reps for a given accession in the high stoich dataset
    # df2 = df1.loc[~df1.duplicated(subset = 'accession_high', keep = 'first')]
    df2 = df1.drop_duplicates(subset = 'accession_high', keep = 'first')
    print('cleaning mann data')
    return df2

def mann_merge(row, high_id_list):
    """this is to merge the high stoichiometry data with the total Mann data
    Sometimes they only list one accession in high stoichiometry, but there are multiple listed in the whole data"""
    acc = list(row['accession'].split(';'))
    for x in acc:
        if x in high_id_list:
            return x
        else:
            pass
        
def choose_uniprot(df, high_df):
    high_id_list = list(high_df.accession_high)
    df['accession_mann'] = df.apply(mann_merge, high_id_list = high_id_list, axis =1) 
    # #this number should equal 829
    # check = df.loc[df.accession_mann.notnull()]
    # len(set(check.accession_mann))
    df1 = df.merge(high_df,left_on = 'accession_mann', right_on = 'accession_high', how = 'outer').drop('accession_high', axis = 1)
    #prep the sequence column to get rid of _ and multiple sequences
    df1['sequence'] = df1.sequence.str.split(';', expand = True)[0]
    df1['original_seq'] = df1['sequence'].copy()
    df1['sequence'] = df1.sequence.replace('_', '', regex = True)
    return df1

def find_accession(row, fasta_id):
    """#we will look for uniprot ID that matches with the ID in our fasta database
    #Mann has multiple Ids separated by ;, and also has isoforms with dashes, which we don't have"""
    acc = list(row['accession'].split(';'))
    if not acc:
        print('ERROR!!', row['accession'])
    for x in acc:
        if '-' in x:
            x = x.split('-')[0]
            if x in fasta_id:
                return x
            else:
                pass
        elif x in fasta_id:
            return x
        else:
            pass

def match_accession(df):
    print('matching uniprot IDs with our FASTA')
    fasta_id = list(set(fasta.accession)  )
    #accession1 is the uniprot ID that is the same between the two datasets
    df['accession1'] = df.apply(find_accession, fasta_id = fasta_id, axis = 1)
    return df

def check_match(df, col):
    """Sometimes the a row will only have one uniprot that doesn't match, but that uniprot ID is shared
    with another row that has other uniprots that match"""
    #all the sites that were able to find a match
    pre_matched = df.loc[df[col].notnull()]
    #all the sites that were not able to find a match
    pre_unassigned = df.loc[df[col].isnull()]
    matched_ids =set(pre_matched.accession_mann)
    #this is the mann accession that is in matched and unassigned
    #happens if the accession didn't have the ID that is in our fasta
    #will append these rows to matched and take them off of unassigned
    confused_ids = list(matched_ids.intersection(list(set(pre_unassigned.accession_mann))))
    confused_ids =  [i for i in confused_ids if i]
    actually_matching = pre_unassigned.loc[pre_unassigned.accession_mann.isin(confused_ids)]
    matched = pre_matched.append(actually_matching)
    unassigned = pre_unassigned.loc[~pre_unassigned.accession_mann.isin(confused_ids)]#.drop('accession1', axis = 1)
    return matched, unassigned

# def find_seq(row, fasta_seqs):
#     seq = str(row['sequence'])
#     for x in fasta_seqs:
#         if seq in x:
#             return x
#         else:
#             pass

# def match_sequence(df):
#     print('matching by sequence with our FASTA')
#     #we'll try to match by sequence now
#     fasta_seqs = list(set(fasta.sequence))
#     df['seq'] = df.apply(find_seq, axis =1, fasta_seqs = fasta_seqs)
#     return df

# def find_gene(row, fasta_genes):
#     gene = str(row['gene'])
#     acc = list(gene.split(';'))
#     for x in acc:
#         if x in fasta_genes:
#             return x
#         else:
#             pass
        
# def match_gene(df):
#     print('matching by gene name with our FASTA')
#     fasta_genes = list(set(fasta.gene))
#     df.replace('2-Sep','SEPT2', inplace = True)
#     df['gene1'] = df.apply(find_gene, fasta_genes = fasta_genes, axis =1)
#     return df

def match_seq(row):
    seq = row['sequence']
    fasta_seq = row['sequence_fasta']
    # total_len = len(seq) - 1
    match = SequenceMatcher(None, seq, fasta_seq, autojunk = False).find_longest_match(0, len(seq), 0, len(fasta_seq))
    start = int(match.b)
    end = int(match.b + match.size - 1)
    #reassigns the position from Mann based on our fasta database
    # original_pos = row['pos']
    #if the sequence window from Mann is not the standard 31
    #because the assigned site is too early or late in seq
    #match.a is where it is in the sequence that matches, #match.b is where it matches in fasta
    #usually the S/T is in the middle of the sequence, but if at end or beginning it's not!!
    if int(match.a) > 16:
        new_pos = -1
    elif int(match.a + match.size - 1) < 16:
        new_pos = -1
    elif row['original_seq'][0] == '_':
        no_blank = row['original_seq'].count('_')
        if no_blank + int(match.a) > 16:
            new_pos = -1
        else:
            mann_seq_pos = 16 - row['original_seq'].count('_')
            new_pos = start - int(match.a) + mann_seq_pos
    else:    
        new_pos = start - int(match.a) + 16
        
    new_seq = seq[match.a: match.a + match.size]
    # seq_len = len(match) - 1
    if match.size >= 15:
        if 'S' in new_seq or 'T' in new_seq:
            return pd.Series([start, end, new_seq, new_pos])
        else:
            return pd.Series([start, end, new_seq, new_pos])
    else:
        return pd.Series([-1, -1, new_seq, -1])
        
def transform(df):
    df[['start','end','seq','new_pos']] = df.apply(match_seq, axis = 1)
    #only get one entry for re-mapped peptides
    df['len'] = df['end'] - df['start']
    return df

def total_combine(mann, mann_high):
    mann2 = mann_cleanup(mann)
    mann_high2 = high_mann_cleanup(mann_high)
    mann3 = choose_uniprot(mann2, mann_high2)
    mann4 = match_accession(mann3)
    matched, unassigned = check_match(mann4, 'accession1')
    # unassigned = match_sequence(unassigned)
    # matched1, unassigned1 = check_match(unassigned, 'seq')
    # unassigned1 = match_gene(unassigned1)
    # matched2, unassigned2 = check_match(unassigned1, 'gene1')
    fasta1 = fasta[['accession','gene','sequence']]
    print('finding the overlap of Mann reported sequence with our FASTA')
    # for df in [matched, matched1]:#, matched2]:
    #     df.drop(['accession','gene'], axis = 1, inplace = True)
    matched.drop(['accession_mann','gene'], axis = 1, inplace = True)
    matched.rename({'accession':'accession_mann'}, axis = 1, inplace = True)
    df1 = matched.merge(fasta1, left_on = 'accession1', right_on = 'accession', how = 'inner', suffixes = ['','_fasta']).drop('accession1', axis =1)
    # df2 = matched1.merge(fasta1, left_on = 'seq', right_on = 'sequence', how = 'inner', suffixes = ['', '_fasta']).drop('seq', axis = 1)
    # df3 = matched2.merge(fasta1, left_on = 'gene1', right_on = 'gene', how = 'inner', suffixes = ['', '_fasta']).drop('gene1', axis = 1)
    # for df in [df1,df2]:#,df3]:
    #     df = transform(df)
    mann_combined = transform(df1)
    # mann_combined = pd.concat([df1,df2])#],df3])
    return mann_combined

# def pure_tryptic(row):
#     fasta = row['sequence_fasta']
#     pos = row['new_pos'] - 1
#     first_seq = fasta[:pos]
#     end_seq = fasta[pos:]
#     start = re.search('K|R', first_seq[::-1])
#     end = re.search('K|R', end_seq)
#     if start:
#         start = start.start()
#     else:
#         start = 0
#     if end:
#         end = end.end()
#         tryptic = fasta[(pos-start):(pos+end)]
#         tryp_start = fasta.find(tryptic)
#         tryp_end = tryp_start + len(tryptic)
#         tryp_range = str(tryp_start) + '-' + str(tryp_end)
#         tryptic1 = row['accession'] + '_' + str(tryp_range) + ':' + tryptic
#     else:
#         tryptic =  fasta[(pos-start):]
#         tryp_start = fasta.find(tryptic)
#         tryp_end = tryp_start + len(tryptic)
#         tryp_range = str(tryp_start) + '-' + str(tryp_end)
#         tryptic1 = row['accession'] + '_' + str(tryp_range) + ':' + tryptic
#     return tryptic1

def mann_quality_check(mann, mann_high):
    mann_combined = total_combine(mann, mann_high)
    mann_unassigned = mann_combined.loc[(mann_combined['start'] == -1) | (mann_combined['new_pos'] < 1)]
    mann_total = mann_combined.loc[~ ((mann_combined['start'] == -1) | (mann_combined['new_pos'] < 1))]
    mann_total = mann_total.assign(mann_accession_pos = mann_total.accession + '_' + mann_total.new_pos.astype(str))
    mann_total = mann_total.assign(mann_gene_pos = mann_total.gene + '_' + mann_total.new_pos.astype(str))
    mann_total = mann_total.reset_index(drop = True)
    
    mann_total['temp'] = mann_total['accession'] + '_' + mann_total['new_pos'].astype(str)
    mann_total1 = mann_total.sort_values('len', ascending = False)
    mann_total1 = mann_total1.drop_duplicates(subset = 'temp', keep = 'first')
    mann_total1.drop('temp', axis =1 , inplace = True)
    
    # mann_total1 = mann_total.assign(tryptic = mann_total.apply(pure_tryptic, axis = 1))
    # print('assigning pure tryptics to Mann data')
    return mann_total1

def kemper_merge(lpp_df):
    lpp_df = lpp_df.sort_values('accession_res', ascending = False)
    lpp_df['temp'] = lpp_df.accession_res.str.rsplit('_', n = 1, expand = True)[0]
    lpp_df1 = lpp_df.drop_duplicates(subset = 'temp', keep = 'first')#loc[lpp_df.temp.duplicated(keep = 'first')]
    lpp_df1.drop('temp', axis = 1, inplace = True)
    fasta1 = fasta[['accession','gene','sequence']]
    slim_kemp = lpp_df1[['accession','accession_res','sequence','gene_res']]
    slim_kemp = slim_kemp.assign(new_pos = slim_kemp.gene_res.str.split('_', expand = True)[1].astype(int))
    ekk_total1 = slim_kemp.merge(fasta1, on = ['accession'], how = 'left', suffixes = ['','_fasta'])
    # ekk_total1 = ekk_total.assign(tryptic = ekk_total.apply(pure_tryptic, axis = 1))
    print('assigning pure tryptics to Kemper data')
    return ekk_total1

# def suppfig2a_peptides_venn(ekk_total1, mann_total1):
#     #source data generated here
#     only_y = set(ekk_total1.tryptic) - set(mann_total1.tryptic)
#     only_x = set(mann_total1.tryptic) - set(ekk_total1.tryptic)    
#     total = list(set(list(mann_total1.tryptic) + list(ekk_total1.tryptic)))
#     both = list( set(total) - set(only_x) - set(only_y))
    
#     ekk_total1.loc[ekk_total1.tryptic.isin(list(only_y)), 'Quantified in'] = 'Only in this study'
#     mann_total1.loc[mann_total1.tryptic.isin(list(only_x)), 'Quantified in'] = 'Only in Sharma et al.'
#     ekk_total1.loc[ekk_total1.tryptic.isin(both), 'Quantified in'] = 'Found in both'
#     mann_total1.loc[mann_total1.tryptic.isin(both), 'Quantified in'] = 'Found in both'
#     sizes = [ len(only_x), len(only_y),len(both)]
#     print('Making venn diagram')
    
#     plt.figure(figsize=(3, 3))
#     out = venn2(subsets = sizes, set_labels = ('Mann et al.', 'This study'), set_colors = ('r','b'), alpha = 0.8)
#     for text in out.set_labels:
#         text.set_fontsize(10)
#     for text in out.subset_labels:
#         text.set_fontsize(10)
#     plt.title('Tryptic phosphopeptides\nquantified from phosphoenrichment', fontsize = 10, y=0.95)
#     plt.savefig('{}/figures/suppfig2/{}_suppfig2a_peptides_venn.pdf'.format(dir,date), bbox_inches= 'tight')
    
def suppfig2a_sites_venn(ekk_total1, mann_total1):
    #source data generated here
    ekk_total1['accession_pos'] = ekk_total1['accession'] + '_' + ekk_total1['new_pos'].astype(str)
    mann_total1['accession_pos'] = mann_total1['accession'] + '_' + mann_total1['new_pos'].astype(str)
    only_y = set(ekk_total1.accession_pos) - set(mann_total1.accession_pos)
    only_x = set(mann_total1.accession_pos) - set(ekk_total1.accession_pos)    
    total = list(set(list(mann_total1.accession_pos) + list(ekk_total1.accession_pos)))
    both = list( set(total) - set(only_x) - set(only_y))
    
    ekk_total1.loc[ekk_total1.accession_pos.isin(list(only_y)), 'Quantified in'] = 'Only in this study'
    mann_total1.loc[mann_total1.accession_pos.isin(list(only_x)), 'Quantified in'] = 'Only in Sharma et al.'
    ekk_total1.loc[ekk_total1.accession_pos.isin(both), 'Quantified in'] = 'Found in both'
    mann_total1.loc[mann_total1.accession_pos.isin(both), 'Quantified in'] = 'Found in both'
    sizes = [ len(only_x), len(only_y),len(both)]
    print('Making venn diagram')
    
    plt.figure(figsize=(3, 3))
    out = venn2(subsets = sizes, set_labels = ('Mann et al.', 'This study'), set_colors = ('r','b'), alpha = 0.8)
    for text in out.set_labels:
        text.set_fontsize(10)
    for text in out.subset_labels:
        text.set_fontsize(10)
    plt.title('Phosphosites quantified\nfrom phosphoenrichment', fontsize = 10, y=0.95)
    plt.savefig('{}/figures/suppfig2/{}_suppfig2a_sites_venn.pdf'.format(dir,date), bbox_inches= 'tight')
    
# def suppfig2a_proteins_venn(ekk_total1, mann_total1):
#     #source data generated here
#     only_y = set(ekk_total1.accession) - set(mann_total1.accession)
#     only_x = set(mann_total1.accession) - set(ekk_total1.accession)    
#     total = list(set(list(mann_total1.accession) + list(ekk_total1.accession)))
#     both = list( set(total) - set(only_x) - set(only_y))
    
#     ekk_total1.loc[ekk_total1.accession.isin(list(only_y)), 'Quantified in'] = 'Only in this study'
#     mann_total1.loc[mann_total1.accession.isin(list(only_x)), 'Quantified in'] = 'Only in Sharma et al.'
#     ekk_total1.loc[ekk_total1.accession.isin(both), 'Quantified in'] = 'Found in both'
#     mann_total1.loc[mann_total1.accession.isin(both), 'Quantified in'] = 'Found in both'
#     sizes = [ len(only_x), len(only_y),len(both)]
#     print('Making venn diagram')
    
#     plt.figure(figsize=(3, 3))
#     out = venn2(subsets = sizes, set_labels = ('Mann et al.', 'This study'), set_colors = ('r','b'), alpha = 0.8)
#     for text in out.set_labels:
#         text.set_fontsize(10)
#     for text in out.subset_labels:
#         text.set_fontsize(10)
#     plt.title('Proteins with phosphosites\nquantified from phosphoenrichment', fontsize = 10, y=0.95)
#     plt.savefig('{}/figures/suppfig2/{}_suppfig2a_proteins_venn.pdf'.format(dir,date), bbox_inches= 'tight')

def col_rename(mann_total1, ekk_total1):
    """this is for source data"""
    ekk_total1['accession_pos'] = ekk_total1['accession'] + '_' + ekk_total1['new_pos'].astype(str)
    mann_total1['accession_pos'] = mann_total1['accession'] + '_' + mann_total1['new_pos'].astype(str)
    mann_total2 = mann_total1[['accession','accession_mann','gene','accession_pos','original_seq','seq','pos','new_pos']]
    mann_dict = {'original_seq':'Reported sequence window, Sharma et al', \
                 'pos':'Position, Sharma et al.', \
                 'new_pos':'position', \
                 'seq':'Portion of Sharma et al. sequence matching our FASTA', \
                 'aa':'amino_acid',\
                'accession_mann': 'accession, Sharma et al.'}
    mann_total2 = mann_total2.rename(mann_dict, axis = 1)
    ekk_total2 = ekk_total1.drop([ 'sequence_fasta','gene_res'],axis = 1)
    ekk_total2 = ekk_total2.rename({'sequence':'Reported sequence window, this paper',\
                                    'new_pos':'position'}, axis = 1)
    combined_total = mann_total2.merge(ekk_total2, on = ['accession_pos','accession','position','gene'], how = 'outer')
    combined_total['gene_res'] = combined_total.gene + '_' + combined_total.position.astype(str)
    combined_total = combined_total[['accession','accession, Sharma et al.','gene','accession_pos','gene_res','position','Position, Sharma et al.',\
                    'Reported sequence window, Sharma et al','Reported sequence window, this paper','Portion of Sharma et al. sequence matching our FASTA']]
    this_paper = combined_total['Reported sequence window, this paper'].notnull()
    mann_paper = combined_total['Reported sequence window, Sharma et al'].notnull()
    combined_total.loc[this_paper, 'Quantified in'] = 'This paper'
    combined_total.loc[mann_paper, 'Quantified in'] = 'Sharma et al.'
    combined_total.loc[(this_paper)&(mann_paper), 'Quantified in'] = 'Both'
    combined_total.to_csv('{}/supplemental_tables/source_data/Kemper_etal_2021_source_suppfig2a_venndiagram.csv'.format(dir), index = False)

    
def main():
    mann_total1 = mann_quality_check(mann, mann_high)
    ekk_total1 = kemper_merge(lpp_df)
    # suppfig2a_peptides_venn(ekk_total1, mann_total1)
    suppfig2a_sites_venn(ekk_total1, mann_total1)
    # suppfig2a_proteins_venn(ekk_total1, mann_total1)
    col_rename(mann_total1, ekk_total1)
    print('Writing to excel')
    

main()

