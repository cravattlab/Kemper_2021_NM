#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 12:20:08 2021

@author: estherkemper
"""
import os
import numpy as np
import pandas as pd
import seaborn as sns, matplotlib.pyplot as plt
sns.set(style="ticks")
from matplotlib import rcParams
import matplotlib.colors as mc
from matplotlib_venn import venn2
from  matplotlib.ticker import PercentFormatter
from difflib import SequenceMatcher

plt.rcParams['pdf.fonttype'] = 42
font = {'family' : 'ARIAL',
        'weight' : 'normal',
        'size'   : 10}         
plt.rc('font', **font) #set the font style created
plt.rcParams.update({'text.color' : "black"})

target_cutoff = 2
lpp_cutoff = 2
asynch_cutoff = 1.6
unchanging_cutoff = 1.6
protein_cutoff = 1.6
adapted_cutoff = 1.8

protein_cols = ['unique','accession','description', 'protein', 'gene_res']
phospho_col = ['pos_seq_mann', 'asynch_stoichiometry_mann','mitosis_stoichiometry_mann', 'stoichiometry','pos_seq_olsen', 'asynch_stoichiometry_olsen','mitosis_stoichiometry_olsen']

dir = os.getcwd()
#set a target cutoff defined in table
date = '20210830'

x_med ='cysteine_Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
x_prot_med ='protein_Mitosis_LPP(-,-)/(+,-)_combined_median'
y_prot_med = 'protein_Mitosis/Asynch_combined_median'
adapted_med ='cysteine_Mitosis_LPP(-,+)/(+,+)_combined_median'

# df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
# df.columns = df.columns.str.replace('Mitosis_LPP','cysteine_Mitosis_LPP')

# lpp_df = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')

fasta = pd.read_csv('{}/inputs/RS_2016_fasta_20210127.csv'.format(dir))
mann = pd.read_csv('{}/other_datasets/mann_all.csv'.format(dir))
mann_high = pd.read_csv('{}/other_datasets/mann_high.csv'.format(dir))
lpp_df = pd.read_excel('{}/supplemental_tables/{}_phosphoenrichment_combined.xlsx'.format(dir, date), sheet_name = 'filtered')
# olsen = pd.read_csv('{}/other_datasets/olsen_all.csv'.format(dir))


mann_high_dict = {'Protein ID':'accession1', 'ProteinName':'description', 'GeneName':'gene', 'Co-existing sites: AminoAcid_Position(Occupancy)':'res_stoichiometry'}
mann_high = mann_high.rename(mann_high_dict, axis = 1)
mann_high1 = mann_high[['accession1','gene','res_stoichiometry']].dropna(how = 'all', axis = 0)

mann_ctrl_cols = [col for col in mann.columns if 'Intensity C' in col ] 
mann_mit_cols = [col for col in mann.columns if 'Intensity N' in col ]
mann_dict = {'Uniprot ID':'accession', 'Protein names':'description', 'Gene names':'gene', 'Sequence window':'sequence', 'Position':'pos', 'Amino acid': 'aa'}
mann = mann.rename(mann_dict, axis = 1)
mann1 = mann[['accession','description','gene','pos','aa','sequence']+mann_ctrl_cols + mann_mit_cols]
sum_ints = np.sum(mann1[mann_ctrl_cols + mann_mit_cols], axis =1)
mann2 = mann1.loc[sum_ints > 0]
mann2 = mann2.loc[~(mann2.aa == 'Y')]
mann2.drop(mann_ctrl_cols + mann_mit_cols, axis = 1, inplace = True)
mann_high1['res_stoichiometry'] = mann_high1.res_stoichiometry.str.split(';')
mann_high2 = mann_high1.groupby('accession1')['res_stoichiometry'].sum().reset_index()

high_id_list = list(mann_high2.accession1)
def mann_merge(row):
    acc = list(row['accession'].split(';'))
    for x in acc:
        if x in high_id_list:
            return x
        else:
            pass

mann2['accession_mann'] = mann2.apply(mann_merge, axis =1) 
mann2 = mann2.merge(mann_high2,left_on = 'accession_mann', right_on = 'accession1', how = 'outer').drop('accession1', axis = 1)
 
# olsen_dict = {'Uniprot ID (Swiss-Prot or TrEMBL)':'accession', 'Phosphosite position in protein':'pos','ENSEMBL ID':'ensembl',\
#              'Gene Names':'gene', 'Protein Names':'description', 'Amino Acid':'aa', 'Raw Sequence':'sequence'}
# olsen = olsen.rename(olsen_dict, axis = 1)
# olsen1 = olsen[['accession','ensembl','description','gene','pos', 'aa','sequence']]
# olsen2 = olsen1.loc[olsen1.accession.notnull()]
# olsen2 = olsen2.loc[~(olsen2.aa == 'Y')]
# olsen2['pos'] = olsen2['pos'].astype(int)


fasta_id = list(fasta.accession)
    
def find_accession(row):
    acc = list(row['accession'].split(';'))
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

mann2['accession1'] = mann2.apply(find_accession, axis = 1)
mann_match = mann2.loc[mann2.accession1.notnull()]
mann_nomatch = mann2.loc[mann2.accession1.isnull()].drop('accession1', axis = 1)
# olsen2['accession1'] = olsen2.apply(find_accession, axis =1)
# olsen_match = olsen2.loc[olsen2.accession1.notnull()]
# olsen_nomatch = olsen2.loc[olsen2.accession1.isnull()].drop('accession1', axis = 1)


# fasta_des = str(list(fasta.description))
# def find_ensembl(row):
#     gene = str(row['ensembl'])
#     acc = list(set(gene.split(';')))
#     for x in acc:
#         if x in fasta_des:
#             return x
#         else:
#             pass

# olsen_nomatch['ensembl1'] = olsen_nomatch.apply(find_ensembl, axis = 1)
# foo = list(set(olsen_nomatch.ensembl1))

# def find_desc(row):
#     desc = str(row['description'])
#     for x in foo:
#         x = str(x)
#         if x in desc:
#             return x
#         else:
#             pass

# fasta['ensembl'] = fasta.apply(find_desc, axis = 1)


# olsen_match1 = olsen_nomatch.loc[olsen_nomatch.ensembl1.notnull()]
# olsen_nomatch1 = olsen_nomatch.loc[olsen_nomatch.ensembl1.isnull()].drop('ensembl1', axis = 1)

fasta_seqs = list(set(fasta.sequence))
def find_seq(row):
    seq = str(row['sequence'])
    for x in fasta_seqs:
        if seq in x:
            return x
        else:
            pass
        
mann_nomatch['seq'] = mann_nomatch.apply(find_seq, axis =1)
mann_match1 = mann_nomatch.loc[mann_nomatch.seq.notnull()]
mann_nomatch1 = mann_nomatch.loc[mann_nomatch.seq.isnull()].drop('seq', axis =1)


fasta_genes = list(set(fasta.gene))
def find_gene(row):
    gene = str(row['gene'])
    acc = list(gene.split(';'))
    for x in acc:
        if x in fasta_genes:
            return x
        else:
            pass

mann_nomatch1.replace('2-Sep','SEPT2', inplace = True)
mann_nomatch1['gene1'] = mann_nomatch1.apply(find_gene, axis =1)
mann_match2 = mann_nomatch1.loc[mann_nomatch1.gene1.notnull()]
mann_nomatch2 = mann_nomatch1.loc[mann_nomatch1.gene1.isnull()].drop('gene1', axis =1)

fasta1 = fasta[['accession','gene','sequence']]
# fasta3 = fasta[['ensembl','gene', 'sequence']]
mann_match.drop(['accession','gene'], axis = 1, inplace = True)
mann_match1.drop(['accession','gene'], axis = 1, inplace = True)
mann_match2.drop(['accession','gene'], axis = 1, inplace = True)

df1 = mann_match.merge(fasta1, left_on = 'accession1', right_on = 'accession', how = 'inner', suffixes = ['','_fasta']).drop('accession1', axis =1)
df2 = mann_match1.merge(fasta1, left_on = 'seq', right_on = 'sequence', how = 'inner', suffixes = ['', '_fasta']).drop('seq', axis = 1)
df3 = mann_match2.merge(fasta1, left_on = 'gene1', right_on = 'gene', how = 'inner', suffixes = ['', '_fasta']).drop('gene1', axis = 1)


# olsen_match.drop(['accession','gene'], axis = 1, inplace = True)
# df3 = olsen_match.merge(fasta1, left_on = 'accession1', right_on = 'accession', how = 'inner', suffixes = ['','_fasta']).drop(['ensembl','accession1'], axis =1)
# olsen_match1.drop(['ensembl','accession','gene'], axis = 1, inplace = True)
# df4 = olsen_match1.merge(fasta3, left_on = 'ensembl1', right_on = 'ensembl', how = 'inner', suffixes = ['','_fasta']).drop(['ensembl','ensembl1'], axis =1)
# olsen_match2.drop(['accession','gene'], axis = 1, inplace = True)
# df5 = olsen_match2.merge(fasta1, left_on = 'gene1', right_on = 'gene', how = 'inner', suffixes = ['', '_fasta']).drop(['ensembl','gene1'], axis = 1)

   
def match_seq(row):
    seq = row['sequence']
    fasta_seq = row['sequence_fasta']
    total_len = len(seq) - 1
    match = SequenceMatcher(None, seq, fasta_seq, autojunk = False).find_longest_match(0, len(seq), 0, len(fasta_seq))
    start = int(match.b)
    #reassigns the position from Mann based on our fasta database
    new_pos = int(int(total_len/2) - int(match.a) + int(match.b) + 1)
    end = int(match.b + match.size - 1)
    new_seq = seq[match.a: match.a + match.size]
    seq_len = len(match) - 1
    if match.size > 10:
        if 'S' in new_seq or 'T' in new_seq:
            return pd.Series([start, end, new_seq, new_pos])
        else:
            return pd.Series([start, end, new_seq, new_pos])
    else:
        return pd.Series([-1, -1, new_seq, -1])
        
def transform(df):
    df['sequence'] = df['sequence'].replace('_', '', regex = True)
    # df['sequence'] = df['sequence'].replace('\(ox\)', '', regex = True)
    # df['sequence'] = df['sequence'].replace('\(ph\)', '', regex = True)
    # df['sequence'] = df['sequence'].replace('\(gl\)', '', regex = True)
    # df['sequence'] = df['sequence'].replace('\(ac\)', '', regex = True)
    df['sequence'] = df.sequence.str.split(';', expand = True)[0]
    df[['start','end','seq','new_pos']] = df.apply(match_seq, axis = 1)
    df['len'] = df['end'] - df['start']
    return df

for df in [df1,df2,df3]:
    df = transform(df)
mann_combined = pd.concat([df1,df2,df3])

mann_unassigned = mann_combined.loc[mann_combined['start'] == -1]
mann_total = mann_combined.loc[~(mann_combined['start'] == -1)]
mann_total = mann_total.assign(mann_accession_pos = mann_total.accession + '_' + mann_total.new_pos.astype(str))
mann_total = mann_total.assign(mann_gene_pos = mann_total.gene + '_' + mann_total.new_pos.astype(str))



slim_mann = mann_total[['accession','start','new_pos','mann_accession_pos','mann_gene_pos','end','seq','res_stoichiometry']]
slim_kemp = lpp_df[['accession','accession_res','sequence','gene_res']]
slim_kemp = slim_kemp.assign(pos = slim_kemp.gene_res.str.split('_', expand = True)[1].astype(int))
slim_kemp = slim_kemp.assign(start = slim_kemp.pos - slim_kemp.sequence.str.find('*') + 1)
slim_kemp = slim_kemp.assign(end = slim_kemp.start + slim_kemp.sequence.str.len() - 2)



def test1(row, mann_df):
    start = int(row.start)
    end = int(row.end)
    matching = mann_df.loc[(mann_df.new_pos <= end) & (mann_df.new_pos >= start)]
    if matching.empty:
        return pd.Series([np.nan])
    else:
        id_prep = (matching.mann_accession_pos + ':' + matching.seq).to_list()
        # id_prep = matching.set_index('mann_accession_pos')['seq'].to_dict()
        return pd.Series([id_prep])
    
# def test1(row, k_df):
#     pos = row.new_pos
#     matching = k_df.loc[(k_df.start <= pos) & (k_df.end >= pos)]
#     if matching.empty:
#         return pd.Series([np.nan])
#     else:
#         id_prep = matching.set_index('gene_res')['sequence'].to_dict()
#         return pd.Series([id_prep])
    
prot_list = list(set(slim_kemp.accession))
def test(df, df1):
    cols = list(df.columns)
    new_df = pd.DataFrame(columns = cols)
    for acc in prot_list:
        k_df = df.loc[df.accession == acc]
        m_df = df1.loc[df1.accession == acc]
        # m_df1 = m_df.assign(mann_assigned = m_df.apply(test1, args = [k_df], axis = 1))
        k_df1 = k_df.assign(mann_assigned = k_df.apply(test1, mann_df = m_df, axis = 1))
        new_df = new_df.append(k_df1)
    return new_df
   
#all mann data
slim_mann1 = slim_mann[['accession','res_stoichiometry']]
slim_mann1 = slim_mann1.drop_duplicates(subset = 'accession')
#mann data that is high stoich
slim_high1 = slim_mann1.loc[slim_mann1.res_stoichiometry.notnull()]
slim_high1 = slim_high1.drop_duplicates(subset = 'accession')

#takes my dataset, then shows mann_unassined column shows if Mann quantified that site
matched_k_df = test(slim_kemp, slim_mann)   
#my data with all of high stoichiometry Mann data
matched_k_df1 = matched_k_df.merge(slim_high1, how = 'outer', on = 'accession')  
#my data with all mann data
all_matched = matched_k_df.merge(slim_mann1, how = 'outer', on ='accession')

#PROTEIN INFO FOR QUANTIFIED IN MINE AND HIGH STOICH IN MANN
both_high = matched_k_df1.loc[(matched_k_df1.res_stoichiometry.notnull()) & (matched_k_df1.accession_res.notnull())]
only_ekk = matched_k_df1.loc[(matched_k_df1.res_stoichiometry.isnull()) & (matched_k_df1.accession_res.notnull())]
only_mann = matched_k_df1.loc[(matched_k_df1.res_stoichiometry.notnull()) & (matched_k_df1.accession_res.isnull())]
#PROTEIN INFO FOR QUANTIFIED IN MINE and QUANTIFIED IN MANN
both_all = all_matched.loc[(all_matched.res_stoichiometry.notnull()) & (all_matched.accession_res.notnull())]
only_ekk_all = all_matched.loc[(all_matched.res_stoichiometry.isnull()) & (all_matched.accession_res.notnull())]
only_mann_all = all_matched.loc[(all_matched.res_stoichiometry.notnull()) & (all_matched.accession_res.isnull())]



#mann IDS that we found overlapping with our data
matched_ids = pd.Series(quant.mann_assigned.sum())
matched_ids1 = list(set(matched_ids.str.split(':', expand = True)[0]))
#the overlap between Mann data and ours
overlap = slim_mann.loc[slim_mann.mann_accession_pos.isin(matched_ids1)]
#Mann data that is NOT in our data
only_m = slim_mann.loc[~slim_mann.mann_accession_pos.isin(matched_ids1)]
overlap_high = overlap.loc[overlap.res_stoichiometry.notnull()]
only_m_high = only_m.loc[only_m.res_stoichiometry.notnull()]
only_k_high = only_k.loc[only_k.res_stoichiometry.notnull()]
#Mann sites that reside on peptides we quantify
len(set(overlap.mann_accession_pos))
#no Mann sites that do not reside on peptides we quantify
len(set(only_m.mann_accession_pos))
#no Mann proteins that we also quant
len(set(overlap.accession))
#Mann proteins that we also quant
len(set(only_m.accession))
#Mann sites that reside on peptides we quantify
len(set(overlap_high.mann_accession_pos))
#no Mann sites that do not reside on peptides we quantify
len(set(only_m_high.mann_accession_pos))
#no Mann proteins that we also quant
len(set(overlap_high.accession))
#Mann proteins that we also quant
len(set(only_m_high.accession))

def pie(sizes,labels,title,sub):
    total = sum(sizes)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.4,2.4)
    ax1.pie(sizes, autopct = (lambda p: '{:.0f}'.format(p * total / 100)), startangle=90)
    ax1.legend(labels = labels, loc='center left', bbox_to_anchor=(0.02, -0.3),frameon=False, handletextpad=0.3, fontsize = 10)
    plt.title(title, fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrich_{}_pie.pdf'.format(dir, date,sub), bbox_inches='tight',transparent = True)

both = list(set(overlap.mann_accession_pos))
only_mann = list(set(only_m.mann_accession_pos))
sizes = [ len(only_mann), len(both)]
labels = ['Only quantified in Mann et al.', 'Quantified in Mann et al. AND\nthese studies']
pie(sizes, labels,'All quantified\nphosphopeptides in mitosis', 'sites')

both2 = list(set(overlap_high.mann_accession_pos))
only_mann2 = list(set(only_m_high.mann_accession_pos))
sizes2 = [ len(only_mann2), len(both2)]
labels2 = ['Only quantified in Mann et al.', 'Quantified in Mann et al. AND\nthese studies']
pie(sizes2, labels2,'All phosphopeptides with high\nstoichiometry phosphorylation in mitosis ', 'high_sites')


both1 = list(set(overlap.accession))
only_mann1 = list(set(only_m.accession))
only_kemp = list
protein_sizes = [ len(only_mann1), len(both1)]
pie(protein_sizes, labels,'All proteins with quantified\nphosphopeptides in mitosis','protein')

both3 = list(set(overlap_high.accession))
only_mann3 = list(set(only_m_high.accession))
protein_sizes3 = [ len(only_mann3), len(both3)]
pie(protein_sizes3, labels2,'All proteins high stoichiometry\nphosphorylation in mitosis','high_protein')

def fig_venn(sizes, labels, title):
    plt.figure(figsize=(3.2, 3.2))
    plt.title(title, fontsize = 10)
    out = venn2(subsets = sizes, set_labels = (labels), set_colors = ('g','b'), alpha = 0.9)
    for text in out.set_labels:
        text.set_fontsize(10)
    for text in out.subset_labels:
        text.set_fontsize(10)
        plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrich_sites_venn.pdf'.format(dir,date), bbox_inches= 'tight')
only_x = list(set(only_k_high.accession))
only_y = list(set(only_m_high.accession))
both = list(set(overlap_high.accession))
sizes3 = [ len(only_y), len(only_x),len(both)]
labels3 = ['Mann et al.', 'Quantified in this study']
fig_venn(sizes3, labels3, 'Proteins with high stoichiometry\nphosphorylation in mitosis by Mann et al.')



#mann IDS w high stoichiometry that we found overlapping with our data

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir, date), sheet_name = 'filtered')


phospho_df = pd.read_csv('{}/other_datasets/final_phospho_mann_olsen.csv'.format(dir))
mann_only = phospho_df.loc[phospho_df.mitosis_stoichiometry_mann.notnull()]


x_med ='Mitosis_LPP(-,-)/(+,-)_combined_median'
y_med = 'cysteine_Mitosis/Asynch_combined_median'
adapted_med ='Mitosis_LPP(-,+)/(+,+)_combined_median'

df = pd.read_excel('{}/supplemental_tables/{}_fig2_combined.xlsx'.format(dir, date), sheet_name = 'filtered')




def fig2e_venn_prep():
    #source data generated here
   
    pdf = phospho_df[['accession', 'mitosis_stoichiometry_mann', 'mitosis_stoichiometry_olsen']]
    pdf.loc[pdf.mitosis_stoichiometry_olsen.notnull(), 'mitosis_stoichiometry_mann'] = pdf['mitosis_stoichiometry_mann'] + '; '
    pdf['phospho'] = pdf['mitosis_stoichiometry_mann'].fillna('') + pdf['mitosis_stoichiometry_olsen'].fillna('')
    pdf['phospho'] = pdf.phospho.str.split('; ')
    drop_cols = [col for col in pdf.columns if 'stoichiometry' in col]
    pdf.drop(drop_cols, axis = 1, inplace = True)
    pdf1 = pdf.explode('phospho')
    pdf1 = pdf1.replace('\)','', regex = True)
    pdf1[['res', 'stoichiometry']] = pdf1.phospho.str.split('(', expand = True)
    pdf1['aa'] = pdf1.res.str[0]
    pdf1['accession_res'] = pdf1.accession + '_' + pdf1.res.str[1:]
    pdf2 = pdf1.drop(['phospho', 'res'], axis = 1)
    pdf3 = pdf2.loc[~pdf2.duplicated(keep = 'first')].reset_index(drop = True)
    
    
    # lpdf1 = lpp_df[['accession','accession_res', x_med, y_med, 'cysteine_Mitosis_LPP(-,-)/(+,-)_combined_no',\
    #         'cysteine_Mitosis/Asynch_combined_no','cysteine_Mitosis_LPP(-,-)/(+,-)_combined_CV',\
    #             'cysteine_Mitosis/Asynch_combined_CV']]
    # lpdf1['accession_res1'] = lpdf1.accession_res.str.split('_', n = 1, expand = True)[1]
    # lpdf1['accession_res1'] = lpdf1.accession_res1.str.rsplit('_', n = 1, expand = True)[0]
    # lpdf1 = lpdf1.sort_values('accession_res', ascending = True).reset_index(drop = True)
    # lpdf2 = lpdf1.loc[lpdf1.accession_res1.duplicated(keep = 'first')]
    
    lpp_df['accession_res1'] = lpp_df.accession_res.str.rsplit('_', n = 1, expand = True)[0]
    
    df1 = pdf3.loc[pdf3.accession.isin(list(set(lpp_df.accession)))]
    df2 = pdf3.loc[pdf3.accession_res.isin(list(set(lpp_df.accession_res1)))]
    
    #df1 the number of proteins that we quantify that mann or olsen have identified as having high stoich
    len(set(df1.accession))
    len(set(phospho_df.accession))
    len(set(lpp_df.accession))

    only_x = list(set(set(lpp_df.accession_res1) - set(pdf3.accession_res)))
    only_y = list(set(pdf3.accession_res) - set(lpp_df.accession_res1))
    total = list(set(list(lpp_df.accession_res1) + list(pdf3.accession_res)))
    both = list( set(total) - set(only_x) - set(only_y))
    sizes = [ len(only_y), len(only_x),len(both)]
    print(sizes)
    
    plt.figure(figsize=(3.2, 3.2))
    out = venn2(subsets = sizes, set_labels = ('High stoichiometry sites\nin Mann et al. or Olsen et al.', 'Any quantified sites\nthis paper'), set_colors = ('g','b'), alpha = 0.9)
    for text in out.set_labels:
        text.set_fontsize(10)
    for text in out.subset_labels:
        text.set_fontsize(10)
        plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrich_sites_venn.pdf'.format(dir,date), bbox_inches= 'tight')
        
    sizes_sub = [ len(only_y),len(both)]
    total = sum(sizes_sub)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.4,2.4)
    ax1.pie(sizes_sub, autopct = (lambda p: '{:.0f}'.format(p * total / 100)), startangle=90)
    ax1.legend(labels = ['Only in Mann et al./Olsen et al.', 'Sites also quantified\nin this paper'], loc='center left', bbox_to_anchor=(0.02, -0.3),frameon=False, handletextpad=0.3, fontsize = 10)
    plt.title('High stoichiometry phosphorylation sites in Mann et al./Olsen et al', fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrich_sites_pie.pdf'.format(dir, date), bbox_inches='tight',transparent = True)
    
    only_x = list(set(set(lpp_df.accession) - set(phospho_df.accession)))
    only_y = list(set(phospho_df.accession) - set(lpp_df.accession))
    total = list(set(list(lpp_df.accession) + list(pdf3.accession)))
    both = list( set(total) - set(only_x) - set(only_y))
    sizes1 = [ len(only_y), len(only_x),len(both)]
    print(sizes1)
    
    plt.figure(figsize=(3.2, 3.2))
    out = venn2(subsets = sizes1, set_labels = ('Proteins w/ high stoichiometry sites\nin Mann et al. or Olsen et al.', 'Proteins w/ quantified\nphosphosites in this paper'), set_colors = ('g','b'), alpha = 0.9)
    for text in out.set_labels:
        text.set_fontsize(10)
    for text in out.subset_labels:
        text.set_fontsize(10)
        plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrich_protein_venn.pdf'.format(dir,date), bbox_inches= 'tight')
    
    sizes_sub1 = [ len(only_y),len(both)]
    total1 = sum(sizes_sub1)
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(2.4,2.4)
    ax1.pie(sizes_sub1, autopct = (lambda p: '{:.0f}'.format(p * total1 / 100)), startangle=90)
    ax1.legend(labels = ['Only found in Mann et al./Olsen et al.', 'Proteins with quantified sites\nin this paper'], loc='center left', bbox_to_anchor=(0.02, -0.3),frameon=False, handletextpad=0.3, fontsize = 10)
    plt.title('Proteins w/ high stoichiometryphosphorylation\nsite in Mann et al/Olsen et als', fontsize = 10)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.savefig('{}/figures/suppfig2/potential/{}_phosphoenrich_protein_pie.pdf'.format(dir, date), bbox_inches='tight',transparent = True)

fig2f_pie(prots_df)


# correct_mann = mann_total.loc[(mann_total['pos'] >= mann_total['start']) & (mann_total['pos'] <= mann_total['end'])]
# incorrect_mann =     mann_total.loc[~((mann_total['pos'] >= mann_total['start']) & (mann_total['pos'] <= mann_total['end']))]
        
    

