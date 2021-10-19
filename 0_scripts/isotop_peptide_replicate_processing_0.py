#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 08:17:16 2021

@author: estherkemper
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed May 8th

@author: EKK
"""

import pandas as pd
import numpy as np
import time
import os

dir = os.getcwd()

#2016 fasta database from radu
fasta = pd.read_csv('{}/inputs/RS_2016_fasta_20210127.csv'.format(dir))
fasta.drop('description', axis =1, inplace = True)

cimage_output_dir = '/dta/output/'
cimage_output_file_url='output_rt_10_sn_2.5.to_excel.txt'
output_dir= dir + '/0_scripts/replicates_isotop'
inputs = '/inputs/isotop_peptide_inputs.csv'

r2_cutoff = 0.8
target_cutoff = 2
date = '20210830'

#final headers for output
stats_col_headers = ['mr','ave','stdev','noqp','ip2','quant_ir','all_ir', 'quant_seq']

def read_cimage(dataset_url):
    """
    Read cimage output_to_excel file. Rename columns ipi->uniprot
    and add entry id column
    """
    df = pd.read_csv(dataset_url, sep='\t', engine='c', skipinitialspace=True)
    #split accession
    df[['gene','description']] = df['description'].str.split(' ', n = 1, expand = True)
    #subset the "hyperlink" crap from cimage
    df.link = df.link.str[14:-2]
    df.link = dataset_url.split(cimage_output_file_url)[0]+df.link
    df.reset_index(inplace=True)
    #rename columns
    df = df.rename(columns={'ipi':'accession'})
    #get rid of filename code
    df.columns = [col.split('.')[0] for col in df.columns]
    #drop 'index' column; namespace
    df = df.drop(['index','mass', 'charge','segment','entry', 'NP', 'INT'], axis = 1)

    return df

def find_residue(dataset_url):
    """find residue information"""
    df = read_cimage(dataset_url)
    #Get residue information
    df['old_sequence'] = df['sequence'].copy()
    df['sequence'] = df.sequence.str.slice(2, -2)
    #find pos of * in sequence, get rid of *
    df['pos'] = df['sequence'].str.find('*', 0) 
    df['sequence'] = df['sequence'].str.replace('*','')
    #merge with fasta data
    df = df.merge(fasta, how = 'left', on = 'accession', suffixes = ['', '_fasta'])
    #Find tryptic peptide within fasta sequence
    # add position of modified residue to position of tryptic peptide in fasta
    def search(r):
        return str(r['sequence_fasta']).find(r['sequence']) + r['pos']
    df['residue'] = df.apply(search, axis = 1)
    #df['residue'] = df.apply(lambda x: str(x['sequence_fasta']).find(x['sequence']) + x['pos']], axis = 1)
	#drop unnecesary columns
    df['gene_res'] = df['gene_fasta'] + '_' + df['residue'].astype(str)
    df['tryptic_position'] = df.apply(lambda r: str(r['sequence_fasta']).find(r['sequence']), axis = 1)
    #drop reverse matches and keratin; these are crap
    df = df.loc[~df.accession.str.contains('Reverse_')]
    df = df.loc[~df.description.str.contains('Keratin')]
    #serialize for easy lookup later

    #Create condition for position to be the 1st or 2nd amino acid
    #at n terminus require either KR before cleavage or cleavage after methionine AND
    #at c-terminus require K/R before cleavage or - at the end (- means it's at the end of the sequence)
    position_condition = (df.tryptic_position <= 2)
    ntryptic_condition = (df.old_sequence.str.slice(0,1).str.contains('-|K|R', regex = True) == True)
    ctryptic_condition = (df.old_sequence.str.slice(-3,-2 ).str.contains('K|R', regex = True))|(df.old_sequence.str.slice(-1, ).str.contains('-', regex = False))  == True
    df = df.loc[position_condition|ntryptic_condition]
    df = df.loc[ctryptic_condition]

    #AND less than 2 missed cleavaged sites, but don't count if it's R.RXXX
    #Make a DF with all the peptides that we've dropped (half_tryp_missed)
    df = df.loc[df.old_sequence.str.slice(3,-3).str.count('K|R') < 2]
    df = df.drop(['sequence','sequence_fasta', 'pos', 'gene', 'gene_fasta', 'residue', 'tryptic_position'], axis = 1)
    df.rename(columns = {'description_y':'description','old_sequence':'sequence'}, inplace = True)
    return df
   
#def nonoutlier_median(series):
    #"""
    #Return the median of values excluding 20s if non_20 values exists
    #Otherwise returns standard median
    #"""
    #series = series.values

    #if 20.0 in series:
        #non_20 = series[series!=20]

        #if non_20.size>0:
            #check if there are many many 20s; i.e. >80% of series. Should be rare
            #if (series[series==20].size > 0.8*series.size):
                #many 20s
                #return np.nanmedian(series)
            #else:
                #most cases
                #return np.nanmedian( series[np.where(series != 20.)] )

   # return np.nanmedian(series)

def my_agg(x):
    names = {
        'IR_median': x['IR'].median(),
        'IR_mean': x['IR'].mean(),
        'IR_std':  x['IR'].std(),
        'No_qp': x['IR'].count()}  
    return pd.Series(names, index = ['IR_median', 'IR_mean', 'IR_std', 'No_qp']) 

def calculate_stats(dataset_url):
    df = find_residue(dataset_url)
    #filter R2. drop zero values because they couldn't be quant by cimage
    df = df.loc[(df.R2 >= r2_cutoff) & (df.IR != 0)]
    #aggregate data and reset index
    unique_R = pd.DataFrame(df.groupby('gene_res').apply(my_agg))
    unique_R = unique_R.reset_index()
    # unique_R['CV'] = unique_R['IR_std']/unique_R['IR_mean']
    # unique_R['CV'] = unique_R.CV.fillna(0)
    
    # #if standard devation is greater than half the value of the median, then get rid of 20s for that gene_res
    # high_std = []
    # for x in unique_R[unique_R.CV > 0.5].gene_res:
    #     high_std.append((df.loc[df.gene_res == x]))
    # high_std = pd.concat(high_std, ignore_index=True)
    # high_std = high_std.loc[high_std.IR != 20]
    # high_std = pd.DataFrame(high_std.groupby('gene_res').apply(my_agg))
    # high_std = high_std.reset_index()
    
    # #aggregate gene_res with std/median < 0.5
    # unique_R = unique_R.loc[unique_R.CV <= 0.5]
    # unique_R = pd.concat([unique_R, high_std], sort = True)
    unique_set = pd.DataFrame(df.groupby(['gene_res'], as_index = False).first())
    unique_set = unique_set.drop(['IR', 'R2'], axis = 1)
    unique_set = unique_set.set_index('gene_res')
    peptides = unique_set.merge(unique_R, on = 'gene_res')
    #get rid of 20s with only 1 quantified MS2
    peptides.drop(peptides[(peptides.IR_median == 20) & (peptides.No_qp <2)].index, inplace = True, axis = 0)
    peptides.sort_values(by =  'IR_median', ascending = 1, inplace = True)
    peptides.reset_index()
    
    return peptides

def main():
    """
    Generate combined files (peptide.csv) for each experiment.
    """
    expt_table = pd.read_csv('{}/{}'.format(dir, inputs))
    #change expt_table to dictionary and set url as the index
    expt_table = expt_table.set_index('url').to_dict()
    expt_table = expt_table['reps']

    #for  url in the url list, get the full url
    #and get the value for a given url, the dataset_name

    for url in [*expt_table]:
        dataset_url = url+cimage_output_dir+cimage_output_file_url
        #infer dataset_name P
        dataset_name = expt_table.get(url)

        # #make scatterplot
        # def plot(df):
        #     df.sort_values(by =  'IR_median', ascending = 1, inplace = True)
        #     y_axis = 'MS1 area log2(' + dataset_name + ')'
        #     mask = df.IR_median < 0.05
        #     df.loc[mask, 'IR_median'] = 0.05
        #     y = np.log2(df.IR_median)
        #     x = np.arange(0,len(df))
        #     y_diff = (df.IR_median > target_cutoff) | (df.IR_median < (1/target_cutoff))
        #     plt.scatter(x[~y_diff], y[~y_diff], lw = 0.5, c = 'black', marker = '.')
        #     plt.scatter(x[y_diff], y[y_diff], lw = 0.5, c = 'red', marker = '.')
        #     plt.xlabel('Peptides', fontsize = 12, family = 'ARIAL')
        #     plt.ylabel(y_axis,fontsize = 12, family = 'ARIAL')
        #     plt.xticks(fontsize = 12, family = 'ARIAL')
        #     plt.yticks(fontsize = 12, family = 'ARIAL')
        #     plt.ylim(-5,5)
        #     plt.savefig('{}/peptides/{}.png'.format(output_dir, dataset_name), bbox_inches='tight', transparent=True)
        #     plt.close()

        #combine per experiment
        print('Working on:', dataset_name, dataset_url, '\n')
        peptides = calculate_stats(dataset_url)
        peptides.drop('level_0', axis = 1, inplace = True)
        peptides.to_csv('{}/{}_{}.csv'.format(output_dir, date, dataset_name), index = False)

if __name__=="__main__":
    main()