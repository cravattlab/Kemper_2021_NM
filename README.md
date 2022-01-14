# Kemper_2021_NM

inputs/ contains the the FASTA in .csv form as well as a sample dataset input for the cysteines TMT-ABPP. This is the output from IP2.
inputs/peptide_inputs.csv contains the details for peptide level TMT-ABPP (not protein level)
inputs/isotop_peptide_inputs.csv contains the details for peptide level isoTOP-ABPP
inputs/protein_inputs.csv contains the details for protein level TMT-ABPP (combines all peptide level for given protein)
modification is the diff mod of the enrichment probe. plex is the multiplex number (6-plex, 10-plex). 
num and denom are the numerator and denominator of the ratio calculations. c1-c10 are the channel assignments.

Scripts are organized by 0_*/ 1_*/ 2_*/ and should be run in that order. 

0_scripts/ 
Takes individual datasets and cleans them up and calculates ratios. 
peptide_* has scripts for peptide level TMT-ABPP
protein_* has scripts for protein level TMT-ABPP
isoTOP_* has scripts for peptide level isoTOP

1_scripts/
Combines replicate datasets for a given experimental group.

2_scripts/
Combines different experimental groups for given figures
fig1_combine_2.py combines TMT-ABPP Mitosis/Asynch and protein expression Mitosis/Asynch
fig2_scripts/ contains combinations for Figs 2-4
fig2_phosphoenrichment_2.py must be run first. This is for the phosphopeptide enrichment in Mitosis +/- LPP vs Asynch
fig2_combine_2.py combines adapted, original, and basal.
fig2_denature_combine_2.py combines denatured and original
fig2_phospho_combine_2.py combines original and basal

3_scripts/
Categorizes reactivity vs expression, finds proteins with changing cysteines

4_fig_scripts/
Scripts used to generate each figure

other_datasets/ 
Datasets from other papers, generally used in 4_fig_scripts/

supplemental_tables/ 
Outputs to .csv or excel for inclusion in supplemental table. Formatting was done in excel.


