# Compare MERFISH or smFISH decoded spots with RNA counts in bulk RNA seq data
# Goal: produce a correlation coefficient to determine the quality of your FISH data based on how well it correlates with the bulk RNA-seq
# the previous version of the file Doug sent has a bunch of extra steps in it to help do the correlation between bulk counts and FISH counts for the Zhuang analysis.
# We loop through all genes, create the scatter plot, and fit for the correlation coefficient.


import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Read in bulk RNA-seq data 
bulkseq_path = Path('/data/bulkRNA/IMG/kallisto/')
bulkseq_1 = bulkseq_path / 'rep_1' / 'bulk_RNA_seq_counts_summary.csv'
bulkseq_df = pd.read_csv(bulkseq_1, index_col=0)
print(bulkseq_df.head())

# Read in full bulk RNA-seq data and get a table
bulkseq_1_full = bulkseq_path / 'rep_1' / 'bulk_RNA_seq_counts.csv'
bulkseq_1_full_df = pd.read_csv(bulkseq_1_full, index_col=0)

gene_names = ['Angpt1', 'Crmp1', 'Dcx', 'Hdac11', 'Itgam', 'Notch2', 'Ptch1', 'Serpine1', 'Shh', 'Slc1a3', 'Tcf12', 'Tek', 'Tgfb1', 'Tgfbi', 'Aif1', 'Nnat']
full_df = bulkseq_1_full_df[bulkseq_1_full_df['gene_name'].isin(gene_names)]


# Create a Series containing only the gene name and the TPM
# For some genes, thare are multiple target_ids meaning multiple transcripts that map to the same gene. 
# These could be isoforms, splice variations, etc
# I summed the counts for all transcripts mapping to the same gene. 
summary = full_df.groupby("gene_name")["tpm"].sum().reset_index()



# Read in smFISH spots from the csv of spots detected by Big-FISH
smfish_path = Path('/data/smFISH/20251028_bartelle_smFISH_mm_microglia_newbuffers/qi2labdatastore/big_fish/results/one_tile_2D/summary_spots.csv')
smfish_df = pd.read_csv(smfish_path, index_col=0)
print(smfish_df.head())

# change order of columns
smfish_df = smfish_df[['gene_name', 'fish_spots']]

# # When spots come from merfish3d_analysis parquet files
# # Obtain the number of counts per gene
# smfish_counts = smfish_df['gene_id'].value_counts().reset_index()
# smfish_counts = smfish_counts.rename(columns={'count': 'fish_spots'})
# # print(smfish_counts)

# Merge the dataframes and make a scatter plot 
merged_df = pd.merge(smfish_df, bulkseq_df, on='gene_name', how='left', indicator=True)
print(merged_df)






plt.scatter(np.log10(merged_df['tpm']),np.log10(merged_df['fish_spots']), alpha=0.5)
# plt.scatter((merged_df['tpm']),(merged_df['fish_spots']), alpha=0.5)
plt.ylabel('Log10(Number of FISH Spots Per Gene)')
plt.xlabel('Log10(TPM from bulk RNA-seq)')
plt.title('Log-Log Correlation Plot')
plt.grid(True)
plt.show()
