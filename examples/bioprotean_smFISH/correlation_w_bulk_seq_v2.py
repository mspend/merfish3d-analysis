# Compare MERFISH or smFISH decoded spots with RNA counts in bulk RNA seq data
# Goal: produce a correlation coefficient to determine the quality of your FISH data based on how well it correlates with the bulk RNA-seq
# the previous version of the file Doug sent has a bunch of extra steps in it to help do the correlation between bulk counts and FISH counts for the Zhuang analysis.
# We loop through all genes, create the scatter plot, and fit for the correlation coefficient.


import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Read in bulk RNA-seq data 
bulkseq_path = Path('/data/bulkRNA/IMG/kallisto/bulk_RNA_seq_counts_tpm.csv')
bulkseq_df = pd.read_csv(bulkseq_path)

# Select the relevant columns
bulkseq_counts = bulkseq_df.drop(columns=['Unnamed: 0'])
bulkseq_counts = bulkseq_counts.rename(columns={'gene_name':'gene_id'})
print(bulkseq_counts.head())

# Read in smFISH spots from the csv of spots detected by Big-FISH
smfish_path = Path('/data/smFISH/20251028_bartelle_smFISH_mm_microglia_newbuffers/qi2labdatastore/big_fish/one_tile_2D/spots.csv')

smfish_df = pd.read_csv(smfish_path)
print(smfish_df.head())

# Obtain the number of counts per gene
smfish_counts = smfish_df['gene_id'].value_counts().reset_index()
smfish_counts = smfish_counts.rename(columns={'count': 'fish_spots'})
# print(smfish_counts)

# Merge the dataframes and make a scatter plot 
merged_df = pd.merge(smfish_counts, bulkseq_counts, on='gene_id', how='left', indicator=True)
print(merged_df)

plt.scatter(np.log10(merged_df['tpm']),np.log10(merged_df['fish_spots']), alpha=0.5)
# plt.scatter((merged_df['tpm']),(merged_df['fish_spots']), alpha=0.5)
plt.ylabel('Log10(Number of FISH Spots Per Gene)')
plt.xlabel('Log10(TPM from bulk RNA-seq)')
plt.title('Log-Log Correlation Plot')
plt.grid(True)
plt.show()
