# Compare MERFISH or smFISH decoded spots with RNA counts in bulk RNA seq data
# Goal: produce a correlation coefficient to determine the quality of your FISH data based on how well it correlates with the bulk RNA-seq
# the previous version of the file Doug sent has a bunch of extra steps in it to help do the correlation between bulk counts and FISH counts for the Zhuang analysis.
# We loop through all genes, create the scatter plot, and fit for the correlation coefficient.


import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Read in bulk RNA-seq data 
bulkseq_path = Path('/data/bulkRNA/IMG/cleanIMGcounts.csv')
bulkseq_df = pd.read_csv(bulkseq_path)
bulkseq_counts = bulkseq_df.rename(columns={'Unnamed: 0':'gene_id', 'IMGbas1': 'counts', 'IMGbas2': 'counts_rep2'})

# Select the relevant columns
# The RNA sequencing data we have had multiple immune stimulants added to IMG cells, but we only want the columns titled 'IMGbas#', base meaning nothing added to these cells.
bulkseq_counts = bulkseq_counts.drop(columns=['width', 'IMGlpc1', 'IMGlpc2', 'IMGlps1', 'IMGlps2', 'IMGpic1', 'IMGpic2', 'IMGpma1', 'IMGpma2'])

# Read in smFISH spots from the parquet file of decoded spots
smfish_path = Path('/data/smFISH/20251028_bartelle_smFISH_mm_microglia_newbuffers/qi2labdatastore/all_tiles_filtered_decoded_features/decoded_features.parquet')

smfish_df = pd.read_parquet(smfish_path)

# Obtain the number of counts per gene
smfish_counts = smfish_df['gene_id'].value_counts().reset_index()
smfish_counts = smfish_counts.rename(columns={'count': 'fish_spots'})

# Merge the dataframes and make a scatter plot 
merged_df = pd.merge(smfish_counts, bulkseq_counts, on='gene_id', how='left', indicator=True)
plt.scatter(np.log10(merged_df['counts']),np.log10(merged_df['fish_spots']), alpha=0.5)
plt.ylabel('Log10(Number of FISH Spots Per Gene)')
plt.xlabel('Log10(Raw RNA Counts)')
plt.title('Log-Log Correlation Plot')
plt.grid(True)
plt.show()


