# Compare MERFISH or smFISH decoded spots with RNA counts in bulk RNA seq data
# Goal: produce a correlation coefficient to determine the quality of your FISH data based on how well it correlates with the bulk RNA-seq
# the previous version of the file Doug sent has a bunch of extra steps in it to help do the correlation between bulk counts and FISH counts for the Zhuang analysis.
# We loop through all genes, create the scatter plot, and fit for the correlation coefficient.


import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Import bulk RNA-seq data 
bulkseq_path = Path('/data/bulkRNA/IMG/cleanIMGcounts.csv')

bulkseq_df = pd.read_csv(bulkseq_path)
bulkseq_df = bulkseq_df.rename(columns={'Unnamed: 0':'Gene_Name'})


# Only select the relevant columns
# The RNA sequencing data we have had multiple immune stimulants added to IMG cells, but we only want the columns titled 'IMGbas#', base meaning nothing added to these cells.
bulkseq_df = bulkseq_df.drop(columns=['IMGlpc1', 'IMGlpc2', 'IMGlps1', 'IMGlps2', 'IMGpic1', 'IMGpic2', 'IMGpma1', 'IMGpma2'])



# smfish_path = Path('/data/smFISH/20251028_bartelle_smFISH_mm_microglia_newbuffers')
smfish_decoded_features = Path('/data/smFISH/20251028_bartelle_smFISH_mm_microglia_newbuffers/qi2labdatastore/all_tiles_filtered_decoded_features/decoded_features.parquet')
# segemented_results = merfish_path / Path('processed_v2') / Path('decoded') / Path('baysor_formatted_genes.csv')

smfish_df = pd.read_parquet(smfish_decoded_features)
print(smfish_df.head)


# merfish_counts_df = merfish_df['gene_id'].value_counts().reset_index()
# merfish_counts_df.columns = ['gene_id','counts']

# seq_df = pd.read_csv(bulkseq_path,sep='\t')
# rnaseq_counts_df = seq_df[['gene_short_name','FPKM']].copy()
# rnaseq_counts_df.rename(columns={'gene_short_name': 'gene_id'},inplace=True)

# merged_df = pd.merge(merfish_counts_df, rnaseq_counts_df, on='gene_id', how='left', indicator=True)
# plt.scatter(np.log10(merged_df['FPKM']),np.log10(merged_df['counts']), alpha=0.5)
# plt.ylabel('Log10(Counts)')
# plt.xlabel('Log10(FPKM)')
# plt.title('Log-Log Correlation Plot')
# plt.grid(True)
# plt.show()