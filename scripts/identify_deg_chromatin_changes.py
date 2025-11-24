# Identify DEG from final_chromatin_changes.bed.gz
# Using 'gene_name' from MANE.GRCh38.v1.4.ensembl_genomic.gff 

import pandas as pd
import subprocess
import os

with open('gene_chromatin_changes.bed', 'w') as outfile:
    subprocess.run(['bedtools', 'intersect', '-a', 'final_chromatin_changes.bed.gz', '-b',
        '../../genomes/mane_transcripts_gene.gff', '-wa', '-wb'], stdout=outfile, check=True)

gene_changes = pd.read_csv('gene_chromatin_changes.bed', sep='\t', header=None)

genes = gene_changes[[4,13]]

# List Chromatin Changes that should be Upregulated and Downregulated Genes/Promoters
# Remember all changes relative from healthy to cancer samples

# This is the "_strict" subset, only the most extreme changes
chromatin_changes = [
    'Tx-ReprPC', 'Tx-ReprPCWk', 'Tx-Het', 'Tx-Quies',
    'ReprPC-Tx', 'ReprPCWk-Tx', 'Het-Tx', 'Quies-Tx',
]

change_scores = {
    'Tx-Het': -4,
    'Tx-ReprPC': -3,
    'Tx-ReprPCWk': -2,
    'Tx-Quies': -1,
    'Het-Tx': 4,
    'ReprPC-Tx': 3,
    'ReprPCWk-Tx': 2,
    'Quies-Tx': 1,
}

# Filter for selected chromatin changes
changed_genes_df = genes[genes[4].isin(chromatin_changes)].copy()
changed_genes_df.columns = ['chromatin_change', 'gene']

# Quantify chromatin change
changed_genes_df['change_score'] = changed_genes_df['chromatin_change'].map(change_scores)
changed_genes_df = changed_genes_df[['gene', 'chromatin_change', 'change_score']]

# Remove genes that appear more than once 
changed_genes_df = changed_genes_df.groupby('gene').filter(lambda x: len(x) == 1)

print(f'There are {changed_genes_df["gene"].nunique()} unique genes with chromatin changes.')

# Save to CSV files
changed_genes_df.to_csv('changed_genes_chromatin.csv', index=False)

# gzip gene_chromatin_changes.bed
subprocess.run(['gzip', 'gene_chromatin_changes.bed'], check=True)
    
print('Saved to changed_genes_chromatin.csv.')

