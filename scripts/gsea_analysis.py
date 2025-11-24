# GSEA on DESeq2 Results, Changes in Chromatin State, and Overlap of the Two

import pandas as pd
import gseapy

# Fix for pandas 2.0 compatibility with gseapy
if not hasattr(pd.DataFrame, 'append'):
    def append_fix(self, other, ignore_index=False, **kwargs):
        return pd.concat([self, other], ignore_index=ignore_index, **kwargs)
    pd.DataFrame.append = append_fix

# Expect to be run from data/

# Import GMT file for GSEA
gmt_file = 'gsea_genesets/combined_gene_sets.gmt'

deseq2_results = pd.read_csv('gdc/deseq2_all_results.csv')
deseq2_results.rename(columns={deseq2_results.columns[0]: 'gene'}, inplace=True)
chromatin_changes = pd.read_csv('chromatin_states/overlap/changed_genes_chromatin.csv')

# Find overlap of deseq2_results and chromatin_changes on 'gene' column
overlap_og = pd.merge(deseq2_results, chromatin_changes, on='gene', how='inner')
print(f'There are {len(overlap_og)} overlapping genes originally.')

# Remove rows were 'log2FoldChange' and change_score are opposite signs
overlap = overlap_og[overlap_og['log2FoldChange'] * overlap_og['change_score'] > 0]
print(f'There are {len(overlap)} overlapping genes after removing genes with conflicting data.')
print(f'Therefore there were {len(overlap_og) - len(overlap)} genes removed.')

# Combine 'log2FoldChange' and 'change_score' into a single column as equal weighted sum
overlap['combined_score'] = 0.5*overlap['log2FoldChange'] + 0.5*overlap['change_score']
print(overlap.head())

# Save overlap to csv
overlap.to_csv('../results/overlap_deseq2_chromatin.csv', index=False)
print('Saved combined dataset to results/overlap_deseq2_chromatin.csv')

# Run GSEA on deseq2_results ordered by 'log2FoldChange'
deseq2_ranked = deseq2_results[['gene', 'log2FoldChange']].copy()
deseq2_ranked = deseq2_ranked.dropna(subset=['log2FoldChange'])
deseq2_ranked = deseq2_ranked.sort_values('log2FoldChange', ascending=False)

# Create a pandas Series with gene names as index and log2FoldChange as values for gseapy.prerank()
deseq2_ranked_genes = deseq2_ranked.set_index('gene')['log2FoldChange']

# Run prerank GSEA using combined GMT file
deseq2_enr = gseapy.prerank(rnk=deseq2_ranked_genes,
                     gene_sets=gmt_file,
                     permutation_num=250, 
                     outdir='../results/gsea/deseq2/',
                     format='pdf',
                     seed=17)

print(f'DESeq2 GSEA prerank results saved to: results/gsea/deseq2/')

# Run GSEA on chromatin_changes ordered by 'change_score'
chromatin_ranked = chromatin_changes[['gene', 'change_score']].copy()
chromatin_ranked = chromatin_ranked.dropna(subset=['change_score'])
chromatin_ranked = chromatin_ranked.sort_values('change_score', ascending=False)

# Create a pandas Series with gene names as index and change_score as values for gseapy.prerank()
chromatin_ranked_genes = chromatin_ranked.set_index('gene')['change_score']

# Run prerank GSEA using combined GMT file
chromatin_enr = gseapy.prerank(rnk=chromatin_ranked_genes,
                     gene_sets=gmt_file,
                     permutation_num=250, 
                     outdir='../results/gsea/chromatin/',
                     format='pdf',
                     seed=17)

print(f'Chromatin GSEA prerank results saved to: results/gsea/chromatin/')

# Run GSEA on overlap ordered by 'combined_score'
overlap_ranked = overlap[['gene', 'combined_score']].copy()
overlap_ranked = overlap_ranked.dropna(subset=['combined_score'])
overlap_ranked = overlap_ranked.sort_values('combined_score', ascending=False)

# Create a pandas Series with gene names as index and combined_score as values for gseapy.prerank()
overlap_ranked_genes = overlap_ranked.set_index('gene')['combined_score']

# Run prerank GSEA using combined GMT file
overlap_enr = gseapy.prerank(rnk=overlap_ranked_genes,
                     gene_sets=gmt_file,
                     permutation_num=250, 
                     outdir='../results/gsea/overlap/',
                     format='pdf',
                     seed=17)

print(f'Overlap GSEA prerank results saved to: results/gsea/overlap/')