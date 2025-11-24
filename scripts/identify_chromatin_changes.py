# Compare feature annotations between merged_healthy_cres.bed and merged_cancer_cres.bed 
# Differences will be changes in chromatin state between healthy and cancer samples


import pandas as pd
import subprocess
from pybedtools import BedTool
import os

# Expect to run from data/chromatin_states/
healthy = BedTool('healthy/merged_healthy_cres.bed.gz')
cancer = BedTool('cancer/merged_cancer_cres.bed.gz')

# Define the chromatin state changes
# Format: (healthy_state, cancer_state, change_name)
chromatin_changes = [
    # Promoters
    ('TssA', 'TssA', 'TssA-TssA'),
    ('ReprPC', 'ReprPC', 'ReprPC-ReprPC'),
    ('ReprPCWk', 'ReprPCWk', 'ReprPCWk-ReprPCWk'),
    ('Het', 'Het', 'Het-Het'),
    ('TssBiv', 'TssBiv', 'TssBiv-TssBiv'),
    ('Quies', 'Quies', 'Quies-Quies'),
    ('TssA', 'ReprPC', 'TssA-ReprPC'),
    ('TssA', 'ReprPCWk', 'TssA-ReprPCWk'),
    ('TssA', 'Het', 'TssA-Het'),
    ('TssA', 'TssBiv', 'TssA-TssBiv'),
    ('TssA', 'Quies', 'TssA-Quies'),
    ('ReprPC', 'TssA', 'ReprPC-TssA'),
    ('ReprPCWk', 'TssA', 'ReprPCWk-TssA'),
    ('Het', 'TssA', 'Het-TssA'),
    ('TssBiv', 'TssA', 'TssBiv-TssA'),
    ('Quies', 'TssA', 'Quies-TssA'),

    # Genes
    ('Tx', 'Tx', 'Tx-Tx'),
    ('Tx', 'Het', 'Tx-Het'),
    ('Tx', 'ReprPC', 'Tx-ReprPC'),
    ('Tx', 'ReprPCWk', 'Tx-ReprPCWk'),
    ('Tx', 'TxWk', 'Tx-TxWk'),
    ('TxWk', 'TxWk', 'TxWk-TxWk'),
    ('TxWk', 'Het', 'TxWk-Het'),
    ('TxWk', 'ReprPC', 'TxWk-ReprPC'),
    ('TxWk', 'ReprPCWk', 'TxWk-ReprPCWk'),
    ('Het', 'Tx', 'Het-Tx'),
    ('ReprPC', 'Tx', 'ReprPC-Tx'),
    ('ReprPCWk', 'Tx', 'ReprPCWk-Tx'),
    ('TxWk', 'Tx', 'TxWk-Tx'),
    ('Het', 'TxWk', 'Het-TxWk'),
    ('ReprPC', 'TxWk', 'ReprPC-TxWk'),
    ('ReprPCWk', 'TxWk', 'ReprPCWk-TxWk'),
    
    # Enhancers
    ('EnhA1', 'EnhA1', 'EnhA1-EnhA1'),
    ('EnhA2', 'EnhA2', 'EnhA2-EnhA2'),
    ('EnhWk', 'EnhWk', 'EnhWk-EnhWk'),
    ('EnhBiv', 'EnhBiv', 'EnhBiv-EnhBiv'),
    ('EnhA1', 'EnhA2', 'EnhA1-EnhA2'),
    ('EnhA1', 'EnhWk', 'EnhA1-EnhWk'),
    ('EnhA1', 'EnhBiv', 'EnhA1-EnhBiv'),
    ('EnhA1', 'Quies', 'EnhA1-Quies'),
    ('EnhA1', 'ReprPCWk', 'EnhA1-ReprPCWk'),
    ('EnhA1', 'ReprPC', 'EnhA1-ReprPC'),
    ('EnhA1', 'Het', 'EnhA1-Het'),
    ('EnhA2', 'EnhA1', 'EnhA2-EnhA1'),
    ('EnhA2', 'EnhWk', 'EnhA2-EnhWk'),
    ('EnhA2', 'EnhBiv', 'EnhA2-EnhBiv'),
    ('EnhA2', 'Quies', 'EnhA2-Quies'),
    ('EnhA2', 'ReprPCWk', 'EnhA2-ReprPCWk'),
    ('EnhA2', 'ReprPC', 'EnhA2-ReprPC'),
    ('EnhA2', 'Het', 'EnhA2-Het'),
    ('EnhWk', 'EnhA1', 'EnhWk-EnhA1'),
    ('EnhBiv', 'EnhA1', 'EnhBiv-EnhA1'),
    ('Quies', 'EnhA1', 'Quies-EnhA1'),
    ('ReprPCWk', 'EnhA1', 'ReprPCWk-EnhA1'),
    ('ReprPC', 'EnhA1', 'ReprPC-EnhA1'),
    ('Het', 'EnhA1', 'Het-EnhA1'),
    ('EnhWk', 'EnhA2', 'EnhWk-EnhA2'),
    ('EnhBiv', 'EnhA2', 'EnhBiv-EnhA2'),
    ('Quies', 'EnhA2', 'Quies-EnhA2'),
    ('ReprPCWk', 'EnhA2', 'ReprPCWk-EnhA2'),
    ('ReprPC', 'EnhA2', 'ReprPC-EnhA2'),
    ('Het', 'EnhA2', 'Het-EnhA2')
    
]

# List to store all results
all_changes = []

# Loop through each chromatin state change
for healthy_state, cancer_state, change_name in chromatin_changes:
    print(f"Processing {change_name}")
    
    # Filter each bed for specific state
    healthy_filtered = healthy.filter(lambda x: x[3] == healthy_state).saveas()
    cancer_filtered = cancer.filter(lambda x: x[3] == cancer_state).saveas()
    
    # Find regions that are healthy_state in healthy and cancer_state in cancer
    intersection = healthy_filtered.intersect(cancer_filtered, u=True)
    
    # Convert to dataframe and add change column
    intersection_df = intersection.to_dataframe()
    intersection_df['change'] = change_name
    
    # Add to the list of all changes
    all_changes.append(intersection_df)

# Combine all changes into one dataframe
master_df = pd.concat(all_changes, ignore_index=True)

# Save as bed
master_df.to_csv('all_chromatin_changes_unsorted.bed', sep='\t', index=False, header=False)

# Use Linux sort command to sort the bed file
subprocess.run(['sort', '-k1,1V', '-k2,2n', '-o', 'all_chromatin_changes.bed', 'all_chromatin_changes_unsorted.bed'], check=True)
subprocess.run(['gzip', 'all_chromatin_changes.bed'], check=True)
subprocess.run(['rm', 'all_chromatin_changes_unsorted.bed'], check=True)
os.makedirs('overlap', exist_ok=True)
subprocess.run(['mv', 'all_chromatin_changes.bed.gz', 'overlap/all_chromatin_changes.bed.gz'], check=True)

print(f"Total changes found: {len(master_df)}")
print("Counts of chromatin changes:")
print(master_df['change'].value_counts())

# Save summary as tsv file
master_df['change'].value_counts().to_csv('overlap/chromatin_changes_summary.tsv', sep='\t', header=False)
print("Done finding chromatin changes!")