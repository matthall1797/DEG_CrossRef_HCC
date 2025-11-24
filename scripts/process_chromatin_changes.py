# Use this to remove false positives by filtering out changes that overlap with no-change regions

import pandas as pd
import subprocess
from pybedtools import BedTool
import os

all_changes = BedTool('all_chromatin_changes.bed.gz')

# Get unique features from column 4
unique_features = set()
for feature in all_changes:
    unique_features.add(feature[3])

# Create output directory for intermediate files
os.makedirs('filtering_temp', exist_ok=True)

# Step 1: Split data by feature type into changes and no_changes
print("\nSplitting data by feature type...")
for feature in unique_features:
    print(f"Processing {feature}...")
    
    # Filter for this feature
    feature_data = all_changes.filter(lambda x: x[3] == feature)
    
    # Split into changes (B != C) and no_changes (B = C)
    changes = []
    no_changes = []
    
    for row in feature_data:
        # Extract B and C from col5 (format: B-C)
        change_parts = row[4].split('-')
        if len(change_parts) == 2:
            B, C = change_parts
            if B == C:
                no_changes.append(row)
            else:
                changes.append(row)
    
    # Save to files
    if changes:
        changes_bed = BedTool(changes)
        changes_bed.saveas(f'filtering_temp/{feature}_changes.bed')
    
    if no_changes:
        no_changes_bed = BedTool(no_changes)
        no_changes_bed.saveas(f'filtering_temp/{feature}_no_changes.bed')
    
# Remove false positives using bedtools intersect
print("\nRemoving false positives...")
all_filtered_changes = []
original_counts = {}
filtered_counts = {}

for feature in unique_features:
    changes_file = f'filtering_temp/{feature}_changes.bed'
    no_changes_file = f'filtering_temp/{feature}_no_changes.bed'
    
    if os.path.exists(changes_file) and os.path.exists(no_changes_file):
        print(f"Filtering {feature}...")
        
        # Use bedtools intersect to remove changes that overlap with no_changes
        changes_bed = BedTool(changes_file)
        no_changes_bed = BedTool(no_changes_file)
        
        # Count original changes by type
        for row in changes_bed:
            change = row[4]
            original_counts[change] = original_counts.get(change, 0) + 1
        
        # Remove overlaps (false positives)
        filtered_changes = changes_bed.intersect(no_changes_bed, v=True)
        
        # Count filtered changes by type
        for row in filtered_changes:
            change = row[4]
            filtered_counts[change] = filtered_counts.get(change, 0) + 1
            all_filtered_changes.append(row)
        
    elif os.path.exists(changes_file):
        # Only changes, no no_changes to filter against
        changes_bed = BedTool(changes_file)
        for row in changes_bed:
            change = row[4]
            original_counts[change] = original_counts.get(change, 0) + 1
            filtered_counts[change] = filtered_counts.get(change, 0) + 1
            all_filtered_changes.append(row)
        
# Save filtered results
print(f"\nSaving filtered results...")
if all_filtered_changes:
    filtered_bed = BedTool(all_filtered_changes)
    filtered_bed.saveas('all_chromatin_changes_filtered.bed')
    
    # Sort the final file
    subprocess.run(['sort', '-k1,1V', '-k2,2n', '-o', 'final_chromatin_changes.bed', 'all_chromatin_changes_filtered.bed'], check=True)
    subprocess.run(['gzip', 'final_chromatin_changes.bed'], check=True)
    subprocess.run(['rm', 'all_chromatin_changes_filtered.bed'], check=True)
    
    print(f"Filtered changes saved: {len(filtered_bed)} total")
    
    # Create summary with removal counts
    print("\nChange counts (original -> filtered -> removed):")
    total_original = 0
    total_filtered = 0
    
    # Write summary file
    with open('filtered_chromatin_changes_summary.tsv', 'w') as f:
        f.write("Change\tOriginal\tFiltered\tRemoved\tRemoval_Rate\n")
        
        # Get all unique changes
        all_changes = set(original_counts.keys()) | set(filtered_counts.keys())
        
        for change in sorted(all_changes):
            original = original_counts.get(change, 0)
            filtered = filtered_counts.get(change, 0)
            removed = original - filtered
            removal_rate = (removed / original * 100) if original > 0 else 0
            
            total_original += original
            total_filtered += filtered
            
            print(f"  {change}: {original} -> {filtered} -> {removed} ({removal_rate:.1f}%)")
            f.write(f"{change}\t{original}\t{filtered}\t{removed}\t{removal_rate:.1f}%\n")
        
        # Add totals row
        total_removed = total_original - total_filtered
        total_removal_rate = (total_removed / total_original * 100) if total_original > 0 else 0
        f.write(f"TOTAL\t{total_original}\t{total_filtered}\t{total_removed}\t{total_removal_rate:.1f}%\n")
        print(f"TOTAL: {total_original} -> {total_filtered} ({total_removal_rate:.1f}%)")
    
    print("Summary saved to filtered_chromatin_changes_summary.tsv")


# Clean up temporary files
print("Cleaning up temporary files...")
subprocess.run(['rm', '-rf', 'filtering_temp'], check=True)

print("Done")