#!/bin/bash

# Fail fast and safer bash
set -euo pipefail

echo "Processing Cancer Chromatin and Enhancer-Gene links..."

# Expect to run from data/chromatin_states/cancer/

# Create temporary directory for intermediate files
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Remove Enh from Chromatin States (col 4, will import them from Enh-Genelinks)
awk '$4 !~ /EnhG1|EnhG2|EnhA1|EnhA2|EnhWk|EnhBiv/ {print}' final_cancer_calls.bed > "$TMPDIR/non_enh_cancer_cres.bed"

# Keep only cols 1-4 to combine into temp_file1.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' "$TMPDIR/non_enh_cancer_cres.bed" > "$TMPDIR/temp_file1.bed"

# Keep cols 1,2,3,6 from cancer_genelinks_hg38.bed to combine into temp_file2.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $6}' ../../enh_gene_links/cancer/cancer_genelinks_hg38.bed > "$TMPDIR/temp_file2.bed"

cat "$TMPDIR/temp_file1.bed" "$TMPDIR/temp_file2.bed" > "$TMPDIR/all_cancer_cres.bed"

# Now merge overlapping features with same state (col 4)
sort -k1,1V -k2,2n -k4,4 "$TMPDIR/all_cancer_cres.bed" > "$TMPDIR/sorted_features.bed"

for feature in $(cut -f4 "$TMPDIR/sorted_features.bed" | sort -u); do
    grep -w "$feature" "$TMPDIR/sorted_features.bed" | bedtools merge -c 4 -o distinct >> "$TMPDIR/temp_merged.bed"
done

sort -k1,1V -k2,2n "$TMPDIR/temp_merged.bed" > merged_cancer_cres.bed

gzip merged_cancer_cres.bed
echo "Done processing Cancer Chromatin and Enhancer-Genelinks"