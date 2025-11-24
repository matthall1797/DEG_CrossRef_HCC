#!/bin/bash

# Fail fast and safer bash
set -euo pipefail

echo "Combining Cancerous Chromatin State beds..."

# Create temporary directory for gunzipped files
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Gunzip input files to temp directory (preserving original inputs for Snakemake)
for gz_file in *_18_CALLS_segments.bed.gz; do
    gunzip -c "$gz_file" > "$TMPDIR/${gz_file%.gz}"
done

# Process files in temp directory
for file in "$TMPDIR"/*.bed; do
    filename=$(basename "$file")
    sample=$(echo "$filename" | sed 's/_18_CALLS_segments.bed//')
    awk -v sample="$sample" 'BEGIN{OFS="\t"} {print $0, sample}' "$file" > "$TMPDIR/${filename%.bed}_with_sample.bed"
done

cat "$TMPDIR"/*_with_sample.bed > "$TMPDIR/cancer_calls.bed"

sort -k1,1V -k2,2n "$TMPDIR/cancer_calls.bed" > "$TMPDIR/cancer_calls_sorted.bed"

rm "$TMPDIR"/*_with_sample.bed

awk 'BEGIN{OFS="\t"} \
    { \
        key=$1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":"$9; \
        samples[key]=(samples[key] ? samples[key] : "") "," $10 \
    } \
    END{ \
        for(k in samples) { \
            split(k,f,":"); \
            print f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],substr(samples[k],2) \
        } \
    }' "$TMPDIR/cancer_calls_sorted.bed" > "$TMPDIR/merged_cancer_calls.bed"

sort -k1,1V -k2,2n "$TMPDIR/merged_cancer_calls.bed" > final_cancer_calls.bed

echo "Done processing Cancer Chromatin States"


