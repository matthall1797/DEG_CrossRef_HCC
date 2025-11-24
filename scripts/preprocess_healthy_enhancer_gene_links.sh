#!/bin/bash

# Fail fast and safer bash
set -euo pipefail

echo "Preprocessing Healthy Enhancer-Gene links..."

# Create temporary directory for processing
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Expect to run from data/enh_gene_links/healthy
# Gunzip input files to temp directory (preserving original inputs)
for gz_file in *_collated_pred.tsv.gz; do
    gunzip -c "$gz_file" > "$TMPDIR/${gz_file%.gz}"
done

# Process files in temp directory
for file in "$TMPDIR"/*.tsv; do 
    filename=$(basename "$file")
    sample=$(echo "$filename" | sed 's/_collated_pred.tsv//') 
    awk  -v sample="$sample" 'BEGIN{OFS="\t"} {print $0, sample}' "$file" > "$TMPDIR/${filename%.tsv}_with_sample.tsv"
done 

cat "$TMPDIR"/*_with_sample.tsv > "$TMPDIR/healthy_genelinks.tsv"

# Convert Enhancer numbers to abbreviations to match
awk 'BEGIN {
    OFS = FS = "\t"
    map["E7"] = "EnhG1"
    map["E8"] = "EnhG2"
    map["E9"] = "EnhA1"
    map["E10"] = "EnhA2"
    map["E11"] = "EnhWk"
    map["E15"] = "EnhBiv"
}
{
    for (i = 1; i <= NF; i++) {
        if ($i in map) $i = map[$i]
    }
    print
}' "$TMPDIR/healthy_genelinks.tsv" > "$TMPDIR/new_healthy_genelinks.tsv"

sort -k1,1V -k2,2n "$TMPDIR/new_healthy_genelinks.tsv" > "$TMPDIR/sorted_healthy_genelinks.tsv"

rm "$TMPDIR/healthy_genelinks.tsv"
rm "$TMPDIR/new_healthy_genelinks.tsv"
rm "$TMPDIR"/*_with_sample.tsv
mv "$TMPDIR/sorted_healthy_genelinks.tsv" "$TMPDIR/healthy_genelinks.tsv"

cut -f1-7 "$TMPDIR/healthy_genelinks.tsv" > "$TMPDIR/trim_healthy_genelinks.bed"

# Use liftOver to convert from hg17 (input) to hg19 to hg38
liftOver -bedPlus=4 "$TMPDIR/trim_healthy_genelinks.bed" ../liftOver/hg17ToHg19.over.chain "$TMPDIR/healthy_hg19.bed" "$TMPDIR/healthy_unmapped19.bed"
liftOver -bedPlus=4 "$TMPDIR/healthy_hg19.bed" ../liftOver/hg19ToHg38.over.chain "$TMPDIR/healthy_hg38.bed" "$TMPDIR/healthy_unmapped38.bed"

# Copy final output to working directory
cp "$TMPDIR/healthy_hg38.bed" healthy_genelinks_hg38.bed 

echo "Done processing Healthy Enhancer-Gene links"