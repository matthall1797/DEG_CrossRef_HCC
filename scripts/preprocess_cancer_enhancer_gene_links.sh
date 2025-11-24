#!/bin/bash

# Fail fast and safer bash
set -euo pipefail

echo "Preprocessing Cancer Enhancer-Gene links..."

# Create temporary directory for processing
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Expect to run from data/enh_gene_links/cancer
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

cat "$TMPDIR"/*_with_sample.tsv > "$TMPDIR/cancer_genelinks.tsv"

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
}' "$TMPDIR/cancer_genelinks.tsv" > "$TMPDIR/new_cancer_genelinks.tsv"

sort -k1,1V -k2,2n "$TMPDIR/new_cancer_genelinks.tsv" > "$TMPDIR/sorted_cancer_genelinks.tsv"

rm "$TMPDIR/cancer_genelinks.tsv"
rm "$TMPDIR/new_cancer_genelinks.tsv"
rm "$TMPDIR"/*_with_sample.tsv
mv "$TMPDIR/sorted_cancer_genelinks.tsv" "$TMPDIR/cancer_genelinks.tsv"

cut -f1-7 "$TMPDIR/cancer_genelinks.tsv" > "$TMPDIR/trim_cancer_genelinks.bed"

# Use liftOver to convert from hg17 (input) to hg19 to hg38
liftOver -bedPlus=4 "$TMPDIR/trim_cancer_genelinks.bed" ../liftOver/hg17ToHg19.over.chain "$TMPDIR/cancer_hg19.bed" "$TMPDIR/cancer_unmapped19.bed"
liftOver -bedPlus=4 "$TMPDIR/cancer_hg19.bed" ../liftOver/hg19ToHg38.over.chain "$TMPDIR/cancer_hg38.bed" "$TMPDIR/cancer_unmapped38.bed"

# Copy final output to working directory
cp "$TMPDIR/cancer_hg38.bed" cancer_genelinks_hg38.bed 

echo "Done processing Cancer Enhancer-Gene links"