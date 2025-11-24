from glob import glob

rule all:
    input:
        "data/chromatin_states/healthy/final_healthy_calls.bed",
        "data/chromatin_states/cancer/final_cancer_calls.bed",
        "data/enh_gene_links/healthy/healthy_genelinks_hg38.bed",
        "data/enh_gene_links/cancer/cancer_genelinks_hg38.bed",
        "data/chromatin_states/healthy/merged_healthy_cres.bed.gz",
        "data/chromatin_states/cancer/merged_cancer_cres.bed.gz",
        "data/chromatin_states/overlap/all_chromatin_changes.bed.gz",
        "data/chromatin_states/overlap/final_chromatin_changes.bed.gz",
        "data/chromatin_states/overlap/changed_genes_chromatin.csv",
        "data/gdc/deseq2_all_results.csv",
        "data/gdc/deseq2_sigs.csv",
        "data/gdc/deseq2_upregulated.csv",
        "data/gdc/deseq2_downregulated.csv",
        "results/overlap_deseq2_chromatin.csv"

rule preprocess_healthy_chromatin: 
    input: 
        glob("data/chromatin_states/healthy/*_18_CALLS_segments.bed.gz")

    output: 
        "data/chromatin_states/healthy/final_healthy_calls.bed"
    
    conda:
        "envs/bedtools.yaml"

    shell: 
        "bash -c 'cd data/chromatin_states/healthy/ && ../../../scripts/preprocess_healthy_chromatin_states.sh'"

rule preprocess_cancer_chromatin:
    input: 
        glob("data/chromatin_states/cancer/*_18_CALLS_segments.bed.gz")

    output: 
        "data/chromatin_states/cancer/final_cancer_calls.bed"

    conda:
        "envs/bedtools.yaml"

    shell: 
        "bash -c 'cd data/chromatin_states/cancer/ && ../../../scripts/preprocess_cancer_chromatin_states.sh'"

rule preprocess_healthy_genelinks: 
    input: 
        glob("data/enh_gene_links/healthy/*_collated_pred.tsv.gz")

    output: 
        "data/enh_gene_links/healthy/healthy_genelinks_hg38.bed"

    conda:
        "envs/bedtools.yaml"

    shell: 
        "bash -c 'cd data/enh_gene_links/healthy/ && ../../../scripts/preprocess_healthy_enhancer_gene_links.sh'"

rule preprocess_cancer_genelinks: 
    input: 
        glob("data/enh_gene_links/cancer/*_collated_pred.tsv.gz")

    output: 
        "data/enh_gene_links/cancer/cancer_genelinks_hg38.bed"

    conda:
        "envs/bedtools.yaml"

    shell: 
        "bash -c 'cd data/enh_gene_links/cancer/ && ../../../scripts/preprocess_cancer_enhancer_gene_links.sh'"

rule process_healthy_chromatin_genelinks: 
    input: 
        chromatin="data/chromatin_states/healthy/final_healthy_calls.bed",
        genelinks="data/enh_gene_links/healthy/healthy_genelinks_hg38.bed"

    output: 
        "data/chromatin_states/healthy/merged_healthy_cres.bed.gz"

    conda:
        "envs/bedtools.yaml"

    shell: 
        "bash -c 'cd data/chromatin_states/healthy/ && ../../../scripts/process_healthy_chromatin_genelinks.sh'"

rule process_cancer_chromatin_genelinks: 
    input: 
        chromatin="data/chromatin_states/cancer/final_cancer_calls.bed",
        genelinks="data/enh_gene_links/cancer/cancer_genelinks_hg38.bed"

    output: 
        "data/chromatin_states/cancer/merged_cancer_cres.bed.gz"

    conda:
        "envs/bedtools.yaml"

    shell: 
        "bash -c 'cd data/chromatin_states/cancer/ && ../../../scripts/process_cancer_chromatin_genelinks.sh'"

rule id_chromatin_changes: 
    input: 
        healthy="data/chromatin_states/healthy/merged_healthy_cres.bed.gz",
        cancer="data/chromatin_states/cancer/merged_cancer_cres.bed.gz"

    output: 
        "data/chromatin_states/overlap/all_chromatin_changes.bed.gz"

    conda:
        "envs/bedtools.yaml"
    
    shell: 
        "bash -c 'cd data/chromatin_states/ && python ../../scripts/identify_chromatin_changes.py'"

rule process_chromatin_changes: 
    input: 
        "data/chromatin_states/overlap/all_chromatin_changes.bed.gz"
    
    output: 
        "data/chromatin_states/overlap/final_chromatin_changes.bed.gz"
    
    conda:
        "envs/bedtools.yaml"
    
    shell: 
        "bash -c 'cd data/chromatin_states/overlap/ && python ../../../scripts/process_chromatin_changes.py'"

rule id_chromatin_deg:
    input: 
        "data/chromatin_states/overlap/final_chromatin_changes.bed.gz",
        "data/genomes/mane_transcripts_gene.gff"
    
    output:
        "data/chromatin_states/overlap/changed_genes_chromatin.csv"
    
    conda:
        "envs/bedtools.yaml"
    
    shell: 
        "bash -c 'cd data/chromatin_states/overlap/ && python ../../../scripts/identify_deg_chromatin_changes.py'"

rule gdc_deseq:
    input: 
        "data/gdc/gdc_download.tar.gz",
        "data/gdc/sample_sheet.csv"

    output:
        "data/gdc/deseq2_all_results.csv",
        "data/gdc/deseq2_sigs.csv",
        "data/gdc/deseq2_upregulated.csv",
        "data/gdc/deseq2_downregulated.csv"

    conda:
        "envs/deseq2.yaml"

    shell:
        "bash -c 'cd data/gdc/ && Rscript ../../scripts/gdc_deseq.r'"

rule gsea_analysis:
    input:
        "data/gdc/deseq2_all_results.csv",
        "data/chromatin_states/overlap/changed_genes_chromatin.csv",
        "data/gsea_genesets/combined_gene_sets.gmt"
    
    output:
        "results/overlap_deseq2_chromatin.csv"

    conda:
        "envs/gsea.yaml"
    
    shell:
        "bash -c 'cd data/ && python ../scripts/gsea_analysis.py'"