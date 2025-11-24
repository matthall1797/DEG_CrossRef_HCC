# downloaded data is in Data/gdc/gdc_download

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("dplyr"))

# Read the sample sheet CSV file
samples <- read.csv("sample_sheet.csv") 

# Clean sample names 
samples$name <- trimws(samples$sample)
samples$name <- gsub("-", "_", samples$name)

# Path to data
data_path <- "gdc_download.tar.gz"

# Build count matrix directly by column binding (much more memory efficient than merging)
# First, get gene names from first sample
first_sample <- TRUE
gene_names <- NULL
count_columns <- list()

for (i in 1:nrow(samples)) {
  sample_name <- samples$name[i]
  
  # Construct the path to each file
  file_path <- paste0(samples$id[i], "/", samples$file[i])
  
  # Extract and read the file from tar.gz
  temp_dir <- tempdir()
  untar(data_path, files = file_path, exdir = temp_dir)
  
  # Full path
  full_file_path <- file.path(temp_dir, file_path)
  
  # Read the TSV file (skip lines starting with #)
  df <- read.delim(full_file_path, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
  
  # Clean up extracted file
  unlink(file.path(temp_dir, samples$id[i]), recursive = TRUE)
  
  # Remove first 4 rows (exp statistics)
  df <- df[-c(1:4), ]
  
  # Extract gene names from first sample
  if (first_sample) {
    gene_names <- df[, 2]  # column 2 is gene_name
    first_sample <- FALSE
  }
  
  # Extract just the counts column (column 4 - unstranded)
  counts <- df[, 4]
  
  # Store counts as a vector in list (name it by sample_name for later column naming)
  count_columns[[sample_name]] <- counts
  
  # Remove dataframe from memory
  rm(df, counts)
  gc()
  
  cat("Added:", sample_name, "\n")
}

# Build matrix by column binding (assuming all samples have same genes in same order, these are)
count_matrix <- do.call(cbind, count_columns)

# Set row names to gene names
rownames(count_matrix) <- gene_names

# Convert to numeric matrix for DESeq2
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "numeric"

# Clean up
rm(count_columns, gene_names)
invisible(gc())

cat("Gene Count matrix:", dim(count_matrix), "\n")
head(count_matrix)

# Create metadata for DESeq2
sample_names <- colnames(count_matrix)

# Samples with 01A = cancer, 11A = healthy
condition <- ifelse(grepl("01A", sample_names), "cancer", "healthy")

metadata <- data.frame(
  condition = condition,
  row.names = sample_names
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

# Filter out lowly-expressed genes 
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes Removed:", (nrow(count_matrix) - nrow(dds)), "\n")

# Run DESeq2
dds <- DESeq(dds)

# Get results: calc change in cancer wrt healthy
res <- results(dds, contrast = c("condition", "cancer", "healthy"))

summary(res)
head(res)

results <- as.data.frame(res)
results <- arrange(results, padj)

sig <- results[which(results$padj < 0.1), ]
upsig <- sig[which(sig$log2FoldChange > 0), ]
downsig <- sig[which(sig$log2FoldChange < 0), ]

write.csv(results, file = "deseq2_all_results.csv", row.names = TRUE)
write.csv(sig, file = "deseq2_sigs.csv", row.names = TRUE)
write.csv(upsig, file = "deseq2_upregulated.csv", row.names = TRUE)
write.csv(downsig, file = "deseq2_downregulated.csv", row.names = TRUE)

