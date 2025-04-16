library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(limma)
library(pheatmap)
library(GenomeInfoDb) 
library(ggrepel)
library(readr)
library(RColorBrewer)

path_wd <- "~/Downloads/SFTP/CLIPOME/"
setwd(path_wd)

# Replace with the paths to your BED files
bed_files <- list.files("~/Downloads/SFTP/CLIPOME/Input", pattern = "\\.bed$", full.names = TRUE)

# Read and normalize each BED file to NCBI seqlevels
bed_list <- lapply(bed_files, function(f) {
  df <- read.delim(f, header = FALSE)
  
  if (ncol(df) == 6) {
    # Standard BED6 format (e.g., iCLIP)
    score_col <- df[[5]]
  } else if (ncol(df) == 10) {
    # MACS2 narrowPeak format (e.g., eCLIP)
    score_col <- df[[7]]  # signalValue
  } else if (ncol(df) == 8) {
    # omniCLIP format
    score_col <- df[[5]]  # signalValue
  } else {
    stop(paste("Unsupported number of columns (", ncol(df), ") in:", f))
  }
  
  #COMMENTED OUT: bring back if memory maxes out
  #score_threshold <- quantile(score_col, 0.8)  # select the top 20% fraction of peaks
  #df <- df[score_col >= score_threshold, ]
  #COMMENTED OUT
  
  score_col <- df[[ifelse(ncol(df) == 10, 7, 5)]]
  
  gr <- GRanges(
    seqnames = df[[1]],
    ranges = IRanges(start = df[[2]] + 1, end = df[[3]]),  # BED is 0-based
    strand = df[[6]],
    score = score_col,
    name = df[[4]]
  )
  
  # Optionally Convert to NCBI-style seqlevels (e.g., "1" instead of "chr1")
  GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"
  
  return(gr)
})

names(bed_list) <- tools::file_path_sans_ext(basename(bed_files))

#BELOW CODE BLOCK ONLY RETAINS ANNOTATED PEAKS AS PER GTF/GFF PROVIDED (OPTIONAL)
# Load annotation file (GTF or GFF3)
annotation_file <- "~/Downloads/SFTP/omniCLIP_results/Homo_sapiens.GRCh38.113.chr.gtf"  # or .gff3
annot <- rtracklayer::import(annotation_file)

# Keep only relevant features (e.g., exons, CDS, genes)
#annot <- annot[annot$type %in% c("exon", "gene", "three_prime_utr")]
annot = annot[annot$type %in% c("exon", "transcript","gene", "three_prime_utr","CDS","five_prime_utr")]

# Standardize chromosome naming (e.g., UCSC-style to NCBI)
GenomeInfoDb::seqlevelsStyle(annot) <- "NCBI"

# Filter bed_list to include only overlaps with annotation
bed_list <- lapply(bed_list, function(gr) {
  gr <- gr[GenomicRanges::countOverlaps(gr, annot, ignore.strand = FALSE) > 0]
  return(gr)
})
#ABOVE CODE BLOCK ONLY RETAINS ANNOTATED PEAKS AS PER GTF/GFF PROVIDED (OPTIONAL)


# Assuming bed_list is a list of GRanges objects
all_gr <- bed_list[[1]]
for (i in 2:length(bed_list)) {
  all_gr <- c(all_gr, bed_list[[i]])
}

#Reduce resolution> Merge peaks that are closer than n bp
all_regions <- GenomicRanges::reduce(all_gr, min.gapwidth = 5)


# Check result
class(all_gr)
# Should return "GRanges"


# Function to extract strand-specific score vector over all_regions
get_scores <- function(gr, ref) {
  hits <- findOverlaps(ref, gr, ignore.strand = FALSE)
  scores <- rep(0, length(ref))  # default score = 0 if no overlap
  scores[queryHits(hits)] <- mcols(gr)$score[subjectHits(hits)]
  return(scores)
}

# Build the score matrix
score_matrix <- sapply(bed_list, get_scores, ref = all_regions)

# Quantile normalize using limma
score_matrix <- limma::normalizeBetweenArrays(score_matrix, method = "quantile")

  #ANNOTATE eclip ENCORE data sets:
  # Load ENCORE annotations
  annotations <- read_tsv("~/Downloads/SFTP/CLIPOME/eCLIP_peaks/ENCORE_annotations.tsv")

  # Create combined label: protein + cell line
  annotations <- annotations %>%
    mutate(Protein_CellLine = paste0(`Target of assay`, "_", `Biosample term name`))

  # Extract ENCODE file IDs from the 'Files' column
  annotations <- annotations %>%
    mutate(File_IDs = str_extract_all(Files, "ENCFF[0-9A-Z]+"))

  # Flatten into a long lookup table: one row per file ID per annotation
  file_lookup <- annotations %>%
    select(Protein_CellLine, File_IDs) %>%
    unnest(File_IDs)

  # Clean sample names to match ENCODE accessions
  clean_names <- colnames(score_matrix)
  clean_names <- sub("\\.bed$", "", clean_names)
  #clean_names <- sub("_genome\\.clippy\\..*", "", clean_names)

  # Join clean names with lookup table
  annotated_df <- tibble(ENCODE_ID = clean_names) %>%
    left_join(file_lookup, by = c("ENCODE_ID" = "File_IDs")) %>%
    mutate(Protein_CellLine = ifelse(is.na(Protein_CellLine), ENCODE_ID, Protein_CellLine))

  # Example: Filter for K562 or HepG2 samples only
  filtered_annotated_df <- annotated_df %>%
    filter(str_detect(Protein_CellLine, "K562|HepG2"))
  
  # Example: Filter for specific protein, e.g., TIA1
   filtered_annotated_df <- annotated_df %>%
     filter(str_detect(Protein_CellLine, "MATR3"))
  
  # Replace clean_names with mapped annotations
  clean_names <- annotated_df$Protein_CellLine
  
  # Get matching columns (these are the ENCODE IDs used as colnames before annotation)
  keep_columns <- annotated_df$ENCODE_ID %in% filtered_annotated_df$ENCODE_ID
  score_matrix <- score_matrix[, keep_columns]
  
  # Update clean_names to match filtered subset
  clean_names <- filtered_annotated_df$Protein_CellLine
  

  
  #ANNOTATE other peak files e.g. CLippy etc
  clean_names <- colnames(score_matrix)
  clean_names <- sub("\\.bed$", "", clean_names)
  clean_names <- sub("_genome_.*", "", clean_names)
  clean_names <- sub("_genome\\.clippy\\..*", "", clean_names)
  
# Transpose matrix for PCA (samples as rows)
score_matrix_t <- t(score_matrix)

# Remove constant columns (zero variance across regions)
var_per_column <- apply(score_matrix_t, 2, var)
score_matrix_t <- score_matrix_t[, var_per_column != 0]

# Perform PCA
pca <- prcomp(score_matrix_t, scale. = TRUE)

# Build plotting data frame
pca_df <- data.frame(
  Sample = clean_names,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

# Get percentage variance
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

lname = length(clean_names)
# Define grouping (replace with your actual condition labels)
sample_group <- c(rep("ARTR-omni", 2),rep("eCLIP", 2), rep("ARTR-clippy", 4),rep("ARTR-nfHB", 2), rep("iCLIP", 1))  # Adjust length as needed
#sample_group <- c(rep("eCLIP",lname))

# Add to PCA dataframe
pca_df$Group <- factor(sample_group)
pca_df$Sample <- factor(clean_names)

# Plot
p = ggplot(pca_df, aes(x = PC1, y = PC2, color = Sample, shape = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, show.legend = FALSE, max.overlaps = 100) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("Strand-specific PCA of BED Score Matrix") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 6),     
    legend.title = element_text(size = 6)     
  )
p
ggsave("pca_PTBP1_CLIPs_vs_ARTR_Annotated.png", plot = p, width = 8, height = 6, dpi = 300)

# Use the PCA-transformed data (samples are rows)
pca_data <- pca$x[, 1:5]  # You can adjust the number of PCs used

rownames(pca_data) <- clean_names  # use ENCORE annotation or raw ID as you prefer

# Compute distance matrix
dist_pca <- dist(pca_data)

# Perform clustering
hc_pca <- hclust(dist_pca, method = "complete")

# Plot with annotated labels
plot(hc_pca, main = "Hierarchical Clustering Based on PCA", cex = 0.5)

# Optional: z-score transform per region (row-wise)
z_score_matrix <- t(scale(t(score_matrix)))  # center and scale each region

# Assign clean sample names again (if not already done)
colnames(z_score_matrix) <- clean_names


# Transpose so samples are rows
dist_matrix <- dist(t(z_score_matrix))  # or t(z_score_matrix) if using z-scores

# Cluster
hc_matrix <- hclust(dist_matrix, method = "complete")

# Plot
plot(hc_matrix, main = "Sample Clustering from Z-Score Matrix",cex = 0.5)

# 1. Remove columns (samples) with all zeros or constant values
zero_var_cols <- apply(z_score_matrix, 2, function(x) var(x, na.rm = TRUE) == 0)
z_score_matrix <- z_score_matrix[, !zero_var_cols]

# 2. Replace NA or NaN values with 0 (or another value if you prefer)
z_score_matrix[is.na(z_score_matrix) | is.nan(z_score_matrix)] <- 0

# 3. Recalculate correlation matrix
cor_matrix <- cor(z_score_matrix)

# 4. Ensure the correlation matrix is valid
stopifnot(!any(is.na(cor_matrix)))

# 5. Distance matrix (1 - correlation)
dist_cor <- as.dist(1 - cor_matrix)

# 6. Plot the heatmap
pheatmap(
  cor_matrix,
  clustering_distance_rows = dist_cor,
  clustering_distance_cols = dist_cor,
  clustering_method = "complete",
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
  main = "Sample-to-Sample Correlation Heatmap",
  fontsize_row = 5,
  fontsize_col = 5,
  border_color = NA
)


# Compute distance matrix
sample_dist <- dist(t(z_score_matrix))  # samples as columns → transpose

# Convert to matrix for pheatmap
sample_dist_matrix <- as.matrix(sample_dist)

# Optional: create a color palette
dist_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)


# Plot sample-to-sample heatmap
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = dist_colors,
  main = "Sample-to-Sample Distance Heatmap",
  fontsize_row = 5,  
  fontsize_col = 5, 
  border_color = NA,
  display_numbers = FALSE
)


# # Calculate variance for each region (row)
# row_vars <- apply(z_score_matrix, 1, var)
# 
# # Determine cutoff for top n%
# top_cutoff <- quantile(row_vars, 0.9)
# 
# # Subset to top 20% most variable regions
# top_matrix <- z_score_matrix[row_vars >= top_cutoff, ]
# 
# # Set color palette
# heatmap_colors <- colorRampPalette(rev(brewer.pal(n = length(clean_names), name = "RdBu")))(100)
# 
# # Plot heatmap on subset
# pheatmap(
#   top_matrix,
#   color = heatmap_colors,
#   show_colnames = TRUE,
#   show_rownames = FALSE,
#   fontsize_col = 10,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   scale = "none",
#   main = "Top 5% Most Variable Regions (Quantile-normalized Z-scores)"
# )



