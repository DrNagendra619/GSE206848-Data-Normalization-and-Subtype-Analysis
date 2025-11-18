# 1. Load Libraries for GSE206848 Analysis
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)

# 2. Load Dataset (GSE206848)
# Fetches the GSE206848 dataset from GEO
gse <- getGEO("GSE206848", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])
expression_data <- exprs(gse[[1]])

# 3. Check for Missing Values in GSE206848
cat("Missing values in GSE206848 expression data:", sum(is.na(expression_data)), "\n")

# 4. Boxplot Before Normalization (GSE206848)
png(file = file.path("GSE206848_boxplot_before_normalization.png"))
boxplot(expression_data, main = "GSE206848 - Before Normalization", las = 2, outline = FALSE)
dev.off()

# 5. Apply Normalization to GSE206848
expression_data <- normalizeBetweenArrays(expression_data, method = "quantile")

# 6. Boxplot After Normalization (GSE206848)
png(file = file.path("GSE206848_boxplot_after_normalization.png"))
boxplot(expression_data, main = "GSE206848 - After Normalization", las = 2, outline = FALSE)
dev.off()

# 7. Filter Low-Expressed Genes from GSE206848
gene_means <- rowMeans(expression_data)
threshold <- quantile(gene_means, probs = 0.25)
expression_data <- expression_data[gene_means > threshold, ]
cat("GSE206848 Dimensions after filtering:", dim(expression_data), "\n")

# 8. Extract Subtypes for GSE206848
# NS = Normal, OAS = Osteoarthritis, RAS = Rheumatoid Arthritis
subtype <- dplyr::case_when(
  grepl("NS", metadata$title) ~ "Normal (Control)",
  grepl("OAS", metadata$title) ~ "Osteoarthritis (OA)",
  grepl("RAS", metadata$title) ~ "Rheumatoid Arthritis (RA)",
  TRUE ~ "Unknown"
)
subtype <- as.factor(subtype)

cat("GSE206848 Subtype distribution:\n")
table(subtype)

# 9. Generate PCA Plot for GSE206848
pca_result <- prcomp(t(expression_data), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Subtype = subtype)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype)) +
  geom_point(size = 3) +
  ggtitle("GSE206848 PCA Plot: OA vs RA vs Normal") +
  theme_minimal()

ggsave(file.path("GSE206848_pca_plot.png"), plot = pca_plot)

# 10. Save Processed Data for GSE206848
save(expression_data, metadata, subtype, file = file.path("GSE206848_processed_data_step1.RData"))
