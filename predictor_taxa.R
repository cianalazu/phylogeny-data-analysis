library(phyloseq)
library(ggplot2)

# 1. Transform counts to relative abundance
physeq_relabund <- transform_sample_counts(physeq, function(x) x / sum(x))

# 2. Extract OTU table and transpose (samples rows, taxa columns)
otu_df <- as.data.frame(t(otu_table(physeq_relabund)))

# 3. Remove taxa with any non-finite values
nonfinite_taxa <- sapply(otu_df, function(x) any(!is.finite(x)))
otu_df_clean <- otu_df[, !nonfinite_taxa]

# 4. Filter taxa present in at least 10 samples
valid_taxa <- vapply(otu_df_clean, function(x) sum(x > 0, na.rm = TRUE) >= 1, logical(1))

otu_df_filtered <- otu_df_clean[, valid_taxa]

cat("Number of taxa retained after cleaning and filtering:", sum(valid_taxa), "\n")

# 5. Extract metadata and combine with filtered OTU data
meta_df <- as.data.frame(sample_data(physeq))
combined_df <- cbind(meta_df, otu_df_filtered)

# Make sure RGR column exists
if(!"RGR" %in% colnames(meta_df)) stop("RGR column not found in sample_data!")

# 6. Initialize results data frame
results <- data.frame(Taxon = character(), R2 = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# 7. Loop over taxa and run linear model predicting RGR
for (taxon in colnames(otu_df_filtered)) {
  try({
    model <- lm(RGR ~ combined_df[[taxon]], data = combined_df)
    summary_model <- summary(model)
    results <- rbind(results, data.frame(
      Taxon = taxon,
      R2 = summary_model$r.squared,
      p_value = coef(summary_model)[2, 4]
    ))
  }, silent = TRUE)
}

# 8. Adjust p-values for multiple comparisons using FDR
results$adj_p_value <- p.adjust(results$p_value, method = "fdr")

# 9. View significant taxa (adjusted p-value < 0.05)
significant_taxa <- subset(results, adj_p_value < 0.05)
print(significant_taxa)

# 10. Optional: Plot the relationship for the top significant taxon
if (nrow(significant_taxa) > 0) {
  top_taxon <- significant_taxa$Taxon[which.min(significant_taxa$adj_p_value)]
  ggplot(combined_df, aes_string(x = top_taxon, y = "RGR")) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = paste("RGR vs", top_taxon), x = top_taxon, y = "RGR")
}
