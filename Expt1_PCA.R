library(randomcoloR)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(vegan)
library(tibble)
library(MCMCglmm)
library(ggtext)
library(Polychrome)

physeq <- readRDS("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/phylogeny-data-analysis/physeq.rds")

# PCA Analysis Section -------------------------------------------------------

otu_data <- otu_table(physeq)
otu_data_log <- log1p(otu_data)



# Perform PCA using prcomp()
pca_results <- prcomp(t(otu_data_log), scale. = TRUE)
eig_vals <- eigenvals(pca_results)  # Eigenvalues
varexpl <- (eig_vals / sum(eig_vals)) * 100  # Variance explained

# Extract PCA scores for the first three components (PC1, PC2, and PC3)
pca_scores <- as.data.frame(pca_results$x[, 1:3])

# If you want to merge with your sample metadata:
pca_scores$SampleID <- rownames(pca_scores)

pca_scores <- merge(pca_scores, samples_df_filtered, by.x = "SampleID", by.y = "row.names", all.x = TRUE)


# separate treatment into components
pca_scores <- pca_scores %>%
    separate(treatment, into = c("ancestral", "derived"), sep = "-", remove = FALSE)

# assign colors to microbes
custom_colors <- c("#AA4499", "#DDCC77", "#88CCEE", "#117733")


pca_scores$lineage <- ifelse(grepl("^ancestral", pca_scores$treatment), "ancestral", "derived")

# Convert to factor
pca_scores$lineage <- factor(pca_scores$lineage, levels = c("ancestral", "derived"))

# Now plot without needing the "Yes"/"No" treatment group
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = lineage)) +
  geom_point(size = 2.5) +
  scale_color_manual(name = "Lineage",
                     values = c("ancestral" = "#1f78b4", "derived" = "#e31a1c")) +
  labs(x = paste("PC1 (", round(varexpl[1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(varexpl[2], 2), "%)", sep = "")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme_cowplot()



pca_plot 

#save plot
ggsave("Expt1_pca_plot.jpg", pca_plot, width = 8, height = 5)
ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/outputs/Expt1_pca_plot.jpg", pca_plot, width = 8, height = 5)









