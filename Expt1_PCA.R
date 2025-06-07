if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")


install.packages("Polychrome")

physeq <- readRDS("physeq.rds")


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

# Data Preparation Section --------------------------------------------------

#set working directory
setwd("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution")

# Step 1: Prepare OTU matrix
otu_mat <- read.table("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/processed_otu_matrix.tsv", header = TRUE, row.names = 1, check.names = FALSE)
otu_mat <- as.matrix(otu_mat)

otu_mat_145 <- read.table("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/filtered_otu_matrix.tsv", header = TRUE, row.names = 1, check.names = FALSE)
otu_mat_145 <- as.matrix(otu_mat)


# Step 2: Prepare taxonomy matrix
tax_mat <- read.table("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/processed_taxonomy_matrix.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
tax_mat <- as.matrix(tax_mat)

tax_mat_145 <- read.table("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/filtered_taxonomy_matrix.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
tax_mat_145 <- as.matrix(tax_mat)

# Step 3: Prepare sample data  

syn_com_treatments <- read.table("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/syn_com_treatments.csv", header = TRUE, stringsAsFactors = FALSE)
samples_df <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/samples_df.csv")
exp1_data <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/exp1_data.csv")# example: data <- read.csv("your_file.csv")

samples_df <- samples_df |> 
  mutate(lineage = case_when(
    grepl("^ancestral", treatment) ~ "ancestral",
    grepl("^derived", treatment) ~ "derived",
    grepl("^blank", treatment) ~ "blank",
    TRUE ~ NA_character_
  ))

# Data Filtering Section ----------------------------------------------------

samples_df_filtered <- samples_df |> 
  filter(lineage %in% c("ancestral", "derived"))


# Remove rows where not ancestral or derived
otu_mat_subset <- otu_mat[, samples_df_filtered$sample.id]
otu_mat_145_subset <- otu_mat_145[, samples_df_filtered$sample.id]

otu_mat_subset_filtered <- otu_mat_subset[rowSums(otu_mat_subset) > 0, ]
otu_mat_145_subset_filtered <- otu_mat_145_subset[rowSums(otu_mat_145_subset) > 0, ]

# Subset the taxonomy matrix
tax_mat_filtered <- tax_mat[rownames(tax_mat) %in% rownames(otu_mat_subset_filtered), ]
tax_mat_145_filtered <- tax_mat_145[rownames(tax_mat_145) %in% rownames(otu_mat_145_subset_filtered), ]


# Ensure unique taxa names
rownames(tax_mat_filtered) <- make.unique(rownames(tax_mat_filtered))
rownames(tax_mat_145_filtered) <- make.unique(rownames(tax_mat_145_filtered))

# Filter OTU matrix based on taxonomy
otu_mat_final <- otu_mat_subset_filtered[rownames(otu_mat_subset_filtered) %in% rownames(tax_mat_filtered), ]
otu_mat_145_final <- otu_mat_145_subset_filtered[rownames(otu_mat_145_subset_filtered) %in% rownames(tax_mat_145_filtered), ]

# Filter sample data to include only ancestral and derived samples
samples_df_filtered <- samples_df |> 
  filter(lineage %in% c("ancestral", "derived"))

# Ensure OTU matrix columns match filtered samples
matched_ids <- intersect(colnames(otu_mat_final), samples_df_filtered$sample.id)
matched_ids_145 <- intersect(colnames(otu_mat_145_final), samples_df_filtered$sample.id)


otu_mat_final <- otu_mat_final[, matched_ids]
otu_mat_145_final <- otu_mat_145_final[, matched_ids_145]

# Reorder samples_df to match OTU matrix column order
samples_df_filtered <- samples_df |> 
  filter(sample.id %in% matched_ids)  |> 
  arrange(match(sample.id, matched_ids))

samples_df_filtered <- samples_df_filtered |> 
  filter(sample.id %in% matched_ids_145)  |> 
  arrange(match(sample.id, matched_ids_145)) 



# Relative Abundance Calculation ----------------------------------------

# Relative abundance for filtering
otu_mat_relative <- sweep(otu_mat_final, 2, colSums(otu_mat_final), `/`) * 100
otu_mat_145_relative <- sweep(otu_mat_145_final, 2, colSums(otu_mat_145_final), `/`) * 100

# Filter taxa based on mean relative abundance (keep only those with >= 1%)
keep_taxa <- rownames(otu_mat_relative)[apply(otu_mat_relative, 1, function(x) any(x >= 1))]
keep_taxa <- keep_taxa[!is.na(keep_taxa)]

keep_taxa_145 <- rownames(otu_mat_145_relative)[apply(otu_mat_145_relative, 1, function(x) any(x >= 1))]
keep_taxa_145 <- keep_taxa_145[!is.na(keep_taxa_145)]

otu_mat_filtered_counts <- otu_mat_final[keep_taxa, ]
otu_mat_145_filtered_counts <- otu_mat_145_final[keep_taxa_145, ]


# Filter taxonomy again based on newly calculated rel abundance in the OTU matrix
tax_mat_filtered_counts <- tax_mat_filtered[rownames(tax_mat_filtered) %in% rownames(otu_mat_filtered_counts), ]
tax_mat_145_filtered_counts <- tax_mat_145_filtered[rownames(tax_mat_145_filtered) %in% rownames(otu_mat_145_filtered_counts), ]
# Phyloseq Object Creation --------------------------------------------------
rownames(samples_df_filtered) <- make.unique(as.character(samples_df_filtered$sample.id))


# Step 12: Create phyloseq object
rownames(samples_df_filtered) <- samples_df_filtered$sample.id
samples_df_filtered$sample.id <- NULL  # optional: drop column after it's used as rowname
samples <- sample_data(samples_df_filtered)



OTU <- otu_table(otu_mat_filtered_counts, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat_filtered_counts)
samples <- sample_data(samples_df_filtered)
physeq <- phyloseq(OTU, TAX, samples)

#i know i overwrote the previous physeq
OTU <- otu_table(otu_mat_145_filtered_counts, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat_145_filtered_counts)
samples <- sample_data(samples_df_filtered)
physeq <- phyloseq(OTU, TAX, samples)


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
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = lineage, shape = lineage)) +
  geom_point(size = 2.5) +
  scale_color_manual(name = "Lineage",
                     values = c("ancestral" = "#1f78b4", "derived" = "#e31a1c")) +
  scale_shape_manual(name = "percent_other",
                     values = c("50" = 16, "75" = 17, )) + # Different shapes for ancestry groups
  labs(x = paste("PC1 (", round(varexpl[1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(varexpl[2], 2), "%)", sep = "")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme_cowplot()



pca_plot 

#save plot
ggsave("Expt1_pca_plot.jpg", pca_plot, width = 8, height = 5)
ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/outputs/Expt1_pca_plot.jpg", pca_plot, width = 8, height = 5)



# -- Bar Plot of Relative Abundance by Genus ----------------------------------

library(ggplot2)
library(dplyr)
library(phyloseq)
library(tidyr)

# Step 1: Collapse to Genus level
physeq_genus <- tax_glom(physeq, taxrank = "Genus")

# Step 2: Transform to relative abundance
physeq_genus_rel <- transform_sample_counts(physeq_genus, function(x) x / sum(x))

# Step 3: Melt to long format
melted_df <- psmelt(physeq_genus_rel)

# Step 4: Replace NA Genus names with "Unclassified"
melted_df$Genus <- as.character(melted_df$Genus)
melted_df$Genus[is.na(melted_df$Genus)] <- "Unclassified"

# Step 5: Summarize mean relative abundance by Treatment and Genus
plot_df <- melted_df %>%
  group_by(treatment, Genus, treatment_plot) %>%
  summarise(MeanRA = mean(Abundance), .groups = "drop")

unique_genera <- unique(plot_df$Genus)
glasbey_colors <- c(
  "#999999",
  Polychrome::createPalette(length(unique_genera) - 1, seedcolors = c("#0000FF", "#00FF00", "#FF0000"))
)
names(glasbey_colors) <- unique_genera


install.packages("randomcoloR")
library(randomcoloR)

# Generate 30 distinct colors
my_colors <- distinctColorPalette(30)


ggplot(plot_df, aes(x = treatment_plot, y = MeanRA, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  labs(x = "Treatment", y = "Mean Relative Abundance", fill = "Genus") +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95, size = 6),
    legend.position = "right"
  ) +
  scale_fill_manual(values = my_colors)
ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/outputs/barplot.png", width = 18, height = 6, dpi = 300)





