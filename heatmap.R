# Load necessary libraries
library(phyloseq)
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(scales)

physeq <- readRDS("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/phylogeny-data-analysis/physeq.rds")

# Remove samples with NA in relevant metadata fields
sample_data_df <- as.data.frame(sample_data(physeq))
valid_samples <- complete.cases(sample_data_df$title_of_graph, sample_data_df$treatment_plot)
physeq <- prune_samples(valid_samples, physeq)

# Compute Bray-Curtis distance matrix
bray_dist <- distance(physeq, method = "bray")
braymat <- as.matrix(bray_dist)

# Define color palette
palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)


## --- pull the two different sample‑level variables --------------------------
titles          <- as.character(sample_data(physeq)$other_host_source)   # for colours
treat_labels    <- as.character(sample_data(physeq)$treatment_plot)   # for axis text
names(titles)   <- names(treat_labels) <- sample_names(physeq)

## --- arrange everything in the order of the distance matrix -----------------
ordered_titles  <- titles[rownames(braymat)]
ordered_labels  <- treat_labels[rownames(braymat)]

annotation_df   <- data.frame(Title = ordered_titles)
rownames(annotation_df) <- rownames(braymat)

## --- pick one distinct colour per title_of_graph ---------------------------
title_levels <- unique(ordered_titles)
title_colors <- scales::hue_pal()(length(title_levels))   # scales::hue_pal
names(title_colors) <- title_levels
ann_colors   <- list(Title = title_colors)

## --- draw the heat‑map ------------------------------------------------------
hm <- pheatmap(
  braymat,
  color                   = palette,
  clustering_distance_rows = bray_dist,
  clustering_distance_cols = bray_dist,
  main  = "Bray‑Curtis Distance Matrix (Clustered)",
  labels_row      = ordered_labels,        # text = treatment_plot
  labels_col      = ordered_labels,
  annotation_row  = annotation_df,         # colors = title_of_graph
  annotation_col  = annotation_df,
  annotation_colors = ann_colors,
  fontsize = 20,
  angle_col = 45
)

ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/outputs/bray_curtis_heatmap_hostsource.png", 
       width = 40, height = 50, dpi = 300, plot = hm, bg = "white", limitsize = FALSE)

