# Load libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(cowplot)
library(ggplotify)
library(ggdendro)
library(reshape2)
library(MCMCglmm)
library(dplyr)

#-----------------------------#
# Bray-Curtis Dissimilarity
#-----------------------------#
bray_dist <- distance(physeq, method = "bray")
braymat <- as.matrix(bray_dist)

# Optional visualization
image(1:nrow(braymat), 1:ncol(braymat), braymat,
      xaxt = "n", yaxt = "n", xlab = "", ylab = "",
      main = "Bray-Curtis Dissimilarity Matrix")
axis(1, at = 1:ncol(braymat), labels = colnames(braymat), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(braymat), labels = rownames(braymat), las = 2, cex.axis = 0.7)

#-----------------------------#
# PCoA and NMDS Ordinations
#-----------------------------#
ordinate_plot <- function(physeq_obj, method, dist_method = "bray", shape_var = NULL, color_var) {
  ord <- ordinate(physeq_obj, method = method, distance = dist_method)
  plot_ordination(physeq_obj, ord, shape = shape_var, color = color_var) +
    geom_point(size = 3) +
    theme_minimal()
}

# Full PCoA plots
ordinate_plot(physeq, "PCoA", "bray", color_var = "lineage")
ordinate_plot(physeq, "PCoA", "bray", color_var = "treatment_plot")

# NMDS for subsets
physeq_subsets <- list(
  `50vs100` = subset_samples(physeq, title_of_graph %in% c("50_riccia", "100_lemna")),
  `50vs75vs100` = subset_samples(physeq, title_of_graph %in% c("50_riccia", "75_riccia", "100_lemna"))
)

for (label in names(physeq_subsets)) {
  ps <- prune_taxa(taxa_sums(physeq_subsets[[label]]) > 0, physeq_subsets[[label]])
  title <- paste("NMDS Ordination -", gsub("vs", " vs ", label), "(Bray-Curtis)")
  print(ordinate_plot(ps, "NMDS", "bray", shape_var = "title_of_graph", color_var = "treatment_plot") +
          ggtitle(title))
}

#-----------------------------#
# BC Distance vs RGR
#-----------------------------#
calculate_and_plot_bc_rgr <- function(phy_obj, reference_label, ref_column, ancestor_id, y_label = "RGR") {
  ref_samples <- sample_names(subset_samples(phy_obj, !!as.name(ref_column) == reference_label))
  bc_mat <- as.matrix(distance(phy_obj, method = "bray"))
  
  bc_to_ref <- rowMeans(bc_mat[, ref_samples, drop = FALSE])
  bc_to_ref[ref_samples] <- NA
  sample_data(phy_obj)$BC_to_ref <- bc_to_ref[rownames(sample_data(phy_obj))]
  
  bc_to_ancestor <- bc_mat[ancestor_id, ]
  bc_to_ancestor[ancestor_id] <- NA
  sample_data(phy_obj)$BC_to_ancestor <- bc_to_ancestor[rownames(sample_data(phy_obj))]
  
  plot_df <- data.frame(
    y_value = sample_data(phy_obj)[[y_label]],
    BC_to_ref = sample_data(phy_obj)$BC_to_ref,
    BC_to_ancestor = sample_data(phy_obj)$BC_to_ancestor,
    title_of_graph = sample_data(phy_obj)$title_of_graph
  )
  
  ggplot(plot_df, aes(x = BC_to_ref, y = y_value, color = title_of_graph)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray40") +
    geom_text(aes(label = title_of_graph), vjust = -0.5, size = 3) +
    theme_minimal() +
    labs(title = paste(y_label, "vs Bray-Curtis Distance to", reference_label),
         x = paste("Bray-Curtis Distance to", reference_label),
         y = y_label)
}
