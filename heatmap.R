# Load necessary libraries
library(phyloseq)
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(cowplot)
library(pheatmap)

# Compute Bray-Curtis distance matrix
bray_dist <- distance(physeq, method = "bray")
braymat <- as.matrix(bray_dist)

# Define color palette
palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)

# --------- Option 1: Corrected Base R Image Plot ---------
# Plot the transposed matrix and reverse y-axis for correct orientation
image(1:nrow(braymat), 1:ncol(braymat), t(braymat)[, nrow(braymat):1],
      col = palette, xaxt = 'n', yaxt = 'n', main = "Bray-Curtis Distance Matrix")

# Add axis labels
axis(1, at = 1:ncol(braymat), labels = colnames(braymat), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(braymat), labels = rev(rownames(braymat)), las = 2, cex.axis = 0.7)

# Add color legend
legend("topright", 
       legend = round(seq(min(braymat), max(braymat), length.out = 5), 2),
       fill = palette[seq(1, 100, length.out = 5)],
       title = "Bray-Curtis",
       cex = 0.8)

# --------- Option 2: Cleaner Heatmap with Clustering using pheatmap ---------
pheatmap(braymat,
         color = palette,
         clustering_distance_rows = bray_dist,
         clustering_distance_cols = bray_dist,
         main = "Bray-Curtis Distance Matrix (Clustered)",
         fontsize = 8,
         angle_col = 45)




# trying to use the treatment_plot labels so it's more accurate
# Extract treatment labels
sample_labels <- as.character(sample_data(physeq)$treatment_plot)
names(sample_labels) <- sample_names(physeq)  # Ensure correct mapping

# Reorder the labels to match matrix order
ordered_labels <- sample_labels[rownames(braymat)]

# Plot with treatment labels
hm <-pheatmap(braymat,
         color = palette,
         clustering_distance_rows = bray_dist,
         clustering_distance_cols = bray_dist,
         main = "Bray-Curtis Distance Matrix (Clustered)",
         labels_row = ordered_labels,
         labels_col = ordered_labels,
         fontsize = 7,
         angle_col = 45)

# for this figure -  I need to change:
# sort based on phylogenetic distance then by 50% or 75% - then can not have derived-x-x in the labels

# Save the heatmap to a file
ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/outputs/bray_curtis_heatmap.png", width = 8, height = 5, dpi = 300, plot=hm)



