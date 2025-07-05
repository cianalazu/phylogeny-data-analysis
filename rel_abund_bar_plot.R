library(ggplot2)
library(phyloseq)
library(cowplot)
library(Polychrome)
library(randomcoloR)

#what is getting dropped from the melted_df/plot_df

physeq <- readRDS("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/phylogeny-data-analysis/physeq.rds")
# -- Bar Plot of Relative Abundance by Genus ----------------------------------


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
  group_by(percent_other,other_host_source, lineage, Genus, treatment_plot) %>%
  summarise(MeanRA = mean(Abundance), .groups = "drop")

unique_genera <- unique(plot_df$Genus)
glasbey_colors <- c(
  "#999999",
  Polychrome::createPalette(length(unique_genera) - 1, seedcolors = c("#0000FF", "#00FF00", "#FF0000"))
)
names(glasbey_colors) <- unique_genera


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
ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/outputs/barplot.png", width = 20, height = 6, dpi = 300)



# attempt to get barplot more organized

unique_genera <- unique(plot_df$Genus)
glasbey_colors <- c(
  "#999999",
  Polychrome::createPalette(length(unique_genera) - 1, seedcolors = c("#0000FF", "#00FF00", "#FF0000"))
)
names(glasbey_colors) <- unique_genera


# Generate 30 distinct colors
my_colors <- distinctColorPalette(30)


# Step: Reorder Genus by total relative abundance
genus_order <- plot_df %>%
  group_by(Genus) %>%
  summarise(TotalRA = sum(MeanRA, na.rm = TRUE)) %>%
  arrange(desc(TotalRA)) %>%
  pull(Genus)

plot_df$Genus <- factor(plot_df$Genus, levels = genus_order)

plot_df <- plot_df %>%
  arrange(lineage, other_host_source, percent_other) %>%
  mutate(
    treatment_plot = factor(treatment_plot, levels = unique(treatment_plot)),
    Genus = factor(Genus, levels = genus_order)
  )


# Step: Generate colors in the same order
# glasbey_colors <- Polychrome::createPalette(length(genus_order), seedcolors = c("#0000FF", "#00FF00", "#FF0000"))
# names(glasbey_colors) <- genus_order

# Step: Plot
ggplot(plot_df, aes(x = treatment_plot, y = MeanRA, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  scale_fill_manual(values = glasbey_colors) +
  labs(x = "Treatment", y = "Mean Relative Abundance", fill = "Genus") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95, size = 6),
    legend.position = "right"
  )


#_____________genus avg +/- s.e.___________


# Helper function to compute mean and SE safely
summary_stats <- function(x) {
  x <- as.numeric(x)
  tibble(
    mean_ra = mean(x, na.rm = TRUE),
    se_ra = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  )
}

# Ensure all combinations of Genus × Lineage are present
melted_df$Genus[is.na(melted_df$Genus)] <- "Unclassified"

# Group and summarize
genus_change_stats <- melted_df %>%
  group_by(Genus, lineage) %>%
  summarise(
    mean_ra = mean(Abundance, na.rm = TRUE),
    se_ra = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = lineage,
    values_from = c(mean_ra, se_ra),
    values_fill = 0  # <- fill missing combinations with 0
  ) %>%
  mutate(
    change_mean = mean_ra_derived - mean_ra_ancestral,
    change_se = sqrt(se_ra_derived^2 + se_ra_ancestral^2)
  )



#____________How many of these taxa were lost in the derived communities_________________
# Split phyloseq object by lineage
phy_ancestral <- subset_samples(physeq, lineage == "ancestral")
phy_derived <- subset_samples(physeq, lineage == "derived")

# Prune taxa not present in each subset
phy_ancestral <- prune_taxa(taxa_sums(phy_ancestral) > 0, phy_ancestral)
phy_derived <- prune_taxa(taxa_sums(phy_derived) > 0, phy_derived)

# Get taxa names
ancestral_taxa <- taxa_names(phy_ancestral)
derived_taxa <- taxa_names(phy_derived)

# Calculate lost and gained
lost_taxa <- setdiff(ancestral_taxa, derived_taxa)
gained_taxa <- setdiff(derived_taxa, ancestral_taxa)
shared_taxa <- intersect(ancestral_taxa, derived_taxa)

# Summarize
cat("Number of taxa in ancestral:", length(ancestral_taxa), "\n")
cat("Number of taxa in derived:", length(derived_taxa), "\n")
cat("Number of shared taxa:", length(shared_taxa), "\n")
cat("Number of lost taxa (present in ancestral, absent in derived):", length(lost_taxa), "\n")
cat("Number of gained taxa (absent in ancestral, present in derived):", length(gained_taxa), "\n")


# Collapse to genus level
physeq_genus <- tax_glom(physeq, taxrank = "Genus")

# Then subset by lineage
phy_ancestral_genus <- subset_samples(physeq_genus, lineage == "ancestral")
phy_derived_genus <- subset_samples(physeq_genus, lineage == "derived")

# Prune zero-sum taxa (i.e., remove absent genera)
phy_ancestral_genus <- prune_taxa(taxa_sums(phy_ancestral_genus) > 0, phy_ancestral_genus)
phy_derived_genus <- prune_taxa(taxa_sums(phy_derived_genus) > 0, phy_derived_genus)

# Get Genus names (row names of tax_table)
ancestral_genera <- as.character(tax_table(phy_ancestral_genus)[, "Genus"])
derived_genera <- as.character(tax_table(phy_derived_genus)[, "Genus"])

# Remove NAs if needed
ancestral_genera <- ancestral_genera[!is.na(ancestral_genera)]
derived_genera <- derived_genera[!is.na(derived_genera)]

# Compare
lost_genera <- setdiff(ancestral_genera, derived_genera)
gained_genera <- setdiff(derived_genera, ancestral_genera)
shared_genera <- intersect(ancestral_genera, derived_genera)

# Summary
cat("Number of genera in ancestral:", length(unique(ancestral_genera)), "\n")
cat("Number of genera in derived:", length(unique(derived_genera)), "\n")
cat("Number of shared genera:", length(shared_genera), "\n")
cat("Number of lost genera:", length(lost_genera), "\n")
cat("Number of gained genera:", length(gained_genera), "\n")


# Output lost genera
cat("Number of lost genera from ancestral to derived:", length(lost_genera), "\n")
cat("Lost genera:\n")
print(lost_genera)



#_________ how much did Pseudomonadaceae increase in derived communities________

physeq_family <- tax_glom(physeq, taxrank = "Family")

physeq_family_rel <- transform_sample_counts(physeq_family, function(x) x / sum(x))

melted_family_df <- psmelt(physeq_family_rel)

pseudomonadaceae_df <- melted_family_df %>%
  filter(Family == "Pseudomonadaceae")

pseudo_summary <- pseudomonadaceae_df %>%
  group_by(lineage) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE))


derived_abundance <- pseudo_summary$mean_abundance[pseudo_summary$lineage == "derived"]
ancestral_abundance <- pseudo_summary$mean_abundance[pseudo_summary$lineage == "ancestral"]

percent_increase <- ((derived_abundance - ancestral_abundance) / ancestral_abundance) * 100

cat("Percent increase in Pseudomonadaceae in derived communities:", round(percent_increase, 2), "%\n")


