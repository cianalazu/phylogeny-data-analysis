#This script is meant to see if any genus level taxa are significantly correlated with RGR

library(broom)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(tibble)
library(MCMCglmm)

final_physeq <- readRDS("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/phylogeny-data-analysis/final_physeq.rds")


physeq_genus <- tax_glom(final_physeq, taxrank = "Genus")
physeq_rel <- transform_sample_counts(physeq_genus, function(x) x / sum(x))
df <- psmelt(physeq_rel)

results <- df |>
  filter(!is.na(Abundance), !is.infinite(Abundance), !is.na(RGR)) |>
  group_by(Genus) |>
  do(tidy(lm(Abundance ~ RGR, data = .))) |>
  filter(term == "RGR") |>
  arrange(p.value)

results |>  filter(p.value < 0.05)
# Plot the significant results
ggplot(results, aes(x = reorder(Genus, p.value), y = estimate)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Genus", y = "Estimate of RGR effect", title = "Significant Genus-level Correlations with RGR") +
  theme_minimal()

ggplot(filter(df, Genus == "Pseudomonas"), aes(x = RGR, y = Abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Pseudomonas Abundance vs. RGR")


ggplot(filter(df, Genus == "Chryseobacterium"), aes(x = RGR, y = Abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Chryseobacterium Abundance vs. RGR")

ggplot(filter(df, Genus == "Stenotrophomonas"), aes(x = RGR, y = Abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Stenotrophomonas Abundance vs. RGR")

#_____________basic graphing of some of the significant taxa___________________________


# Filter the data
df_pseudo <- filter(df, Genus == "Pseudomonas")

# Fit the model
model <- lm(Abundance ~ RGR, data = df_pseudo)
pval <- summary(model)$coefficients["RGR", "Pr(>|t|)"]

# Format p-value for display
pval_label <- paste0("p = ", signif(pval, 3))

# Plot with p-value
ggplot(df_pseudo, aes(x = RGR, y = Abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = pval_label, hjust = 1.1, vjust = 1.5, size = 5) +
  labs(title = "Pseudomonas Abundance vs. RGR")



# Filter the data
df_Stenotrophomonas <- filter(df, Genus == "Stenotrophomonas")

# Fit the model
model <- lm(Abundance ~ RGR, data = df_Stenotrophomonas)
pval <- summary(model)$coefficients["RGR", "Pr(>|t|)"]

# Format p-value for display
pval_label <- paste0("p = ", signif(pval, 3))

# Plot with p-value
ggplot(df_Stenotrophomonas, aes(x = RGR, y = Abundance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = pval_label, hjust = 1.1, vjust = 1.5, size = 5) +
  labs(title = "Pseudomonas Abundance vs. RGR")



#______________________mcmcglmm to see if Stenotrophomonas is significant_________________________________
# Load necessary packages
library(phyloseq)
library(MCMCglmm)
library(tidyverse)


# Sample metadata
meta <- data.frame(sample_data(physeq_end))

# Add Stenotrophomonas abundance (assuming you've already extracted `steno_abundance`)
steno_taxa <- taxa_names(physeq_end)[
  tax_table(physeq_end)[, "Genus"] == "Stenotrophomonas"
]
steno_taxa <- steno_taxa[!is.na(steno_taxa)]
steno_abundance <- otu_table(physeq_end)[steno_taxa, ]

# Collapse multiple ASVs if needed
steno_counts <- if (length(steno_taxa) > 1) {
  colSums(steno_abundance)
} else {
  as.numeric(steno_abundance)
}

# Add to metadata
meta$Steno <- steno_counts[rownames(meta)]

model <- MCMCglmm(
  fixed = RGR ~ Steno,
  data = meta,
  nitt = 13000, burnin = 3000, thin = 10,
  verbose = FALSE
)
# Check model summary
summary(model)



