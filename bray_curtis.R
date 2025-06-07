# 1. Calculate Bray-Curtis dissimilarity on relative abundance data
library(vegan)
library(ggplot2)
# install.packages("phyloseq")
library(phyloseq)
library(cowplot)

bray_dist <- distance(physeq, method = "bray")
braymat <- as.matrix(bray_dist)

#_______________________________________________________________________________
# side mission here:
image(braymat)

# Calculate Bray-Curtis distance and convert to matrix
bray_dist <- distance(physeq, method = "bray")
braymat <- as.matrix(bray_dist)

# Set up the image plot with labeled axes
image(
  1:nrow(braymat), 1:ncol(braymat), braymat,
  xaxt = "n", yaxt = "n", xlab = "", ylab = "",
  main = "Bray-Curtis Dissimilarity Matrix"
)

# Add axis labels
axis(1, at = 1:ncol(braymat), labels = colnames(braymat), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(braymat), labels = rownames(braymat), las = 2, cex.axis = 0.7)

#_______________________________________________________________________________

ordinate <- ordinate(physeq, method = "PCoA", distance = bray_dist)

plot_ordination(physeq, ordinate, color = "lineage") + 
  geom_point(size = 2) +
  theme_cowplot() 
  
#theme(legend.position = "none")
plot_ordination(physeq, ordinate, color = "treatment_plot") + 
  geom_point(size = 2) +
  theme_cowplot() #log transformed? 

#make a heatmap
library(ggplotify)
library(ggdendro)
library(ggdendro)
library(ggplot2)
library(reshape2)


#___________________ messing with bc stuff more______________________
sample_data(physeq)$title_of_graph  # Check variable

physeq_sub <- subset_samples(physeq, title_of_graph %in% c("50_riccia", "100_lemna"))

# Remove taxa not present in the subset
physeq_sub <- prune_taxa(taxa_sums(physeq_sub) > 0, physeq_sub)


# NMDS ordination using Bray-Curtis
ordu <- ordinate(physeq_sub, method = "NMDS", distance = "bray")

# Plot NMDS with treatment labels
plot_ordination(physeq_sub, ordu, shape = "title_of_graph", color = "treatment_plot") +
  geom_point(size = 4) +
  theme_minimal() +
  ggtitle("NMDS Ordination - 50% Riccia vs 100% Lemna (Bray-Curtis)")




physeq_sub2 <- subset_samples(physeq, title_of_graph %in% c("50_riccia", "100_lemna", "75_riccia"))

# Remove taxa not present in the subset
physeq_sub2 <- prune_taxa(taxa_sums(physeq_sub2) > 0, physeq_sub2)


# NMDS ordination using Bray-Curtis
ordu <- ordinate(physeq_sub2, method = "NMDS", distance = "bray")

# Plot NMDS with treatment labels
plot_ordination(physeq_sub2, ordu, shape = "title_of_graph", color = "treatment_plot") +
  geom_point(size = 4) +
  theme_minimal() +
  ggtitle("NMDS Ordination - 50% Riccia, 75% Riccia, and 100% Lemna (Bray-Curtis)")


#__________trying to compare bc distance and rgr__________
# Subset only derived communities
phy_bc <- subset_samples(physeq, lineage == "derived")

# Prune taxa not present in the subset
phy_bc <- prune_taxa(taxa_sums(phy_bc) > 0, phy_bc)

# Include 100% Lemna samples for comparison
phy_bc <- merge_phyloseq(
  phy_bc,
  subset_samples(physeq, title_of_graph == "100_lemna")
)

# Re-prune taxa after merge
phy_bc <- prune_taxa(taxa_sums(phy_bc) > 0, phy_bc)

# Compute Bray-Curtis distance matrix
bc_dist <- phyloseq::distance(phy_bc, method = "bray")
bc_mat <- as.matrix(bc_dist)

# Get sample names of 100% Lemna
lemna_samples <- sample_names(subset_samples(phy_bc, title_of_graph == "100_lemna"))

# Compute mean Bray-Curtis distance to 100% Lemna for each sample
bc_to_lemna <- rowMeans(bc_mat[, lemna_samples, drop = FALSE])

# Assign NA to Lemna samples (optional cleanup)
bc_to_lemna[lemna_samples] <- NA

# Add distances to sample_data
sample_data(phy_bc)$BC_to_lemna <- bc_to_lemna[rownames(sample_data(phy_bc))]

bcLemnaAnc <- bc_mat["ancestral-1",]  #
bcLemnaAnc["ancestral-1"] <- NA  # Set the distance to itself as NA
sample_data(phy_bc)$BC_to_lemnaanc <- bcLemnaAnc
  
library(ggplot2)

# Convert sample_data to data frame for plotting
plot_df <- data.frame(RGR = sample_data(phy_bc)$RGR, 
                         BC_to_lemna = sample_data(phy_bc)$BC_to_lemna,
                         BC_to_lemnaanc = sample_data(phy_bc)$BC_to_lemnaanc,
                         title_of_graph = sample_data(phy_bc)$title_of_graph)

ggplot(plot_df, aes(x = BC_to_lemna, y = RGR, color = title_of_graph)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray40") +
  geom_text(aes(label = title_of_graph), vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(title = "RGR vs Bray-Curtis Distance to 100% Lemna",
       x = "Bray-Curtis Distance to 100% Lemna",
       y = "Relative Growth Rate (RGR)")
library(MCMCglmm)

ggplot(plot_df, aes(x = BC_to_lemnaanc, y = RGR, color = title_of_graph)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray40") +
  geom_text(aes(label = title_of_graph), vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(title = "RGR vs Bray-Curtis Distance to 100% Lemna",
       x = "Bray-Curtis Distance to 100% Lemna",
       y = "Relative Growth Rate (RGR)")

plot_df2 <- plot_df[!is.na(plot_df$BC_to_lemnaanc),]
summary(MCMCglmm(RGR ~ BC_to_lemnaanc , data = plot_df2, nitt = 10000, burnin = 5000, thin = 10))

plot_df3 <- plot_df[!is.na(plot_df$BC_to_lemna),]
summary(MCMCglmm(RGR ~ BC_to_lemna , data = plot_df3, nitt = 100000, burnin = 5000, thin = 10))
#ex[;pre to see if there is a positive relationship for some


#more results with bray-curtis




#__________________________trying to compare bc distance and rgr with other_host_type__________

# Define your target host sources
target_hosts <- c("ceratophyllum", "wolffia", "spirodela", "riccia", "elodea", "nymphaea")

# Subset derived communities
phy_derived <- subset_samples(physeq, lineage == "derived")

# Subset target host source samples (and remove NAs)
phy_hosts <- subset_samples(physeq,
                            !is.na(other_host_source) &
                              other_host_source %in% target_hosts)

# Merge derived and host source samples
phy_bc <- merge_phyloseq(phy_derived, phy_hosts)

# Prune taxa
phy_bc <- prune_taxa(taxa_sums(phy_bc) > 0, phy_bc)

# Compute Bray-Curtis distance
bc_dist <- phyloseq::distance(phy_bc, method = "bray")
bc_mat <- as.matrix(bc_dist)

# Get sample names of the reference host source group
ref_samples <- sample_names(phy_hosts)

# Compute mean Bray-Curtis distance to the host group
bc_to_hosts <- rowMeans(bc_mat[, ref_samples, drop = FALSE])

# Optional: set distances to NA for the reference samples themselves
bc_to_hosts[ref_samples] <- NA

# Add to sample_data
sample_data(phy_bc)$BC_to_hosts <- bc_to_hosts[rownames(sample_data(phy_bc))]

sample_data_df <- as.data.frame(sample_data(phy_bc))

# Filter to only derived samples
df_derived <- subset(sample_data_df, lineage == "derived" & !is.na(BC_to_hosts))

ggplot(df_derived, aes(x = BC_to_hosts, y = RGR, color = title_of_graph)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "gray40", linetype = "dashed") +
  theme_minimal() +
  labs(title = "RGR vs Bray-Curtis Distance to Host-Derived Communities",
       x = "Bray-Curtis Distance to Host Source Communities",
       y = "Relative Growth Rate (RGR)")







#___________________________________________________________________________________________________________
# Goal 1: one plot comparing the avg BC distance between the ancestral and derived communities, 
# the avg BC distance between ancestral and derived within the same host source (column "other_host_source" needs to match),
# when all treatments match (column "title_of_graph" needs to match)


# 1. Calculate Bray-Curtis distances
bray_dist <- distance(physeq, method = "bray")
braymat <- as.matrix(bray_dist)


get_avg_bc_between_groups <- function(dist_mat, meta, group1_label, group2_label, match_col = NULL) {
  samples <- rownames(meta)
  dists <- c()
  count <- 0
  for (i in 1:(length(samples) - 1)) {
    for (j in (i + 1):length(samples)) {
      s1 <- samples[i]
      s2 <- samples[j]
      
      cond1 <- (meta[s1, "lineage"] == group1_label && meta[s2, "lineage"] == group2_label)
      cond2 <- (meta[s1, "lineage"] == group2_label && meta[s2, "lineage"] == group1_label)
      lineage_match <- cond1 || cond2
      
      match_group <- TRUE
      if (!is.null(match_col)) {
        match_group <- !is.na(meta[s1, match_col]) &&
          !is.na(meta[s2, match_col]) &&
          meta[s1, match_col] == meta[s2, match_col]
      }
      
      if (lineage_match && match_group) {
        dists <- c(dists, dist_mat[s1, s2])
        count <- count + 1
      }
    }
  }
  message("Number of pairs: ", count)
  return(mean(dists, na.rm = TRUE))
}


### 2. Subset and prune
anc_derived <- subset_samples(physeq, lineage %in% c("ancestral", "derived"))
anc_derived <- prune_taxa(taxa_sums(anc_derived) > 0, anc_derived)
anc_derived_dist <- distance(anc_derived, method = "bray")
anc_derived_mat <- as.matrix(anc_derived_dist)
meta <- as(sample_data(anc_derived), "data.frame")

# Avg Bray-Curtis: Ancestral vs Derived
avg_anc_derived <- get_avg_bc_between_groups(anc_derived_mat, meta, "ancestral", "derived")

# Avg Bray-Curtis: Ancestral vs Derived (Same Host)
avg_anc_derived_same_host <- get_avg_bc_between_groups(anc_derived_mat, meta, "ancestral", "derived", match_col = "other_host_source")

# Avg Bray-Curtis: Ancestral vs Derived (Same Treatment)
avg_anc_derived_same_treatment <- get_avg_bc_between_groups(anc_derived_mat, meta, "ancestral", "derived", match_col = "title_of_graph")

### 3. Plotting data
avg_bc_df <- data.frame(
  Comparison = c("Ancestral vs Derived", 
                 "Ancestral vs Derived (Same Host)", 
                 "Ancestral vs Derived (Same Treatment)"),
  Avg_BC = c(avg_anc_derived, avg_anc_derived_same_host, avg_anc_derived_same_treatment)
)

### 4. Plot
ggplot(avg_bc_df, aes(x = Comparison, y = Avg_BC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Average Bray-Curtis Dissimilarity",
       x = "Comparison",
       y = "Average Bray-Curtis Dissimilarity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("avg_bc_plot.png", width = 8, height = 6)

library(dplyr)


#visualize how different the distributions are - the figure i saved
get_bc_pairs_df <- function(dist_mat, meta, group1_label, group2_label, match_col = NULL, label = "Comparison") {
  samples <- rownames(meta)
  pairs <- list()
  for (i in 1:(length(samples) - 1)) {
    for (j in (i + 1):length(samples)) {
      s1 <- samples[i]
      s2 <- samples[j]
      
      cond1 <- (meta[s1, "lineage"] == group1_label && meta[s2, "lineage"] == group2_label)
      cond2 <- (meta[s1, "lineage"] == group2_label && meta[s2, "lineage"] == group1_label)
      lineage_match <- cond1 || cond2
      
      match_group <- TRUE
      if (!is.null(match_col)) {
        match_group <- !is.na(meta[s1, match_col]) &&
          !is.na(meta[s2, match_col]) &&
          meta[s1, match_col] == meta[s2, match_col]
      }
      
      if (lineage_match && match_group) {
        pairs <- append(pairs, list(data.frame(
          Sample1 = s1,
          Sample2 = s2,
          Distance = dist_mat[s1, s2],
          Group = label
        )))
      }
    }
  }
  do.call(rbind, pairs)
}

df_all <- get_bc_pairs_df(anc_derived_mat, meta, "ancestral", "derived", label = "Ancestral vs. Derived")
df_host <- get_bc_pairs_df(anc_derived_mat, meta, "ancestral", "derived", match_col = "other_host_source", label = "Same Host")
df_treat <- get_bc_pairs_df(anc_derived_mat, meta, "ancestral", "derived", match_col = "title_of_graph", label = "Same Treatment")

plot_df <- bind_rows(df_all, df_host, df_treat)

ggplot(plot_df, aes(x = Group, y = Distance)) +
  geom_boxplot() +
  theme_cowplot() +
  labs(title = "Bray-Curtis Distance by Group", y = "BC Distance")

ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/outputs/avg_bc_plot.png", width = 8, height = 6)



# Goal 2: one plot with the avg BC within Ancestral, avg BC within derived, avg BC within derived same host, 
# avg BC within derived same treatment (title_of_graph matching)

#step 1: 
get_avg_bc_within_group <- function(dist_mat, meta, group_label, match_col = NULL) {
  samples <- rownames(meta)
  dists <- c()
  count <- 0
  for (i in 1:(length(samples) - 1)) {
    for (j in (i + 1):length(samples)) {
      s1 <- samples[i]
      s2 <- samples[j]
      
      if (meta[s1, "lineage"] == group_label && meta[s2, "lineage"] == group_label) {
        match_group <- TRUE
        if (!is.null(match_col)) {
          match_group <- !is.na(meta[s1, match_col]) &&
            !is.na(meta[s2, match_col]) &&
            meta[s1, match_col] == meta[s2, match_col]
        }
        
        if (match_group) {
          dists <- c(dists, dist_mat[s1, s2])
          count <- count + 1
        }
      }
    }
  }
  message("Number of pairs: ", count)
  return(mean(dists, na.rm = TRUE))
}

#calculating the averages
avg_within_anc <- get_avg_bc_within_group(anc_derived_mat, meta, "ancestral")
avg_within_der <- get_avg_bc_within_group(anc_derived_mat, meta, "derived")
avg_within_der_same_host <- get_avg_bc_within_group(anc_derived_mat, meta, "derived", match_col = "other_host_source")
avg_within_der_same_treatment <- get_avg_bc_within_group(anc_derived_mat, meta, "derived", match_col = "title_of_graph")

#making the df so it can be plotted
avg_bc_within_df <- data.frame(
  Comparison = c("Within Ancestral", 
                 "Within Derived", 
                 "Within Derived (Same Host)", 
                 "Within Derived (Same Treatment)"),
  Avg_BC = c(avg_within_anc, avg_within_der, avg_within_der_same_host, avg_within_der_same_treatment)
)

#actually plotting
ggplot(avg_bc_within_df, aes(x = Comparison, y = Avg_BC)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  theme_minimal() +
  labs(title = "Average Bray-Curtis Dissimilarity Within Lineages",
       x = "Comparison",
       y = "Average Bray-Curtis Dissimilarity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


df_within_anc <- get_bc_pairs_df(anc_derived_mat, meta, "ancestral", "ancestral", label = "Within Ancestral")
df_within_der <- get_bc_pairs_df(anc_derived_mat, meta, "derived", "derived", label = "Within Derived")
df_within_der_host <- get_bc_pairs_df(anc_derived_mat, meta, "derived", "derived", match_col = "other_host_source", label = "Within Derived (Same Host)")
df_within_der_treat <- get_bc_pairs_df(anc_derived_mat, meta, "derived", "derived", match_col = "title_of_graph", label = "Within Derived (Same Treatment)")

plot_within_df <- bind_rows(df_within_anc, df_within_der, df_within_der_host, df_within_der_treat)

ggplot(plot_within_df, aes(x = Group, y = Distance)) +
  geom_boxplot() +
  theme_cowplot() +
  labs(title = "Bray-Curtis Distance Within Groups", y = "BC Distance")

ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/outputs/bc_within_boxplot.png", width = 6, height = 8, bg = "white")

#statistical tests
modtmp <- MCMCglmm(Distance ~ Group -1 , data = plot_within_df, nitt = 10000, burnin = 5000, thin = 10,verbose=F)
#modtmp <- MCMCglmm(Distance ~ Group , data = plot_within_df, nitt = 10000, burnin = 5000, thin = 10,verbose=F)
summary(modtmp)
