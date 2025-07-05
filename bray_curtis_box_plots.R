# 1. Calculate Bray-Curtis dissimilarity on relative abundance data
library(vegan)
library(ggplot2)
library(phyloseq)
library(cowplot)

#FLAGGED - another table: within derived same host vs different host, within same treatment vs different treatment

#import physeq
physeq <- readRDS("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/phylogeny-data-analysis/physeq.rds")

bray_dist <- distance(physeq, method = "bray")
braymat <- as.matrix(bray_dist)


#___________________________________________________________________________________________________________
# Goal 1: one plot comparing the avg BC distance between the ancestral and derived communities, 
# the avg BC distance between ancestral and derived within the same host source (column "other_host_source" needs to match),
# when all treatments match (column "title_of_graph" needs to match)



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

# #mcmcglmm
# modtmp_plotdf <- MCMCglmm(Distance ~ Group -1 , data = plot_df, nitt = 10000, burnin = 5000, thin = 10,verbose=F)
# #modtmp <- MCMCglmm(Distance ~ Group , data = plot_within_df, nitt = 10000, burnin = 5000, thin = 10,verbose=F)
# summary(modtmp_plotdf)






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

test_df<- bind_rows(df_within_anc, df_within_der)
test <- MCMCglmm(Distance ~ Group , data = test_df, nitt = 10000, burnin = 5000, thin = 10,verbose=F)
summary(test)


ggplot(plot_within_df, aes(x = Group, y = Distance)) +
  geom_boxplot() +
  theme_cowplot() +
  labs(title = "Bray-Curtis Distance Within Groups", y = "BC Distance")

ggsave("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/outputs/bc_within_boxplot.png", width = 6, height = 8, bg = "white")

#statistical tests

modtmpnointercept <- MCMCglmm(Distance ~ Group -1 , data = plot_within_df, nitt = 1000, burnin = 500, thin = 10,verbose=F)
summary(modtmpnointercept)


modtmp <- MCMCglmm(Distance ~ Group , data = plot_within_df, nitt = 10000, burnin = 5000, thin = 10,verbose=F)
summary(modtmp)
tapply(plot_within_df$Distance, plot_within_df$Group, mean)
