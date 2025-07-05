# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(vegan)
library(tibble)
library(MCMCglmm)
library(phyloseq)
library(tibble)
library(stringr)

# Set working directory
setwd("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution")

# ----- Step 1: Prepare OTU matrix -----
otu_mat <-  read.table("feature_table_no_contam.tsv", header = TRUE, sep = "\t", check.names = FALSE)

cat("Number of ASVs before filtering:", nrow(otu_mat), "\n")

colnames(otu_mat)[1] <- "otu"

rownames(otu_mat) <- NULL  

otu_mat <- column_to_rownames(otu_mat, var = "otu")

otu_mat <- as.matrix(otu_mat)


write.table(otu_mat, "processed_otu_matrix.tsv", quote = FALSE)

# ----- Step 2: Prepare taxonomy matrix -----
tax_mat <- read.table("taxonomy_no_contam.tsv", header = TRUE, sep = "\t")
colnames(tax_mat)[1] <- "otu"
tax_mat <- tax_mat %>%
  column_to_rownames("otu") %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species),
                ~ sub("^.*?__", "", .))) %>%
  as.matrix()

write.table(tax_mat, "processed_taxonomy_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)


# Relative abundance for filtering

otu_mat_relative <- sweep(otu_mat, 2, colSums(otu_mat), `/`) * 100

# Keep taxa with relative abundance 
keep_taxa <- rownames(otu_mat_relative)[apply(otu_mat_relative, 1, function(x) any(x >= .1))]
keep_taxa <- keep_taxa[!is.na(keep_taxa)]


otu_filtered <- otu_mat[keep_taxa, ]
tax_filtered <- tax_mat[keep_taxa, ]


cat("Number of ASVs after relative abundance filtering:", nrow(otu_filtered), "\n")

write.table(otu_filtered, "filtered_otu_matrix.tsv", sep = "\t", quote = FALSE)
write.table(tax_filtered, "filtered_taxonomy_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)

# ----- Step 3: Prepare sample metadata -----
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t")
write.table(metadata, "metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

samples_df <- read.table("metadata.tsv", header = TRUE)
colnames(samples_df)[colnames(samples_df) == "sampleid"] <- "sample"
samples_df$sample <- rownames(samples_df)
write.table(samples_df, "processed_sample_metadata.tsv", quote = FALSE)

# ----- Load filtered matrices -----
otu_mat <- as.matrix(read.table("filtered_otu_matrix.tsv", header = TRUE, row.names = 1, check.names = FALSE))
tax_mat <- as.matrix(read.table("filtered_taxonomy_matrix.tsv", header = TRUE, row.names = 1, sep = "\t"))



# ----- Load additional metadata -----
samples_df <- read.csv("samples_df.csv")
samples_df <- samples_df %>%
  mutate(lineage = case_when(
    grepl("^ancestral", treatment) ~ "ancestral",
    grepl("^derived", treatment) ~ "derived",
    grepl("^blank", treatment) ~ "blank",
    TRUE ~ NA_character_
  )) %>%
  filter(lineage %in% c("ancestral", "derived"))


print(colnames(samples_df))




# Subset OTU matrix to selected samples
otu_subset <- otu_mat[, samples_df$sample.id, drop = FALSE]
print(head(rownames(otu_subset)))  # Should not be NULL


# Now subset taxonomy matrix
common_otus <- intersect(rownames(otu_subset), rownames(tax_mat))
tax_subset <- tax_mat[common_otus, , drop = FALSE]
print(head(rownames(tax_subset)))


# Keep only taxa present in both otu_subset and tax_subset
keep_taxa <- intersect(keep_taxa, rownames(otu_subset))
keep_taxa <- intersect(keep_taxa, rownames(tax_subset))


# Subset OTU and taxonomy matrices using filtered keep_taxa
otu_final <- otu_subset[keep_taxa, , drop = FALSE]
tax_final <- tax_subset[keep_taxa, , drop = FALSE]

cat("Number of ASVs otu_final:", nrow(otu_final), "\n")

# Match sample IDs

# Step 1: Find the matching sample IDs between OTU matrix columns and samples_df sample.id column
matched_ids <- intersect(colnames(otu_final), samples_df$sample.id)

samples_df <- samples_df %>%
  filter(!duplicated(sample.id))

print(colnames(samples_df))


# Step 2: Filter samples_df by matched_ids and arrange to match order
samples_df <- samples_df %>%
  filter(sample.id %in% matched_ids) %>%
  arrange(match(sample.id, matched_ids))

print(colnames(samples_df))

# Step 3: Set sample.id as rownames (only now)
samples_df <- samples_df %>% column_to_rownames("sample.id")

print(colnames(samples_df))

# Step 4: Subset otu_final columns using matched_ids
otu_final <- otu_final[, matched_ids, drop = FALSE]

# Check results
print(head(rownames(samples_df)))  # Should print matched_ids
print(head(colnames(otu_final)))   # Should print matched_ids, matching sample_df rownames
print(head(colnames(samples_df)))  # Should print matched_ids, matching otu_final columns)))

# ----- Create phyloseq object -----
OTU <- otu_table(otu_final, taxa_are_rows = TRUE)
TAX <- tax_table(tax_final)
samples <- sample_data(samples_df)
print(head(colnames(samples_df)))

physeq <- phyloseq(OTU, TAX, samples)

# Save phyloseq object
saveRDS(physeq, "C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/phylogeny-data-analysis/physeq.rds")



#_______________________________________________________________________________


#______________creation of final experiment phyloseq object_____________________

#read csv
enddat <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/final_experiment/enddat.csv")


#Preparing enddat to be used in phyloseq object

enddat <- enddat %>%
  filter(str_detect(treatment, "^Start|^End")) %>%  # Keep only Start* or End* entries
  mutate(
    sample.id = case_when(
      str_starts(treatment, "Start") ~ str_replace(treatment, "^Start", "ancestral-"),
      str_starts(treatment, "End") ~ str_replace(treatment, "^End", "derived-")
    )
  )


# Match sample IDs

# Step 1: Find matching sample IDs between OTU matrix columns and samples_df sample.id column
matched_ids <- intersect(colnames(otu_final), enddat$sample.id)

# Step 2: Subset enddat rows using matched_ids
enddat_filtered <- enddat[enddat$sample.id %in% matched_ids, , drop = FALSE]
enddat_filtered <- enddat_filtered[!duplicated(enddat_filtered$sample.id), ]

# Set rownames of enddat_filtered to sample.id
rownames(enddat_filtered) <- enddat_filtered$sample.id

# Step 3: Subset otu_final columns using matched_ids
otu_final <- otu_final[, matched_ids, drop = FALSE]

# Check results
print(head(rownames(enddat_filtered)))  # Now should print the matched sample IDs
print(head(colnames(otu_final)))        # Should be the same sample IDs




# ----- Create phyloseq object -----
OTU <- otu_table(otu_final, taxa_are_rows = TRUE)
TAX <- tax_table(tax_final)
samples <- sample_data(enddat_filtered)

final_physeq <- phyloseq(OTU, TAX, samples)



#save phyloseq object
saveRDS(final_physeq, "C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/phylogeny-data-analysis/final_physeq.rds")



#Filtering amounts
cat("Step 0 - Original number of taxa (ASVs):", nrow(otu_mat), "\n")

cat("Step 1 - Taxa retained after >0.1% relative abundance filter:", nrow(otu_filtered), "\n")
cat("         Taxa removed:", nrow(otu_mat) - nrow(otu_filtered), "\n")


cat("tax_final length", nrow(tax_final), "\n")


#How many asvs in ancestral vs derived

# Extract OTU matrix
otu_mat <- as.matrix(otu_table(physeq))

# Extract sample metadata
sample_data_df <- data.frame(sample_data(physeq))

# Identify samples by lineage
ancestral_samples <- rownames(sample_data_df)[sample_data_df$lineage == "ancestral"]
derived_samples <- rownames(sample_data_df)[sample_data_df$lineage == "derived"]

# Subset OTU matrix for each lineage
otu_ancestral <- otu_mat[, ancestral_samples, drop = FALSE]
otu_derived <- otu_mat[, derived_samples, drop = FALSE]

# Count number of ASVs (rows with reads > 0) per sample (columns)
asvs_per_sample_ancestral <- colSums(otu_ancestral > 0)
asvs_per_sample_derived <- colSums(otu_derived > 0)

# Get range (min and max) of ASVs per sample in each lineage
range_ancestral <- range(asvs_per_sample_ancestral)
range_derived <- range(asvs_per_sample_derived)
cat("Range of ASVs per sample in ancestral lineage:", range_ancestral, "\n")
cat("Range of ASVs per sample in derived lineage:", range_derived, "\n")

mean_ancestral <- mean(asvs_per_sample_ancestral)
se_mean_ancestral <- sd(asvs_per_sample_ancestral) / sqrt(length(asvs_per_sample_ancestral))
cat("Mean number of ASVs per sample in ancestral lineage:", mean_ancestral, "\n")
cat("Standard error of mean ASVs in ancestral lineage:", se_mean_ancestral, "\n")

mean_derived <- mean(asvs_per_sample_derived)
se_mean_derived <- sd(asvs_per_sample_derived) / sqrt(length(asvs_per_sample_derived))
cat("Mean number of ASVs per sample in derived lineage:", mean_derived, "\n")
cat("Standard error of mean ASVs in derived lineage:", se_mean_derived, "\n")



