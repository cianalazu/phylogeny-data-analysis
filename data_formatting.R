library(ggplot2)

library(dplyr)
library(tidyr)

library(cowplot)

library(vegan)

library(tibble)
library(MCMCglmm)

# Data Preparation Section --------------------------------------------------

#set working directory
setwd("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution")

# Step 1: Prepare OTU matrix
otu_mat <- read.table("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/feature_table_no_contam.tsv", header = TRUE, sep = "\t", check.names = FALSE)
colnames(otu_mat)[colnames(otu_mat) == "OTU_ID"] <- "otu"
rownames(otu_mat) <- NULL            # Remove existing row names
otu_mat <- otu_mat %>% 
    tibble::column_to_rownames(var="otu") %>% 
    as.matrix

# Write out processed OTU matrix
write.table(otu_mat, "processed_otu_matrix.tsv", quote = FALSE)
write.table(otu_mat, "C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/processed_otu_matrix.tsv", quote = FALSE)


# Step 2: Prepare taxonomy matrix
tax_mat <- read.table("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/taxonomy_no_contam.tsv", header = TRUE, sep = "\t")

colnames(tax_mat)[colnames(tax_mat) == "Feature.ID"] <- "otu"
tax_mat <- tax_mat |> 
    tibble::column_to_rownames("otu")

# Split and clean up taxonomy
tax_mat <- tax_mat %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species),
                ~sub("^.*?__", "", .)))

tax_mat <- as.matrix(tax_mat)

# Write out processed taxonomy matrix
write.table(tax_mat, file = "processed_taxonomy_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(tax_mat, file = "C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/processed_taxonomy_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA) 


#interrupted step, filter taxa based on number of reads across all samples
relabund <- sapply(2:ncol(otu_mat), function(z) otu_mat[,z] / sum(otu_mat[,z]))
sum(sapply(1:nrow(relabund), function(z) any(relabund[z,]>0.001,na.rm=T))) #na.rm because some samples have no reads



#only 145 taxa
keeptaxa <- which(sapply(1:nrow(relabund), function(z) any(relabund[z,]>0.001,na.rm=T)))

keepotu_tab <- relabund[keeptaxa,]
keepotu_tax <- tax_mat[keeptaxa,]#otu mat and tax mat are organized identically: sum(otu_mat$otu == tax_mat$otu) == nrow(otu_mat)


# Save filtered OTU matrix
write.table(keepotu_tab, "filtered_otu_matrix.tsv", sep = "\t", quote = FALSE)

# Save filtered taxonomy
write.table(keepotu_tax, "filtered_taxonomy_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)



#barplot at genus level OR at species level, colored by genus.

# Step 3: Prepare sample data
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t") #I only had a txt file of metadata so converted it to a tsv
write.table(metadata, "metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


samples_df <- read.table("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/metadata.tsv", header = TRUE)
colnames(samples_df)[colnames(samples_df) == "sampleid"] <- "sample"
samples_df$sample <- rownames(samples_df)

# Write out processed sample metadata
write.table(samples_df, "processed_sample_metadata.tsv", quote = FALSE)
write.table(samples_df, "C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/processed_sample_metadata.tsv", quote = FALSE)






# Data Preparation Section --------------------------------------------------

# Step 1: Prepare OTU matrix
otu_mat_145 <- read.table("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution/filtered_otu_matrix.tsv", header = TRUE, row.names = 1, check.names = FALSE)
otu_mat_145 <- as.matrix(otu_mat)


# Step 2: Prepare taxonomy matrix
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

otu_mat_145_subset <- otu_mat_145[, samples_df_filtered$sample.id]

otu_mat_145_subset_filtered <- otu_mat_145_subset[rowSums(otu_mat_145_subset) > 0, ]

# Subset the taxonomy matrix

tax_mat_145_filtered <- tax_mat_145[rownames(tax_mat_145) %in% rownames(otu_mat_145_subset_filtered), ]


# Ensure unique taxa names

rownames(tax_mat_145_filtered) <- make.unique(rownames(tax_mat_145_filtered))


otu_mat_145_final <- otu_mat_145_subset_filtered[rownames(otu_mat_145_subset_filtered) %in% rownames(tax_mat_145_filtered), ]

# Filter sample data to include only ancestral and derived samples
samples_df_filtered <- samples_df |> 
  filter(lineage %in% c("ancestral", "derived"))

# Ensure OTU matrix columns match filtered samples

matched_ids_145 <- intersect(colnames(otu_mat_145_final), samples_df_filtered$sample.id)


otu_mat_145_final <- otu_mat_145_final[, matched_ids_145]

# Reorder samples_df to match OTU matrix column order


samples_df_filtered <- samples_df_filtered |> 
  filter(sample.id %in% matched_ids_145)  |> 
  arrange(match(sample.id, matched_ids_145)) 



# Relative Abundance Calculation ----------------------------------------

# Relative abundance for filtering

otu_mat_145_relative <- sweep(otu_mat_145_final, 2, colSums(otu_mat_145_final), `/`) * 100



keep_taxa_145 <- rownames(otu_mat_145_relative)[apply(otu_mat_145_relative, 1, function(x) any(x >= 1))]
keep_taxa_145 <- keep_taxa_145[!is.na(keep_taxa_145)]

otu_mat_145_filtered_counts <- otu_mat_145_final[keep_taxa_145, ]
tax_mat_145_filtered_counts <- tax_mat_145_filtered[keep_taxa_145, ]


# Filter taxonomy again based on newly calculated rel abundance in the OTU matrix
tax_mat_filtered_counts <- tax_mat_145_filtered[rownames(tax_mat_145_filtered) %in% rownames(otu_mat_145_filtered_counts), ]
tax_mat_145_filtered_counts <- tax_mat_145_filtered[rownames(tax_mat_145_filtered) %in% rownames(otu_mat_145_filtered_counts), ]

# Phyloseq Object Creation --------------------------------------------------
rownames(samples_df_filtered) <- make.unique(as.character(samples_df_filtered$sample.id))


# Step 12: Create phyloseq object
rownames(samples_df_filtered) <- samples_df_filtered$sample.id
samples_df_filtered$sample.id <- NULL  # optional: drop column after it's used as rowname
samples <- sample_data(samples_df_filtered)



OTU <- otu_table(otu_mat_145_filtered_counts, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat_145_filtered_counts)
samples <- sample_data(samples_df_filtered)
physeq <- phyloseq(OTU, TAX, samples)



#creation of phyloseq object
# First, keep only rows of enddat that match columns in the OTU matrix
common_ids <- intersect(enddat$sample.id, colnames(otu_mat_145_filtered_counts))

# Filter both metadata and OTU matrix
enddat_filtered <- enddat[enddat$sample.id %in% common_ids, ]
otu_mat_145_matched <- otu_mat_145_filtered_counts[, common_ids]

# Reorder enddat to match OTU matrix column order
enddat_filtered <- enddat_filtered[match(common_ids, enddat_filtered$sample.id), ]

rownames(enddat_filtered) <- enddat_filtered$sample.id
enddat_filtered$sample.id <- NULL  # optional to remove column

samples <- sample_data(enddat_filtered)


OTU <- otu_table(otu_mat_145_matched, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat_145_filtered_counts)

physeq_end <- phyloseq(OTU, TAX, samples)





