library(tidyr)
library(dplyr)
library(decontam)
library(tibble)

# ───────────────────────────────────────────────────────────────
# 1. Setup
# ───────────────────────────────────────────────────────────────
setwd("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/16s_experimental_evolution")

# ───────────────────────────────────────────────────────────────
# 2. Load Data
# ───────────────────────────────────────────────────────────────
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
feature_table <- read.table("feature-table.tsv", header = TRUE, check.names = FALSE, sep = "\t")
taxonomy <- read.table("taxonomy.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
chordata_taxa <- read.csv("ncbi_taxid_Chordata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
streptophyta_taxa <- read.csv("ncbi_taxid_Streptophyta.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
family_taxa_df <- read.csv("family_taxa.csv", header = TRUE, stringsAsFactors = FALSE)

# ───────────────────────────────────────────────────────────────
# 3. Taxonomy Cleanup
# ───────────────────────────────────────────────────────────────
taxonomy <- taxonomy |> 
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"), sep = ";", remove = TRUE) |> 
  mutate(across(Kingdom:Genus, ~sub("^.*?__", "", .))) |> 
  select(-Confidence)

# ───────────────────────────────────────────────────────────────
# 4. Identify & Remove Blanks 
# ───────────────────────────────────────────────────────────────

otu_table <- feature_table |>
  column_to_rownames("OTU ID") |>  # OTU ID becomes row names
  as.matrix()

# Convert to numeric matrix
otu_table <- apply(otu_table, 2, as.numeric)
rownames(otu_table) <- feature_table$`OTU ID`

# Calculate relative abundance
relabund <- sweep(otu_table, 2, colSums(otu_table), `/`)

# Subset blanks
relabund_blanks <- relabund[, grepl("blank", colnames(relabund), ignore.case = TRUE)]

# Identify features with >5% abundance in blanks
poss_blank <- rownames(relabund)[
  rowSums(relabund[, c("blank-1", "blank-2"), drop = FALSE]) > 0.05
]

# Remove blank contaminants
feature_no_blanks <- feature_table |> 
  filter(!`OTU ID` %in% poss_blank)

taxonomy_no_blanks <- taxonomy |> 
  filter(!Feature.ID %in% poss_blank)

# ───────────────────────────────────────────────────────────────
# 5. Remove Host Contaminants (Chordata, Streptophyta, mito/chloro/unassigned)
# ───────────────────────────────────────────────────────────────

combined_taxa <- bind_rows(chordata_taxa, streptophyta_taxa)

# Get indices
Ch_St <- which(taxonomy_no_blanks$Family %in% combined_taxa$Family.name)
Chloro <- which(taxonomy_no_blanks$Family == "Chloroplast")
Mitos <- which(taxonomy_no_blanks$Family == "Mitochondria")
Unassigned <- which(taxonomy_no_blanks$Kingdom == "Unassigned")

# Remove those features
taxfilt <- taxonomy_no_blanks[-c(Ch_St, Chloro, Mitos, Unassigned), ]

featfilt <- feature_no_blanks |>
  filter(!`OTU ID` %in% taxonomy_no_blanks$Feature.ID[c(Ch_St, Chloro, Mitos, Unassigned)])

# ───────────────────────────────────────────────────────────────
# 6. Save Cleaned Files
# ───────────────────────────────────────────────────────────────
colnames(featfilt)[1] <- "OTU_ID"
write.table(featfilt, "feature_table_no_contam.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
write.table(taxfilt, "taxonomy_no_contam.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# ───────────────────────────────────────────────────────────────
# 7. Report Summary
# ───────────────────────────────────────────────────────────────


cat("Number of OTUs before filtering:", nrow(feature_table), "\n")
#how many removed just from blanks filtering
cat("Number of taxa removed from blanks filtering:", length(poss_blank), "\n")
#how many removed from host contaminants
cat("Number of taxa removed from host contaminants:", 
    length(Ch_St) + length(Chloro) + length(Mitos) + length(Unassigned), "\n")
cat("Number of OTUs after filtering:", nrow(featfilt), "\n")
