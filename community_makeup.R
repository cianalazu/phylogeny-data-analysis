# Load the dplyr package
library(dplyr)

# Read in the data
community_makeup <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/community_makeup.csv", stringsAsFactors = FALSE)
exp1_data <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/exp1_data.csv", stringsAsFactors = FALSE)
DNA_extractions <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/DNA_extractions.csv", stringsAsFactors = FALSE)
wgs_dilutions <- read.csv("C:/Users/Ciana/OneDrive - USNH/Desktop/mac-files/R/exp_1/wgs_dilutions.csv", stringsAsFactors = FALSE)

#combine the 


# Set encoding to UTF-8 and remove case sensitivity by converting strain_ID to lowercase
community_makeup <- community_makeup %>%
  mutate(strain_ID = iconv(strain_ID, to = "UTF-8") %>% tolower())

DNA_extractions <- DNA_extractions %>%
  mutate(strain_ID = iconv(strain_ID, to = "UTF-8") %>% tolower())

# Join the two data frames on "strain_ID"
merged_df <- community_makeup %>%
  left_join(DNA_extractions %>%
              select(strain_ID, DNA_concentration, DNA_concentration.1), by = "strain_ID")

# View the merged data
print(merged_df)

# Keep only unique strain values
unique_strains_df <- merged_df[!duplicated(merged_df$strain), ]

#in unique_strains_df, what is the prefix of each character in strain_id



# Summarize data by strain_ID (after standardizing)
result_df <- merged_df %>%
  group_by(strain_ID) %>%
  summarise(
    communities = {
      vals <- unique(na.omit(syn_com))
      vals <- vals[vals != ""]
      if (length(vals) == 0) "" else paste(vals, collapse = ", ")
    },
    DNA_concentration = {
      vals <- na.omit(DNA_concentration)
      vals <- vals[vals != ""]
      if (length(vals) == 0) "" else paste(vals, collapse = ", ")
    },
    DNA_concentration_1 = {
      vals <- na.omit(DNA_concentration.1)
      vals <- vals[vals != ""]
      if (length(vals) == 0) "" else paste(vals, collapse = ", ")
    }
  )


# View the result
print(n=137,result_df)

write.csv(result_df, "strains_in_syncoms.csv", row.names = FALSE)

# Filter result_df to get strains with NA or blank values for DNA_concentration
na_or_blank_dna_df <- result_df %>%
  filter(is.na(DNA_concentration) | DNA_concentration == "")

# View the new dataframe with NA or blank values for DNA_concentration
print(na_or_blank_dna_df)


library(dplyr)

library(dplyr)

# Step 1: Remove duplicates in result_df
result_df_clean <- result_df %>%
  distinct(strain_ID, DNA_concentration, .keep_all = TRUE)

# Step 2: Merge wgs_dilutions with DNA_extractions
merged_data <- wgs_dilutions %>%
  left_join(DNA_extractions, by = c("sampleID" = "strain_ID"))

# Step 3: Merge with cleaned result_df to get DNA_concentration
final_data <- merged_data %>%
  left_join(result_df_clean %>% select(strain_ID, DNA_concentration), 
            by = c("sampleID" = "strain_ID"))


# View result
head(final_data)

write.csv(final_data, "dilutions.csv", row.names = FALSE)
