library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg19)
library(MutationalPatterns)
library(RColorBrewer)
library(ggrepel)


setwd("C:/Moein/projects/prostate_stages/PRAD_files")

###Preparing data

gleason <- read.delim("C:/Moein/projects/prostate_stages/PRAD_files/gleason.txt")

MAF1 <- preload_maf(maf = "armenia_final.maf.txt", refset = ces.refset.hg19)
MAF4 <- preload_maf(maf = "SU2C_final.maf.txt", refset = ces.refset.hg19)


# removing samples where column Problem is equal to NA
MAF1 <- MAF1[is.na(problem)]
MAF4 <- MAF4[is.na(problem)]


# keeping only samples that do not occur at germline variant sites
MAF1 <- MAF1[germline_variant_site == F]
MAF4 <- MAF4[germline_variant_site == F]


# keeping only samples that do not occur in repetitive regions 
MAF1 <- MAF1[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF4 <- MAF4[(repetitive_region == F | cosmic_site_tier %in% 1:3)]



##Create CESAnalysis and load data:
cesa <- CESAnalysis(refset = "ces.refset.hg19")
cesa <- load_maf(cesa = cesa, maf = MAF1, maf_name = "armenia")
cesa <- load_maf(cesa = cesa, maf = MAF4, maf_name = "SU2C",, coverage = "exome",
                 covered_regions = "SureSelect_All_Exon_covered_regions.bed",
                 covered_regions_name = "SureSelect_V4", covered_regions_padding = 100)


cesa <- load_sample_data(cesa, gleason)

# Load dplyr package
library(dplyr)

# Count unique values in "Unique_Patient_Identifier"


prevalence <- cesa@maf

prevalence <- prevalence %>%
  inner_join(gleason, by = "Unique_Patient_Identifier")


lower_grade <- prevalence %>% filter(Gleason == "Early")
higher_grade <- prevalence %>% filter(Gleason == "Late")
metastasis <- prevalence %>% filter(Gleason == "Metastasis") 

# Count unique patients
unique_patient_count_lower <- lower_grade %>%
  summarise(Unique_Patient_Count = n_distinct(Unique_Patient_Identifier))

unique_patient_count_higher <- higher_grade %>%
  summarise(Unique_Patient_Count = n_distinct(Unique_Patient_Identifier))

unique_patient_count_metastasis <- metastasis %>%
  summarise(Unique_Patient_Count = n_distinct(Unique_Patient_Identifier))

# Define genes of interest
genes_of_interest <- c("SPOP", "FOXA1", "AR", "PIK3CA", "PIK3CB", "TP53", "ROCK1", "RHOA", "AKT1", "ATM", "CUL3",
                       "APC", "CTNNB1", "MUC16", "KMT2C", "KMT2D")

# Filter data for only the genes of interest

filtered_lower_grade <- lower_grade %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)  # Remove duplicate patient-gene entries

filtered_higher_grade <- higher_grade %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)  

filtered_metastasis <- metastasis %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)  

# Count unique patients per gene

lower_grade_gene_frequencies <- filtered_lower_grade %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n())

higher_grade_gene_frequencies <- filtered_higher_grade %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n())

metastasis_gene_frequencies <- filtered_metastasis %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n())

# Count unique patients per gene
lower_grade_gene_frequencies <- filtered_lower_grade %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n(), .groups = "drop") %>%
  mutate(Frequency_Percentage = as.numeric((Unique_Patient_Count / 293) * 100))  # Normalize by total patients

# Ensure all genes are represented, even if missing
lower_grade_gene_frequencies <- tibble(genes = genes_of_interest) %>%
  left_join(lower_grade_gene_frequencies, by = "genes") %>%
  mutate(Unique_Patient_Count = replace_na(Unique_Patient_Count, 0),
         Frequency_Percentage = replace_na(Frequency_Percentage, 0))

higher_grade_gene_frequencies <- filtered_higher_grade %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n()) %>%
  mutate(Frequency_Percentage = as.numeric((Unique_Patient_Count / 320) * 100))  # Calculate frequency

metastasis_gene_frequencies <- filtered_metastasis %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n()) %>%
  mutate(Frequency_Percentage = as.numeric((Unique_Patient_Count / 480) * 100))  # Calculate frequency

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Combine all datasets into one for grouped plotting
lower_grade_gene_frequencies$Category <- "Lower Grade"
higher_grade_gene_frequencies$Category <- "Higher Grade"
metastasis_gene_frequencies$Category <- "Metastasis"

# Ensure 'genes' column is character in all datasets
lower_grade_gene_frequencies <- lower_grade_gene_frequencies %>%
  mutate(genes = as.character(genes))

higher_grade_gene_frequencies <- higher_grade_gene_frequencies %>%
  mutate(genes = as.character(genes))

metastasis_gene_frequencies <- metastasis_gene_frequencies %>%
  mutate(genes = as.character(genes))

# Now combine datasets
combined_gene_frequencies <- bind_rows(lower_grade_gene_frequencies,
                                       higher_grade_gene_frequencies,
                                       metastasis_gene_frequencies)

# Convert Category to a factor for better ordering
combined_gene_frequencies$Category <- factor(combined_gene_frequencies$Category, 
                                             levels = c("Lower Grade", "Higher Grade", "Metastasis"))

# Plot grouped bar chart
ggplot(combined_gene_frequencies, aes(x = genes, y = Frequency_Percentage, fill = Category)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    labs(title = "Gene Mutation Frequency in Different Grades",
         x = "Genes",
         y = "Frequency Percentage (%)") +
    scale_fill_manual(values = c("#8dd3c7", "#ffffb3", "#bebada")) +  # Lighter colors for better visibility
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


