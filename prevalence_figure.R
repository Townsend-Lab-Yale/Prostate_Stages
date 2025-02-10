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
  summarise(Unique_Patient_Count = n()) %>%
  mutate(Frequency_Percentage = as.numeric((Unique_Patient_Count / 293) * 100))  # Calculate frequency

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
library(RColorBrewer)


# Save as high-resolution PNG with lighter colors
ggsave("gene_prevalence_plot_lighter.png", width = 10, height = 6, dpi = 300)

# Save as high-resolution TIFF
ggsave("gene_prevalence_plot_lighter.tiff", width = 10, height = 6, dpi = 600, compression = "lzw")



