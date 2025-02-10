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
library(ggpubr)  # For better publication-style themes

# Ensure 'genes' is character in all datasets
lower_grade_gene_frequencies <- lower_grade_gene_frequencies %>%
  mutate(genes = as.character(genes), Category = "Lower Grade")

higher_grade_gene_frequencies <- higher_grade_gene_frequencies %>%
  mutate(genes = as.character(genes), Category = "Higher Grade")

metastasis_gene_frequencies <- metastasis_gene_frequencies %>%
  mutate(genes = as.character(genes), Category = "Metastasis")

# Combine all datasets
combined_gene_frequencies <- bind_rows(lower_grade_gene_frequencies,
                                       higher_grade_gene_frequencies,
                                       metastasis_gene_frequencies)

# Convert Category to a factor for proper ordering in facets
combined_gene_frequencies$Category <- factor(combined_gene_frequencies$Category, 
                                             levels = c("Lower Grade", "Higher Grade", "Metastasis"))

# Set a professional color palette (colorblind-friendly)
color_palette <- c("#1b9e77", "#d95f02", "#7570b3")  # Distinct & visually appealing

# Create the three-panel faceted bar plot with professional styling
ggplot(combined_gene_frequencies, aes(x = genes, y = Frequency_Percentage, fill = Category)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # Thin black border for contrast
  facet_grid(Category ~ ., scales = "fixed") +  # Keeps the same y-axis scale for all panels
  labs(title = "Mutation Frequency Across Tumor Grades",
       x = "Genes",
       y = "Mutation Frequency (%)") +  # Shared y-axis label
  scale_fill_manual(values = color_palette) +  # High-contrast, professional color scheme
  theme_pubr(base_size = 16) +  # Professional journal-quality theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"),  # Readable x-axis
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 18, face = "bold"),  # Panel labels in bold
    panel.grid.major = element_blank(),  # No distracting grid lines
    panel.grid.minor = element_blank(),
    legend.position = "none"  # No legend needed since we have facets
  )

# Save the high-resolution figure (600 DPI for publication)
ggsave("Gene_Mutation_Frequency_High_Impact.png", width = 10, height = 12, dpi = 600)
ggsave("Gene_Mutation_Frequency_High_Impact.tiff", width = 10, height = 12, dpi = 600, compression = "lzw")



