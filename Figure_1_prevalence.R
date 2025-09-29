library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg19)
library(MutationalPatterns)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)  

###Preparing data

gleason <- read.delim("C:/Moein/projects/prostate_stages/PRAD_files/gleason.txt")

MAF1 <- preload_maf(maf = "armenia_final.maf.txt", refset = ces.refset.hg19)
MAF2 <- preload_maf(maf = "boutros_final.maf.txt", refset = ces.refset.hg19)
MAF3 <- preload_maf(maf = "tcga_final.maf.txt", refset = ces.refset.hg19)
MAF4 <- preload_maf(maf = "SU2C_final.maf.txt", refset = ces.refset.hg19)
MAF5 <- preload_maf(maf = "MSK_341_final.maf.txt", refset = ces.refset.hg19)
MAF6 <- preload_maf(maf = "MSK_410_final.maf.txt", refset = ces.refset.hg19)
MAF7 <- preload_maf(maf = "MSK_468_final.maf.txt", refset = ces.refset.hg19)

# removing samples where column Problem is equal to NA
MAF1 <- MAF1[is.na(problem)]
MAF2 <- MAF2[is.na(problem)]
MAF3 <- MAF3[is.na(problem)]
MAF4 <- MAF4[is.na(problem)]
MAF5 <- MAF5[is.na(problem)]
MAF6 <- MAF6[is.na(problem)]
MAF7 <- MAF7[is.na(problem)]

# keeping only samples that do not occur at germline variant sites
MAF1 <- MAF1[germline_variant_site == F]
MAF2 <- MAF2[germline_variant_site == F]
MAF3 <- MAF3[germline_variant_site == F]
MAF4 <- MAF4[germline_variant_site == F]
MAF5 <- MAF5[germline_variant_site == F]
MAF6 <- MAF6[germline_variant_site == F]
MAF7 <- MAF7[germline_variant_site == F]

# keeping only samples that do not occur in repetitive regions 
MAF1 <- MAF1[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF2 <- MAF2[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF3 <- MAF3[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF4 <- MAF4[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF5 <- MAF5[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF6 <- MAF6[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF7 <- MAF7[(repetitive_region == F | cosmic_site_tier %in% 1:3)]

# keeping snv:
MAF1 <- subset(MAF1, variant_type == "snv")
MAF2 <- subset(MAF2, variant_type == "snv")
MAF3 <- subset(MAF3, variant_type == "snv")
MAF4 <- subset(MAF4, variant_type == "snv")
MAF5 <- subset(MAF5, variant_type == "snv")
MAF6 <- subset(MAF6, variant_type == "snv")
MAF7 <- subset(MAF7, variant_type == "snv")

##Create CESAnalysis and load data:
cesa <- CESAnalysis(refset = "ces.refset.hg19")
cesa <- load_maf(cesa = cesa, maf = MAF1, maf_name = "armenia")
cesa <- load_maf(cesa = cesa, maf = MAF2, maf_name = "boutros", coverage = "genome")
cesa <- load_maf(cesa = cesa, maf = MAF3, maf_name = "tcga", coverage = "genome")
cesa <- load_maf(cesa = cesa, maf = MAF4, maf_name = "SU2C",, coverage = "exome",
                 covered_regions = "SureSelect_All_Exon_covered_regions.bed",
                 covered_regions_name = "SureSelect_V4", covered_regions_padding = 100)
cesa <- load_maf(cesa = cesa, maf = MAF5, maf_name = "341", coverage = "targeted", 
                 covered_regions = "msk_341_exons.bed", 
                 covered_regions_name = "MSK_IMPACT_341", covered_regions_padding = 10)
cesa <- load_maf(cesa = cesa, maf = MAF6, maf_name = "410", coverage = "targeted", 
                 covered_regions = "msk_410_exons.bed", 
                 covered_regions_name = "MSK_IMPACT_410", covered_regions_padding = 10)
cesa <- load_maf(cesa = cesa, maf = MAF7, maf_name = "468", coverage = "targeted", 
                 covered_regions = "msk_468_exons.bed", 
                 covered_regions_name = "MSK_IMPACT_468", covered_regions_padding = 10)


cesa <- load_sample_data(cesa, gleason)

# Extract prevalence data
prevalence <- cesa@maf %>%
  inner_join(gleason, by = "Unique_Patient_Identifier")

# Classify data into different Gleason grades
low_grade <- prevalence %>% filter(Gleason == "Early")
high_grade <- prevalence %>% filter(Gleason == "Late")
mCRPC <- prevalence %>% filter(Gleason == "Metastasis")

#Count the number of samples in each group:
low_grade_samples <- low_grade %>%
  summarise(Unique_Count = n_distinct(Unique_Patient_Identifier))

high_grade_samples <- high_grade %>%
  summarise(Unique_Count = n_distinct(Unique_Patient_Identifier))

mCRPC_samples <- mCRPC %>%
  summarise(Unique_Count = n_distinct(Unique_Patient_Identifier))

# Define genes of interest with correct order (Removed "MUC16")
genes_of_interest <- c("SPOP", "TP53", "KMT2C", "FOXA1", "PIK3CA", "ATM", "CTNNB1",
                       "KMT2D", "PTEN", "CUL3", "AR", "AKT1", "ROCK1", "PIK3CB", "APC", "RHOA")

# Filter data for genes of interest
filtered_low_grade <- low_grade %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)

filtered_high_grade <- high_grade %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)

filtered_mCRPC <- mCRPC %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)

# Count unique patients per gene and normalize frequency
low_grade_gene_frequencies <- filtered_low_grade %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n(), .groups = "drop") %>%
  mutate(Frequency_Percentage = (Unique_Patient_Count / 479) * 100, Category = "Low-Grade")

high_grade_gene_frequencies <- filtered_high_grade %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n(), .groups = "drop") %>%
  mutate(Frequency_Percentage = (Unique_Patient_Count / 406) * 100, Category = "High-Grade")

mCRPC_gene_frequencies <- filtered_mCRPC %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n(), .groups = "drop") %>%
  mutate(Frequency_Percentage = (Unique_Patient_Count / 1024) * 100, Category = "mCRPC")

# Combine datasets
combined_gene_frequencies <- bind_rows(low_grade_gene_frequencies,
                                       high_grade_gene_frequencies,
                                       mCRPC_gene_frequencies)

# Convert Category to a factor for correct grouping
combined_gene_frequencies$Category <- factor(combined_gene_frequencies$Category,
                                             levels = c("Low-Grade", "High-Grade", "mCRPC"))

# Convert genes to a factor with the correct order 
combined_gene_frequencies$genes <- factor(combined_gene_frequencies$genes, levels = genes_of_interest)

# Define new distinct color scheme for each tumor stage
stage_colors <- c("Low-Grade" = "#1b9e77",  # Green
                  "High-Grade" = "#d95f02",  # Orange
                  "mCRPC" = "#7570b3")   # Blue

# Generate one-panel grouped bar chart
ggplot(combined_gene_frequencies, aes(x = genes, y = Frequency_Percentage, fill = Category)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, position = position_dodge(width = 0.8)) +  
  labs(x = "Genes", y = "Mutation Frequency (%)") +  
  scale_fill_manual(values = stage_colors) +  
  theme_pubr(base_size = 16) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.position = "top",  
    legend.title = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Figure_1.png", width = 12, height = 6, dpi = 600)

#End


