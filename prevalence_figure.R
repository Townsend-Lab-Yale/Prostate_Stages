# Load necessary libraries
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

# Set working directory
setwd("C:/Moein/projects/prostate_stages/PRAD_files")

### Preparing data
gleason <- read.delim("C:/Moein/projects/prostate_stages/PRAD_files/gleason.txt")

# Load MAF files
MAF1 <- preload_maf(maf = "armenia_final.maf.txt", refset = ces.refset.hg19)
MAF4 <- preload_maf(maf = "SU2C_final.maf.txt", refset = ces.refset.hg19)

# Removing samples where column Problem is NA
MAF1 <- MAF1[is.na(problem)]
MAF4 <- MAF4[is.na(problem)]

# Keeping only samples that do not occur at germline variant sites
MAF1 <- MAF1[germline_variant_site == FALSE]
MAF4 <- MAF4[germline_variant_site == FALSE]

# Keeping only samples that do not occur in repetitive regions
MAF1 <- MAF1[(repetitive_region == FALSE | cosmic_site_tier %in% 1:3)]
MAF4 <- MAF4[(repetitive_region == FALSE | cosmic_site_tier %in% 1:3)]

# Create CESAnalysis and load data
cesa <- CESAnalysis(refset = "ces.refset.hg19")
cesa <- load_maf(cesa = cesa, maf = MAF1, maf_name = "armenia")
cesa <- load_maf(cesa = cesa, maf = MAF4, maf_name = "SU2C", coverage = "exome",
                 covered_regions = "SureSelect_All_Exon_covered_regions.bed",
                 covered_regions_name = "SureSelect_V4", covered_regions_padding = 100)

cesa <- load_sample_data(cesa, gleason)

# Extract prevalence data
prevalence <- cesa@maf %>%
  inner_join(gleason, by = "Unique_Patient_Identifier")

# Classify data into different Gleason grades
lower_grade <- prevalence %>% filter(Gleason == "Early")
higher_grade <- prevalence %>% filter(Gleason == "Late")
metastasis <- prevalence %>% filter(Gleason == "Metastasis")

# Define genes of interest with correct order
genes_of_interest <- c("SPOP", "TP53", "MUC16", "KMT2C", "FOXA1", "PIK3CA", "ATM", "CTNNB1",
                       "KMT2D", "CUL3", "AR", "AKT1", "ROCK1", "PIK3CB", "APC", "RHOA")

# Filter data for genes of interest
filtered_lower_grade <- lower_grade %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)

filtered_higher_grade <- higher_grade %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)

filtered_metastasis <- metastasis %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)

# Count unique patients per gene
lower_grade_gene_frequencies <- filtered_lower_grade %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n(), .groups = "drop") %>%
  mutate(Frequency_Percentage = (Unique_Patient_Count / 293) * 100)

higher_grade_gene_frequencies <- filtered_higher_grade %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n(), .groups = "drop") %>%
  mutate(Frequency_Percentage = (Unique_Patient_Count / 320) * 100)

metastasis_gene_frequencies <- filtered_metastasis %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n(), .groups = "drop") %>%
  mutate(Frequency_Percentage = (Unique_Patient_Count / 480) * 100)

# Add Category column
lower_grade_gene_frequencies <- lower_grade_gene_frequencies %>%
  mutate(Category = "Lower Grade")

higher_grade_gene_frequencies <- higher_grade_gene_frequencies %>%
  mutate(Category = "Higher Grade")

metastasis_gene_frequencies <- metastasis_gene_frequencies %>%
  mutate(Category = "Metastasis")

# Combine datasets
combined_gene_frequencies <- bind_rows(lower_grade_gene_frequencies,
                                       higher_grade_gene_frequencies,
                                       metastasis_gene_frequencies)

# Convert Category to a factor for correct facet order
combined_gene_frequencies$Category <- factor(combined_gene_frequencies$Category,
                                             levels = c("Lower Grade", "Higher Grade", "Metastasis"))

# Convert genes to a factor with the correct order
combined_gene_frequencies$genes <- factor(combined_gene_frequencies$genes, levels = genes_of_interest)

# Define highlighted genes (darker colors)
highlighted_genes <- list(
  "Lower Grade" = c("ATM", "SPOP", "PIK3CA", "FOXA1"),
  "Higher Grade" = c("SPOP", "PIK3CA", "ATM", "KMT2D", "FOXA1", "APC"),
  "Metastasis" = c("AR")
)

# Define color palettes (light and dark)
color_schemes <- list(
  "Lower Grade" = c("lighter" = "#a7d7c5", "darker" = "#11634d"),  # Green shades
  "Higher Grade" = c("lighter" = "#ffcf99", "darker" = "#a53f00"),  # Orange shades
  "Metastasis" = c("lighter" = "#b3a3d6", "darker" = "#4e3d7a")   # Blue shades
)

# Assign colors to genes based on category
combined_gene_frequencies <- combined_gene_frequencies %>%
  rowwise() %>%
  mutate(
    Fill_Color = ifelse(
      genes %in% highlighted_genes[[Category]],  # Check if gene is highlighted
      color_schemes[[Category]]["darker"],  # Assign darker color if highlighted
      color_schemes[[Category]]["lighter"]   # Assign lighter color otherwise
    )
  ) %>%
  ungroup()

# Generate three-panel faceted bar plot (without title)
ggplot(combined_gene_frequencies, aes(x = genes, y = Frequency_Percentage, fill = Fill_Color)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # Thin black border for contrast
  facet_grid(Category ~ ., scales = "fixed") +  # Ensures consistent Y-axis scaling
  labs(x = "Genes", y = "Mutation Frequency (%)") +  # Only X and Y labels, no title
  scale_fill_identity() +  # Apply precomputed colors
  theme_pubr(base_size = 16) +  # Professional journal-quality theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    strip.text = element_text(size = 18, face = "bold"),  # Bold panel labels
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"  
  )

ggsave("Gene_Mutation_Frequency_High_Impact.png", width = 10, height = 12, dpi = 600)
ggsave("Gene_Mutation_Frequency_High_Impact.tiff", width = 10, height = 12, dpi = 600, compression = "lzw")
ggsave("Gene_Mutation_Frequency_High_Impact.pdf", width = 10, height = 12)
