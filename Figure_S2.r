# Load libraries
library(cancereffectsizeR)
library(ces.refset.hg19)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)

setwd("C:/Moein/projects/prostate_stages/PRAD_files")

### Preparing data
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

# Create CESAnalysis and load data:
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

# Mutational processes and relative mutation rates:
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

cesa = trinuc_mutation_rates(cesa, ces.refset.hg19$signatures$COSMIC_v3.2,
                                  signature_exclusions = signature_exclusions)

# Extract biological weights
bio_weights <- cesa$mutational_signatures$biological_weights
sbs_cols <- grep("^SBS", names(bio_weights), value = TRUE)
sbs_data <- bio_weights[, ..sbs_cols]
nonzero_cols <- sbs_cols[colSums(sbs_data != 0, na.rm = TRUE) > 0]
final_weights <- cbind(bio_weights[, .(Unique_Patient_Identifier)], bio_weights[, ..nonzero_cols])

# Merge with Gleason data
merged_data <- final_weights %>%
  inner_join(gleason, by = "Unique_Patient_Identifier")

# Reshape to long format
long_data <- melt(merged_data, id.vars = c("Unique_Patient_Identifier", "Gleason"),
                  variable.name = "Signature", value.name = "Weight")

# Assign biological signature groups
signature_groupings <- list(
  "Deamination with age, clock-like (1)" = "SBS1",
  "Unknown, clock-like (5)" = "SBS5",
  "APOBEC (2, 13)" = c("SBS2", "SBS13"),
  "Defective homologous recombination (3)" = "SBS3",
  "Tobacco (4, 29)" = c("SBS4", "SBS29"),
  "UV light (7a–d, 38)" = c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38"),
  "Prior treatment (11, 31, 32, 35)" = c("SBS11", "SBS31", "SBS32", "SBS35"),
  "Mutagenic chemical exposure (22, 24, 42, 88)" = c("SBS22", "SBS24", "SBS42", "SBS88"),
  "Alcohol-associated (16)" = "SBS16"
)

signature_to_group <- data.frame(
  Signature = unlist(signature_groupings),
  Group = rep(names(signature_groupings), lengths(signature_groupings)),
  stringsAsFactors = FALSE
)

long_data <- long_data %>%
  left_join(signature_to_group, by = "Signature") %>%
  mutate(Group = ifelse(is.na(Group), "Non-actionable and unknown signatures", Group))

# Filter and recode
signature_order <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS6", "SBS8", "SBS10c", "SBS10d",
                     "SBS12", "SBS13", "SBS18", "SBS33", "SBS37", "SBS39", "SBS40",
                     "SBS41", "SBS86", "SBS87")

long_data <- long_data %>%
  filter(Signature %in% signature_order) %>%
  mutate(Signature = factor(gsub("^SBS", "", Signature), levels = gsub("^SBS", "", signature_order)),
         Gleason = recode(Gleason, "Early" = "Low-grade", "Late" = "High-grade", "Metastasis" = "mCRPC"),
         Gleason = factor(Gleason, levels = c("Low-grade", "High-grade", "mCRPC")))

# Plot settings
cannataro_colors <- c(
  "Deamination with age, clock-like (1)" = "gray40",
  "Unknown, clock-like (5)" = "gray60",
  "APOBEC (2, 13)" = "#7570b3",
  "Defective homologous recombination (3)" = "#e7298a",
  "Tobacco (4, 29)" = "#a6761d",
  "UV light (7a–d, 38)" = "#e6ab02",
  "Prior treatment (11, 31, 32, 35)" = "#1b9e77",
  "Mutagenic chemical exposure (22, 24, 42, 88)" = "#66a61e",
  "Alcohol-associated (16)" = "#d95f02",
  "Non-actionable and unknown signatures" = "black"
)

clean_theme <- theme_classic(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.line = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

y_scale <- scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))

# Violin plot function
make_violin <- function(gleason_label, title) {
  ggplot(filter(long_data, Gleason == gleason_label),
         aes(x = Signature, y = Weight, fill = Group)) +
    geom_violin(scale = "width", color = "black", adjust = 1) +
    scale_fill_manual(values = cannataro_colors) +
    y_scale +
    labs(x = "COSMIC signature", y = "Signature weight") +
    ggtitle(title) +
    clean_theme
}

# Generate and combine plots
p_low <- make_violin("Low-grade", "Low-grade")
p_high <- make_violin("High-grade", "High-grade")
p_mcrpc <- make_violin("mCRPC", "mCRPC")

combined_plot <- plot_grid(p_low, p_high, p_mcrpc,
                           labels = c("A", "B", "C"),
                           ncol = 1, label_fontface = "bold")

# Save plot
ggsave("Supplementary_Figure_S2.png", combined_plot, width = 6, height = 9, dpi = 300, bg = "white")
