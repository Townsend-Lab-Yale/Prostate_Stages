library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg19)
library(MutationalPatterns)
library(RColorBrewer)
library(ggrepel)
library(readr)
library(scales)
library(stringr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)

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

###creating CESAnalysis and loading data 

cesa <- CESAnalysis(refset = "ces.refset.hg19")

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

#defining groups:
Early_groups_all <- cesa$samples[Gleason == "Early", unique(Unique_Patient_Identifier)]
Late_groups_all <- cesa$samples[Gleason == "Late", unique(Unique_Patient_Identifier)]
Metastasis_groups_all <- cesa$samples[Gleason == "Metastasis", unique(Unique_Patient_Identifier)]

#defining groups for gene mutation rate using exome:
Early_groups <- cesa$samples[Gleason == "Early" & coverage == "exome", unique(Unique_Patient_Identifier)]
Late_groups <- cesa$samples[Gleason == "Late" & coverage == "exome", unique(Unique_Patient_Identifier)]
Metastasis_groups <- cesa$samples[Gleason == "Metastasis" & coverage == "exome", unique(Unique_Patient_Identifier)]

cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa, covariates = "PRAD", samples = Early_groups, save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "PRAD", samples = Late_groups, save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "PRAD", samples = Metastasis_groups, save_all_dndscv_output = T)


selected_genes <- c("SPOP", "FOXA1", "AR", "PIK3CA", "PIK3CB", "TP53", "ROCK1", "RHOA", "AKT1", "ATM", "CUL3",
                    "APC", "CTNNB1", "PTEN", "KMT2C", "KMT2D")

RefCDS = ces.refset.hg19$RefCDS
dndscv_gene_names <- cesa_samples_by_groups$gene_rates$gene
nsyn_sites = sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

# selecting mutation rate data for samples in Early_groups 
samples_in_Early_groups <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_1$annotmuts$sampleID ))

# selecting mutation rate data for samples in Late_groups
samples_in_Late_groups <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_2$annotmuts$sampleID ))

# selecting mutation rate data for samples in Metastasis_groups
samples_in_Metastasis_groups <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_3$annotmuts$sampleID ))

### creating a data frame with mutation rate data for Early_group, Late_groups and Metastasis_groups
mut_rate_df <- tibble(gene = cesa_samples_by_groups$dNdScv_results$rate_grp_2$genemuts$gene_name,
                      exp_Early_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_Late_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv,
                      exp_Metastasis_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_3$genemuts$exp_syn_cv)

# Add n_syn_sites column to mut_rate_df:
mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df %>% 
  mutate(Early_mu = (exp_Early_mu / n_syn_sites) / samples_in_Early_groups) %>%
  mutate(Late_mu = (exp_Late_mu / n_syn_sites) / samples_in_Late_groups) %>%
  mutate(Metastasis_mu = (exp_Metastasis_mu / n_syn_sites) / samples_in_Metastasis_groups) %>%
  mutate(cancer_greater = (Metastasis_mu > Late_mu) & (Metastasis_mu > Early_mu) & (Late_mu > Early_mu)) -> 
  mut_rate_df

# change in mutation rate across stages
mut_rate_df <- mut_rate_df %>% 
  select(gene, Early_mu, Late_mu, Metastasis_mu) %>%
  mutate(p_1 = Early_mu / Metastasis_mu) %>%
  mutate(p_2 = (Late_mu - Early_mu)/Metastasis_mu) %>%
  mutate(p_3 = 1 - (p_1+p_2))

# saving gene mutation rates into separate data frame:
Early_rate <- mut_rate_df %>%
  select(gene, rate = Early_mu) %>%
  data.table::setDT()
Late_rate <- mut_rate_df %>%
  select(gene, rate = Late_mu) %>%
  data.table::setDT()
Meta_rate <- mut_rate_df %>%
  select(gene, rate = Metastasis_mu) %>%
  data.table::setDT()

# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# setting gene rates for Early, Late and Metastasis:
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = Meta_rate, missing_genes_take_nearest = T, samples = cesa$samples[Gleason %in% c("Early", "Late", "Metastasis")]) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
cesa_samples_by_groups <- trinuc_mutation_rates(cesa = cesa_samples_by_groups, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)

# defining compound variants
compound <- define_compound_variants(cesa = cesa_samples_by_groups, 
                                     variant_table = cesa_samples_by_groups$variants |>
                                       filter(intergenic == F, gene %in% selected_genes),
                                     by = "gene", merge_distance = Inf)

source("new_sequential_lik.R")

for(comp_ind in 1:length(compound)){
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  these_props <- mut_rate_df[mut_rate_df$gene == this_gene,c("p_1","p_2","p_3")]
  these_props <- c(these_props$p_1, these_props$p_2, these_props$p_3)
  
  cesa_samples_by_groups <- ces_variant(cesa = cesa_samples_by_groups, variants = this_comp, model = sequential_lik_dev, 
                                        ordering_col = 'Gleason', ordering = c('Early', 'Late', 'Metastasis'), 
                                        lik_args = list(sequential_mut_prop = these_props), run_name = this_gene)
  
}

# selecting necessary data
selection_data_Early_Late_Metastasis <- rbindlist(cesa_samples_by_groups$selection)


# custom analysis for AR since the assumption for O->L->H->M doesn't hold from L->H
# instead we just consider organogenesis to primary tumor to metastasis

rm(cesa_samples_by_groups)

pri_met <- gleason %>%
  mutate(Pri_Met = ifelse(Gleason == "Metastasis", "Metastasis", "Primary")) %>%
  select(-Gleason)
cesa <- load_sample_data(cesa, pri_met)

#defining groups:
Primary_groups_all <- cesa$samples[Pri_Met == "Primary", unique(Unique_Patient_Identifier)]
Metastasis_groups_all <- cesa$samples[Pri_Met == "Metastasis", unique(Unique_Patient_Identifier)]

#defining groups for gene mutation rate using exome:
Primary_groups <- cesa$samples[Pri_Met == "Primary" & coverage == "exome", unique(Unique_Patient_Identifier)]
Metastasis_groups <- cesa$samples[Pri_Met == "Metastasis" & coverage == "exome", unique(Unique_Patient_Identifier)]

cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa, covariates = "PRAD", samples = Primary_groups, save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "PRAD", samples = Metastasis_groups, save_all_dndscv_output = T)


RefCDS = ces.refset.hg19$RefCDS
dndscv_gene_names <- cesa_samples_by_groups$gene_rates$gene
nsyn_sites = sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

# selecting mutation rate data for samples in Primary_groups and Metastasis_groups
samples_in_Primary_groups <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_1$annotmuts$sampleID ))
samples_in_Metastasis_groups <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_2$annotmuts$sampleID ))

### creating a data frame with mutation rate data for Primary_groups and Metastasis_groups
mut_rate_df <- tibble(gene = cesa_samples_by_groups$dNdScv_results$rate_grp_2$genemuts$gene_name,
                      exp_Primary_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_Metastasis_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv)

# Add n_syn_sites column to mut_rate_df:
mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df <- mut_rate_df %>% 
  mutate(Primary_mu = (exp_Primary_mu / n_syn_sites) / samples_in_Primary_groups) %>%
  mutate(Metastasis_mu = (exp_Metastasis_mu / n_syn_sites) / samples_in_Metastasis_groups) %>%
  mutate(cancer_greater = (Metastasis_mu > Primary_mu))

# change in mutation rate across stages
mut_rate_df <- mut_rate_df %>% 
  select(gene, Primary_mu, Metastasis_mu) %>%
  mutate(p_1 = Primary_mu / Metastasis_mu) %>%
  mutate(p_2 = 1 - p_1)

# saving gene mutation rates into separate data frame:
Primary_rate <- mut_rate_df %>%
  select(gene, rate = Primary_mu) %>%
  data.table::setDT()
Meta_rate <- mut_rate_df %>%
  select(gene, rate = Metastasis_mu) %>%
  data.table::setDT()

# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# setting gene rates for Primary and Metastasis:
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = Meta_rate, missing_genes_take_nearest = T, samples = cesa$samples[Pri_Met %in% c("Primary", "Metastasis")]) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
cesa_samples_by_groups <- trinuc_mutation_rates(cesa = cesa_samples_by_groups, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)

# defining compound variants
compound <- define_compound_variants(cesa = cesa_samples_by_groups, 
                                     variant_table = cesa_samples_by_groups$variants |>
                                       filter(intergenic == F, gene %in% selected_genes),
                                     by = "gene", merge_distance = Inf)

for(comp_ind in 1:length(compound)){
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  these_props <- mut_rate_df[mut_rate_df$gene == this_gene,c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  
  cesa_samples_by_groups <- ces_variant(cesa = cesa_samples_by_groups, variants = this_comp, model = sequential_lik_dev, 
                                        ordering_col = 'Pri_Met', ordering = c('Primary', 'Metastasis'), 
                                        lik_args = list(sequential_mut_prop = these_props), run_name = this_gene)
  
}

# selecting necessary data
AR_selection_data_Primary_Metastasis <- cesa_samples_by_groups$selection$AR

library(tidyverse)

# reformatting data set
selection_data_Early_Late_Metastasis <- selection_data_Early_Late_Metastasis |> 
  select(variant_name, starts_with("si"), starts_with("ci")) |>
  pivot_longer(cols = -variant_name, names_to = "data_type") |>
  mutate(stage = stringr::word(string = data_type, sep = "_",start = -1)) |>
  mutate(variant_name = stringr::str_remove(variant_name, "\\.1")) |>
  mutate(si_or_ci = stringr::word(string = data_type, sep = "_",start = 1, end=3)) |>
  mutate(si_or_ci = case_when(is.na(si_or_ci) ~ "si", TRUE ~ si_or_ci)) |>
  mutate (value = case_when (is.na(value)~0, TRUE~value))

# pivoting data set to create columns for gene, stage, si, and CIs
selection_data_Early_Late_Metastasis <- selection_data_Early_Late_Metastasis|> 
  select(-data_type) |>
  pivot_wider(values_from = value, names_from = si_or_ci)

# defining stages to be plotted
selection_data_Early_Late_Metastasis$stage <- factor(selection_data_Early_Late_Metastasis$stage, levels = c("Early","Late","Metastasis"))

### Making the Figure:

# Rename a specific word in the Name column
selection_data_Early_Late_Metastasis$stage <- sub("Early", "O → L", selection_data_Early_Late_Metastasis$stage)
selection_data_Early_Late_Metastasis$stage <- sub("Late", "L → H", selection_data_Early_Late_Metastasis$stage)
selection_data_Early_Late_Metastasis$stage <- sub("Metastasis", "H → M", selection_data_Early_Late_Metastasis$stage)

variant_order <- c("SPOP", "AKT1", "KMT2D", "CTNNB1", "CUL3", "PIK3CA", "TP53", "FOXA1", "ATM", "KMT2C", "PTEN", "APC", "RHOA", "ROCK1", "PIK3CB", "AR") 
stage_order <- c("O → L", "L → H", "H → M")
selection_data_Early_Late_Metastasis$stage <- factor(selection_data_Early_Late_Metastasis$stage, levels = stage_order)

distance1 <- 1.1  # Distance between "Lower-grade" and "Higher-grade"
distance2 <- 1.0  # Distance between "Metastasis" and "Higher-grade"


selection_data_Early_Late_Metastasis$stage_adjusted <- ifelse(selection_data_Early_Late_Metastasis$stage == "O → L", 0.1,
                                                              ifelse(selection_data_Early_Late_Metastasis$stage == "L → H", 0.1 + distance1,
                                                                     ifelse(selection_data_Early_Late_Metastasis$stage == "H → M", 0.1 + distance1 + distance2,
                                                                            NA)))										 

AR_selection_data_Primary_Metastasis <- AR_selection_data_Primary_Metastasis %>% 
  select(variant_name, starts_with("si"), starts_with("ci")) %>% 
  pivot_longer(cols = -variant_name, names_to = "data_type") %>% 
  mutate(stage = stringr::word(string = data_type, sep = "_",start = -1)) %>% 
  mutate(variant_name = stringr::str_remove(variant_name, "\\.1")) %>% 
  mutate(si_or_ci = stringr::word(string = data_type, sep = "_",start = 1, end=3)) %>% 
  mutate(si_or_ci = case_when(is.na(si_or_ci) ~ "si", TRUE ~ si_or_ci)) %>% 
  mutate (value = case_when (is.na(value)~0, TRUE~value)) %>%
  select(-data_type) %>% 
  pivot_wider(values_from = value, names_from = si_or_ci) %>%
  mutate(stage_adjusted = ifelse(stage == "Primary", 0.65, 2.2)) %>%
  mutate(stage = ifelse(stage == "Primary", "O → P", "P → M"))
primary_stage_order <- c("O → P", "P → M")
AR_selection_data_Primary_Metastasis$stage <- factor(AR_selection_data_Primary_Metastasis$stage, levels = primary_stage_order)

scientific <- function(x){ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))}

# Safely wrap variant names in italic() for main dataset
selection_data_Early_Late_Metastasis <- selection_data_Early_Late_Metastasis %>%
  mutate(
    gene_clean = str_remove_all(variant_name, "italic\\(|\\)|'"),
    variant_label = paste0("italic('", gene_clean, "')")
  )

# For AR-specific data
AR_selection_data_Primary_Metastasis <- AR_selection_data_Primary_Metastasis %>%
  mutate(
    gene_clean = str_remove_all(variant_name, "italic\\(|\\)|'"),
    variant_label = paste0("italic('", gene_clean, "')")
  )

# Custom AR plot — already adjusted as in previous step
plot_ar <- ggplot(AR_selection_data_Primary_Metastasis, 
                  aes(x = stage_adjusted, y = si, color = stage)) + 
  geom_point(size = 1.5) + 
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) +
  facet_wrap(~ variant_label, scales = "free_y", ncol = 1, labeller = label_parsed) + 
  theme_bw() + 
  xlab("") + 
  ylab("Scaled selection coefficient") +
  theme(
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.spacing.x = unit(0.5, "cm"),
    legend.margin = margin(t = 10, l = 40), 
    legend.box.spacing = unit(2, "cm"),
    axis.text.x = element_blank(), 
    axis.ticks.x.bottom = element_blank(),
    text = element_text(size = 12),
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank(),  
    panel.grid.major.y = element_line(),   
    panel.grid.minor.y = element_line()) +
  scale_y_continuous(labels = scientific) +
  coord_cartesian(xlim = c(-0.15, 2.45)) +
  scale_color_manual(values = c("purple", "green")) +
  scale_x_continuous(
    breaks = AR_selection_data_Primary_Metastasis$stage_adjusted,
    labels = AR_selection_data_Primary_Metastasis$stage) + 
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

# Helper to generate individual plots for each gene
get_gene_plot <- function(gene_name) {
  this_data <- selection_data_Early_Late_Metastasis %>% filter(variant_name == gene_name)
  ggplot(this_data, aes(x = stage_adjusted, y = si, color = stage)) + 
    geom_point(size = 1.5) + 
    geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) +
    facet_wrap(~ variant_label, scales = "free_y", ncol = 1, labeller = label_parsed) + 
    theme_bw() + 
    xlab("") + ylab("Scaled selection coefficient") +
    theme(
      legend.position = "none", 
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.spacing.x = unit(0.5, "cm"),
      legend.margin = margin(t = 10),
      legend.box.spacing = unit(2, "cm"),
      axis.text.x = element_blank(), 
      axis.ticks.x.bottom = element_blank(),
      text = element_text(size = 12),
      panel.grid.major.x = element_blank(),  
      panel.grid.minor.x = element_blank(),  
      panel.grid.major.y = element_line(),   
      panel.grid.minor.y = element_line()
    ) +
    scale_y_continuous(labels = scientific) +
    coord_cartesian(xlim = c(-0.15, 2.45)) +
    scale_color_manual(values = c("red", "blue", "green")) +
    scale_x_continuous(
      breaks = this_data$stage_adjusted,
      labels = this_data$stage) + 
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
}

# Generate 15 plots (excluding AR)
plots_15 <- lapply(variant_order[variant_order != "AR"], get_gene_plot)

# Add the custom AR plot
all_16 <- c(plots_15, list(plot_ar))

# Combine into 4×4 grid with common y-axis label
combined <- wrap_plots(all_16, ncol = 4, guides = "collect") +
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    theme = theme(plot.margin = margin(5.5, 5.5, 5.5, 20),  # widen left margin for space
                  text = element_text(size = 12))) & 
  theme(axis.title.y = element_blank(),  # remove individual y labels
        plot.tag.position = "bottom",
        legend.position = "bottom")

g <- patchwork::patchworkGrob(combined)

# Open PNG device (can also use pdf(), svg(), etc.)
png("Fig_S3_custom_AR.png", width = 7, height = 9, units = "in", res = 300)

# Draw the plot
grid.draw(g)

# Overlay the y-axis label
grid.text("Scaled selection coefficient", 
          x = unit(0.02, "npc"), 
          y = unit(0.5, "npc"),
          rot = 90, 
          just = "center",
          gp = gpar(fontsize = 14))

# Close the device
dev.off()

#END


