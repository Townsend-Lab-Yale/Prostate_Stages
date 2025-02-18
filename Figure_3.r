library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg19)
library(MutationalPatterns)
library(RColorBrewer)
library(ggrepel)
library(readr)


setwd("C:/Moein/projects/prostate_stages/PRAD_files")

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
Early_groups <- cesa$samples[Gleason == "Early" & coverage == "exome", unique(Unique_Patient_Identifier)]
Late_groups <- cesa$samples[Gleason == "Late" & coverage == "exome", unique(Unique_Patient_Identifier)]
Metastasis_groups <- cesa$samples[Gleason == "Metastasis" & coverage == "exome", unique(Unique_Patient_Identifier)]

cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa, covariates = "PRAD", samples = Early_groups, save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "PRAD", samples = Late_groups, save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "PRAD", samples = Metastasis_groups, save_all_dndscv_output = T)


selected_genes <- c("SPOP", "FOXA1", "AR", "PIK3CA", "PIK3CB", "TP53", "ROCK1", "RHOA", "AKT1", "ATM", "CUL3",
                    "APC", "CTNNB1", "MUC16", "KMT2C", "KMT2D")
# Create compound variant table ----
# Get consensus coverage across whichever samples you want to include.
# Here, we use all WES/TGS, but you could choose to exclude some if they don't cover the genes of interest well.
all_cov = c(cesa_samples_by_groups$coverage_ranges$exome, cesa_samples_by_groups$coverage_ranges$targeted, cesa_samples_by_groups$coverage_ranges$genome)

# Exclude "exome", since typically "exome+" is what's applicable
all_cov = all_cov[! names(all_cov) == 'exome'] 
all_cov = Reduce(GenomicRanges::intersect, all_cov)


variants <- select_variants(cesa_samples_by_groups, genes = selected_genes, gr = all_cov)


# Further filter variants table based on COSMIC oncogene/TSG classification (exclude nonrecurrent except nonsense for TSGs).
top_TP53 <- variants[gene == "TP53" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_SPOP <- variants[gene == "SPOP" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_ATM <- variants[gene == "ATM" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_CUL3 <- variants[gene == "CUL3" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_APC <- variants[gene == "APC" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_MUC16 <- variants[gene == "MUC16" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_KMT2C <- variants[gene == "KMT2C" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_KMT2D <- variants[gene == "KMT2D" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]

top_PIK3CA <- variants[gene == "PIK3CA" & (maf_prevalence >1)]
top_FOXA1 <- variants[gene == "FOXA1" & (maf_prevalence >1)]
top_ROCK1 <- variants[gene == "ROCK1" & (maf_prevalence >1)]
top_RHOA <- variants[gene == "RHOA" & (maf_prevalence >1)]
top_CTNNB1 <- variants[gene == "CTNNB1" & (maf_prevalence > 1)]
top_PIK3CB <- variants[gene == "PIK3CB" & (maf_prevalence >1)]
top_AR <- variants[gene == "AR" & (maf_prevalence >1)]
top_AKT1 <- variants[gene == "AKT1" & (maf_prevalence >1)]

for_comp <- rbind(top_TP53, 
                  top_SPOP, 
                  top_ATM, 
                  top_CUL3, 
                  top_APC, 
                  top_MUC16, 
                  top_KMT2C, 
                  top_KMT2D, 
                  top_PIK3CA,
                  top_FOXA1,
                  top_ROCK1,
                  top_RHOA,
                  top_CTNNB1,
                  top_PIK3CB,
                  top_AR,
                  top_AKT1)

# Define compound variants to find cancer effect sizes at the gene level and not for individual variants
compound <- define_compound_variants(cesa = cesa_samples_by_groups, variant_table = for_comp, by = "gene", merge_distance = Inf)


RefCDS = ces.refset.hg19$RefCDS
dndscv_gene_names <- cesa_samples_by_groups$gene_rates$gene
nsyn_sites = sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

# selecting mutation rate data for samples in Early_groups 
samples_in_Early_groups <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_1$annotmuts$sampleID ))

# selecting mutation rate data for samples in Late_groups
samples_in_Late_groups <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_2$annotmuts$sampleID ))

# selecting mutation rate data for samples in Metastasis_groups
samples_in_Metastasis_groups <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_3$annotmuts$sampleID ))

library(tidyverse)
### creating a data frame with mutation rate data for Late_groups and Metastasis_groups
mut_rate_df <- tibble(gene = cesa_samples_by_groups$dNdScv_results$rate_grp_2$genemuts$gene_name,
                      exp_Late_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv,
                      exp_Metastasis_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_3$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df %>% 
  mutate(Late_mu = (exp_Late_mu / n_syn_sites) / samples_in_Late_groups) %>%
  mutate(Metastasis_mu = (exp_Metastasis_mu / n_syn_sites) / samples_in_Metastasis_groups) %>%
  mutate(cancer_greater = Metastasis_mu > Late_mu) -> 
  mut_rate_df

# defining rate 1 and rate 2 as mutation rates for Late_groups and Metastasis_groups
rate_1 <- mut_rate_df|>
  select(gene, Late_mu)
rate_2 <- mut_rate_df|>
  select(gene, Metastasis_mu)

# change in mutation rate across stages
mut_rate_df <- mut_rate_df %>% 
  select(gene, Late_mu, Metastasis_mu) %>% 
  mutate(p_1 = Late_mu / Metastasis_mu) %>% 
  mutate(p_2 = 1 - p_1)
  
# saving "last" gene mutation rates into separate data frame, "last" rates meaning from last stage Metastasis_mu
set_cancer_rates <- mut_rate_df %>%
  select(gene, Metastasis_mu) %>%
  data.table::setDT()

# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# setting gene rates to highest rates from Metastasis_mu
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = set_cancer_rates, missing_genes_take_nearest = T) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
cesa_samples_by_groups <- trinuc_mutation_rates(cesa = cesa_samples_by_groups, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)


source("new_sequential_lik.R")

for(comp_ind in 1:length(compound)){
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  these_props <- mut_rate_df[mut_rate_df$gene == this_gene,c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  
  cesa_samples_by_groups <- ces_variant(cesa = cesa_samples_by_groups, variants = this_comp, model = sequential_lik_dev, 
                                        ordering_col = 'Gleason', ordering = c('Late', 'Metastasis'), 
                                        lik_args = list(sequential_mut_prop = these_props), run_name = this_gene)
  
}


scientific <- function(x){ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))}

# selecting necessary data

selection_data <- rbindlist(cesa_samples_by_groups$selection)

# reformatting data set
selection_data <- selection_data |> 
  select(variant_name, starts_with("si"), starts_with("ci")) |>
  pivot_longer(cols = -variant_name, names_to = "data_type") |>
  mutate(stage = stringr::word(string = data_type, sep = "_",start = -1)) |>
  mutate(variant_name = stringr::str_remove(variant_name, "\\.1")) |>
  mutate(si_or_ci = stringr::word(string = data_type, sep = "_",start = 1, end=3)) |>
  mutate( si_or_ci = case_when(is.na(si_or_ci) ~ "si", TRUE ~ si_or_ci)) |>
  mutate (value = case_when (is.na(value)~0, TRUE~value))

# pivoting data set to create columns for gene, stage, si, and CIs
selection_data <- selection_data|> 
  select(-data_type) |>
  pivot_wider(values_from = value, names_from = si_or_ci)

# defining stages to be plotted
selection_data$stage <- factor(selection_data$stage, levels = c("Late","Metastasis"))


### creating a data frame with mutation rate data for Early_groups and Metastasis_groups
rm (cesa_samples_by_groups)

cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa, covariates = "PRAD", samples = Early_groups, save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "PRAD", samples = Late_groups, save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "PRAD", samples = Metastasis_groups, save_all_dndscv_output = T)


mut_rate_df <- tibble(gene = cesa_samples_by_groups$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_Early_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_Metastasis_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_3$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df %>% 
  mutate(Early_mu = (exp_Early_mu / n_syn_sites) / samples_in_Early_groups) %>%
  mutate(Metastasis_mu = (exp_Metastasis_mu / n_syn_sites) / samples_in_Metastasis_groups) %>%
  mutate(cancer_greater = Metastasis_mu > Early_mu) -> 
  mut_rate_df

# defining rate 1 and rate 2 as mutation rates for Early_groups and Metastasis_groups
rate_1 <- mut_rate_df|>
  select(gene, Early_mu)
rate_2 <- mut_rate_df|>
  select(gene, Metastasis_mu)

# change in mutation rate across stages
mut_rate_df <- mut_rate_df %>% 
  select(gene, Early_mu, Metastasis_mu) %>% 
  mutate(p_1 = Early_mu / Metastasis_mu) %>% 
  mutate(p_2 = 1 - p_1)
  
# saving "last" gene mutation rates into separate data frame, "last" rates meaning from last stage Metastasis_mu
set_cancer_rates <- mut_rate_df %>%
  select(gene, Metastasis_mu) %>%
  data.table::setDT()

# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# setting gene rates to highest rates from Metastasis_mu
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = set_cancer_rates, missing_genes_take_nearest = T) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
cesa_samples_by_groups <- trinuc_mutation_rates(cesa = cesa_samples_by_groups, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)

source("new_sequential_lik.R")

for(comp_ind in 1:length(compound)){
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  these_props <- mut_rate_df[mut_rate_df$gene == this_gene,c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  
  cesa_samples_by_groups <- ces_variant(cesa = cesa_samples_by_groups, variants = this_comp, model = sequential_lik_dev, 
                                        ordering_col = 'Gleason', ordering = c('Early', 'Metastasis'), 
                                        lik_args = list(sequential_mut_prop = these_props), run_name = this_gene)
  
}

# selecting necessary data
selection_data_Early_Metastasis <- rbindlist(cesa_samples_by_groups$selection)

# reformatting data set
selection_data_Early_Metastasis <- selection_data_Early_Metastasis |> 
  select(variant_name, starts_with("si"), starts_with("ci")) |>
  pivot_longer(cols = -variant_name, names_to = "data_type") |>
  mutate(stage = stringr::word(string = data_type, sep = "_",start = -1)) |>
  mutate(variant_name = stringr::str_remove(variant_name, "\\.1")) |>
  mutate(si_or_ci = stringr::word(string = data_type, sep = "_",start = 1, end=3)) |>
  mutate( si_or_ci = case_when(is.na(si_or_ci) ~ "si", TRUE ~ si_or_ci)) |>
  mutate (value = case_when (is.na(value)~0, TRUE~value))

# pivoting data set to create columns for gene, stage, si, and CIs
selection_data_Early_Metastasis <- selection_data_Early_Metastasis|> 
  select(-data_type) |>
  pivot_wider(values_from = value, names_from = si_or_ci)

# defining stages to be plotted
selection_data_Early_Metastasis$stage <- factor(selection_data_Early_Metastasis$stage, levels = c("Early","Metastasis"))


### Making the Figure:

# Rename a specific word in the Name column
selection_data_Early_Metastasis$stage <- sub("Early", "E → L", selection_data_Early_Metastasis$stage)
selection_data_Early_Metastasis$stage <- sub("Metastasis", "L → M", selection_data_Early_Metastasis$stage)
selection_data$stage <- sub("Late", "E → H", selection_data$stage)
selection_data$stage <- sub("Metastasis", "H → M", selection_data$stage)

combined_selection <- rbind(selection_data_Early_Metastasis, selection_data)

variant_order <- c("CUL3", "SPOP", "PIK3CA", "AKT1", "ATM", "KMT2C", "KMT2D", "FOXA1", "APC", "ROCK1", "RHOA", "MUC16", "TP53", "CTNNB1", "PIK3CB", "AR") 
stage_order <- c("E → L", "E → H", "L → M", "H → M")
combined_selection$stage <- factor(combined_selection$stage, levels = stage_order)

library(scales)
library(stringr)
library(dplyr)
library(ggplot2)

distance1 <- 1.1  # Distance between "Lower-risk" and "Higher-risk"
distance2 <- 1.0  # Distance between "Metastasis_Lower-risk" and "Lower-risk"
distance3 <- 1.6  # Distance between "Metastasis_Higher-risk" and "Higher-risk"


combined_selection$stage_adjusted <- ifelse(combined_selection$stage == "E → L", 0.1,
                        ifelse(combined_selection$stage == "E → H", 0.1 + distance1,
                               ifelse(combined_selection$stage == "L → M", 0.1 + distance1 + distance2,
                                      ifelse(combined_selection$stage == "H → M", 0.1 + distance1 + distance3,
                                             NA))))										 
  
Figure_3 <- ggplot(combined_selection, aes(x = stage_adjusted, y = si, color = stage, linetype = stage)) + 
    geom_point(size = 1.5) + 
	geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) +
    facet_wrap(~ factor(variant_name, levels = variant_order), scales = "free_y", ncol = 4) + 
	theme_bw() + xlab("") + ylab("Cancer effect size") +
	theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank(), legend.text = element_text(size = 12)) +
	  scale_y_continuous(labels = scientific) +
	  theme(text = element_text(size = 12)) +  
	  expand_limits(y = 0) +
	  scale_linetype_manual(values = c(rep("solid", 1), rep("solid", 1), rep("twodash", 1), rep("twodash", 1))) +
	  scale_color_manual(values = c("red", "blue", "red", "blue")) +
	  scale_x_continuous(breaks = combined_selection$stage_adjusted, labels = combined_selection$stage)
 
ggsave("Figure_3.png", plot = Figure_3, width = 8, height = 9)

#End
