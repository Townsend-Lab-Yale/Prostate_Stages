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


scientific <- function(x){ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))}


# selecting necessary data
selection_data_Early_Late_Metastasis <- rbindlist(cesa_samples_by_groups$selection)

# reformatting data set
selection_data_Early_Late_Metastasis <- selection_data_Early_Late_Metastasis |> 
  select(variant_name, starts_with("si"), starts_with("ci")) |>
  pivot_longer(cols = -variant_name, names_to = "data_type") |>
  mutate(stage = stringr::word(string = data_type, sep = "_",start = -1)) |>
  mutate(variant_name = stringr::str_remove(variant_name, "\\.1")) |>
  mutate(si_or_ci = stringr::word(string = data_type, sep = "_",start = 1, end=3)) |>
  mutate( si_or_ci = case_when(is.na(si_or_ci) ~ "si", TRUE ~ si_or_ci)) |>
  mutate (value = case_when (is.na(value)~0, TRUE~value))

# pivoting data set to create columns for gene, stage, si, and CIs
selection_data_Early_Late_Metastasis <- selection_data_Early_Late_Metastasis|> 
  select(-data_type) |>
  pivot_wider(values_from = value, names_from = si_or_ci)

# defining stages to be plotted
selection_data_Early_Late_Metastasis$stage <- factor(selection_data_Early_Late_Metastasis$stage, levels = c("Early","Late","Metastasis"))

### Making the Figure:

# Rename a specific word in the Name column
selection_data_Early_Late_Metastasis$stage <- sub("Early", "Lower-grade", selection_data_Early_Late_Metastasis$stage)
selection_data_Early_Late_Metastasis$stage <- sub("Late", "Higher-grade", selection_data_Early_Late_Metastasis$stage)
selection_data_Early_Late_Metastasis$stage <- sub("Metastasis_Lower-grade", "Metastasis", selection_data_Early_Late_Metastasis$stage)

variant_order <- c("CUL3", "SPOP", "PIK3CA", "AKT1", "ATM", "KMT2C", "KMT2D", "FOXA1", "APC", "ROCK1", "RHOA", "MUC16", "TP53", "CTNNB1", "PIK3CB", "AR") 
stage_order <- c("Lower-grade", "Higher-grade", "Metastasis")
selection_data_Early_Late_Metastasis$stage <- factor(selection_data_Early_Late_Metastasis$stage, levels = stage_order)

library(scales)
library(stringr)
library(dplyr)
library(ggplot2)

distance1 <- 1.1  # Distance between "Lower-grade" and "Higher-grade"
distance2 <- 1.0  # Distance between "Metastasis" and "Higher-grade"


selection_data_Early_Late_Metastasis$stage_adjusted <- ifelse(selection_data_Early_Late_Metastasis$stage == "Lower-grade", 0.1,
                        ifelse(selection_data_Early_Late_Metastasis$stage == "Higher-grade", 0.1 + distance1,
                               ifelse(selection_data_Early_Late_Metastasis$stage == "Metastasis", 0.1 + distance1 + distance2,
                                             NA)))										 
  
Figure_3_stage <- ggplot(selection_data_Early_Late_Metastasis, aes(x = stage_adjusted, y = si, color = stage, linetype = stage)) + 
    geom_point(size = 1.5) + 
	geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), width = .5) +
    facet_wrap(~ factor(variant_name, levels = variant_order), scales = "free_y", ncol = 4) + 
	theme_bw() + xlab("") + ylab("Cancer effect size") +
	theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_blank(), legend.text = element_text(size = 12)) +
	  scale_y_continuous(labels = scientific) +
	  theme(text = element_text(size = 12)) +  
	  expand_limits(y = 0) +
	  scale_linetype_manual(values = c(rep("solid", 1), rep("solid", 1), rep("solid", 1))) +
	  scale_color_manual(values = c("red", "blue", "green")) +
	  scale_x_continuous(breaks = selection_data_Early_Late_Metastasis$stage_adjusted, labels = selection_data_Early_Late_Metastasis$stage)
 
ggsave("Figure_3_stage.png", plot = Figure_3_stage, width = 8, height = 9)


#END
