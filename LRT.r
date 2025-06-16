library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg19)
library(MutationalPatterns)
library(RColorBrewer)
library(ggrepel)
library(readr)

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

Early_Late_groups <- cesa$samples[Gleason %in% c("Early", "Late") & coverage == "exome"]
Early_Metastasis_groups <- cesa$samples[Gleason %in% c("Early", "Metastasis") & coverage == "exome"]
Late_Metastasis_groups <- cesa$samples[Gleason %in% c("Late", "Metastasis") & coverage == "exome"]

# Estimate neutral gene mutation rates with dNdScv using tissue-specific mutation rate covariates:

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

library(tidyverse)
### creating a data frame with mutation rate data for Late_groups and Metastasis_groups
mut_rate_df <- tibble(gene = cesa_samples_by_groups$dNdScv_results$rate_grp_2$genemuts$gene_name,
                      exp_Late_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv,
                      exp_Metastasis_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_3$genemuts$exp_syn_cv)

# Add n_syn_sites column to mut_rate_df:
mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

# Mu = ((expected mu)/(number of synonymous sites))/(samples)
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

# saving gene mutation rates into separate data frame:
Late_rate <- mut_rate_df %>%
  select(gene, rate = Late_mu) %>%
  data.table::setDT()
Meta_rate <- mut_rate_df %>%
  select(gene, rate = Metastasis_mu) %>%
  data.table::setDT()

# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# setting gene rates for Late and Metastasis:
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = Late_rate, missing_genes_take_nearest = T, samples = cesa$samples[Gleason=="Late"]) 
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = Meta_rate, missing_genes_take_nearest = T, samples = cesa$samples[Gleason=="Metastasis"]) 

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
  these_props <- mut_rate_df[mut_rate_df$gene == this_gene,c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  
  cesa_samples_by_groups <- ces_variant(cesa = cesa_samples_by_groups, variants = this_comp, model = sequential_lik_dev, 
                                        ordering_col = 'Gleason', ordering = c('Late', 'Metastasis'), samples = c(Late_groups_all, Metastasis_groups_all), 
                                        lik_args = list(sequential_mut_prop = these_props), run_name = this_gene)
  
}


scientific <- function(x){ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))}

# selecting necessary data
selection_data_Late_Meta_seq_model <- rbindlist(cesa_samples_by_groups$selection)


#Clear variant effect output
cesa_samples_by_groups <- clear_effect_output(cesa_samples_by_groups)

for (comp_ind in 1:length(compound)) {
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  
  cesa_samples_by_groups <- ces_variant(
    cesa = cesa_samples_by_groups,
    variants = this_comp, samples = c(Late_groups_all, Metastasis_groups_all),
    run_name = this_gene
  )
}


# selecting necessary data
selection_data_Late_Meta_default_model <- rbindlist(cesa_samples_by_groups$selection)


# Extract relevant columns from both datasets:
selection_data_Late_Meta_seq_model <- selection_data_Late_Meta_seq_model %>%
  select(variant_name, loglikelihood) %>%
  rename(loglikelihood_seq_model = loglikelihood)

selection_data_Late_Meta_default_model <- selection_data_Late_Meta_default_model %>%
  select(variant_name, loglikelihood) %>%
  rename(loglikelihood_default_model = loglikelihood)

# Combine the dataset:
loglik_Late_Meta <- selection_data_Late_Meta_seq_model %>%
  left_join(selection_data_Late_Meta_default_model, by = "variant_name")


loglik_Late_Meta <- loglik_Late_Meta %>%
  mutate(
    loglikelihood_seq_model = ifelse(loglikelihood_seq_model > 1e5, NA, loglikelihood_seq_model),
    loglikelihood_default_model = ifelse(loglikelihood_default_model > 1e5, NA, loglikelihood_default_model),
    LRT_stat = 2 * (loglikelihood_seq_model - loglikelihood_default_model),
    p_value = pchisq(LRT_stat, df = 1, lower.tail = FALSE),
    significant = ifelse(p_value < 0.05, TRUE, FALSE)
  )


### Early_groups & Metastasis_groups

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

# saving gene mutation rates into separate data frame:
Early_rate <- mut_rate_df %>%
  select(gene, rate = Early_mu) %>%
  data.table::setDT()
Meta_rate <- mut_rate_df %>%
  select(gene, rate = Metastasis_mu) %>%
  data.table::setDT()

# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# setting gene rates for Late and Metastasis:
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = Early_rate, missing_genes_take_nearest = T, samples = cesa$samples[Gleason=="Early"]) 
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = Meta_rate, missing_genes_take_nearest = T, samples = cesa$samples[Gleason=="Metastasis"]) 


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
  these_props <- mut_rate_df[mut_rate_df$gene == this_gene,c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  
  cesa_samples_by_groups <- ces_variant(cesa = cesa_samples_by_groups, variants = this_comp, model = sequential_lik_dev, 
                                        ordering_col = 'Gleason', ordering = c('Early', 'Metastasis'), samples = c(Early_groups_all, Metastasis_groups_all),
                                        lik_args = list(sequential_mut_prop = these_props), run_name = this_gene)
  
}


selection_data_Early_Meta_seq_model <- rbindlist(cesa_samples_by_groups$selection)


#Clear variant effect output
cesa_samples_by_groups <- clear_effect_output(cesa_samples_by_groups)

for (comp_ind in 1:length(compound)) {
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  
  cesa_samples_by_groups <- ces_variant(
    cesa = cesa_samples_by_groups,
    variants = this_comp, samples = c(Early_groups_all, Metastasis_groups_all),
    run_name = this_gene
  )
}


# selecting necessary data
selection_data_Early_Meta_default_model <- rbindlist(cesa_samples_by_groups$selection)


# Extract relevant columns from both datasets:
selection_data_Early_Meta_seq_model <- selection_data_Early_Meta_seq_model %>%
  select(variant_name, loglikelihood) %>%
  rename(loglikelihood_seq_model = loglikelihood)

selection_data_Early_Meta_default_model <- selection_data_Early_Meta_default_model %>%
  select(variant_name, loglikelihood) %>%
  rename(loglikelihood_default_model = loglikelihood)

# Combine the dataset:
loglik_Early_Meta <- selection_data_Early_Meta_seq_model %>%
  left_join(selection_data_Early_Meta_default_model, by = "variant_name")


loglik_Early_Meta <- loglik_Early_Meta %>%
  mutate(
    loglikelihood_seq_model = ifelse(loglikelihood_seq_model > 1e5, NA, loglikelihood_seq_model),
    loglikelihood_default_model = ifelse(loglikelihood_default_model > 1e5, NA, loglikelihood_default_model),
    LRT_stat = 2 * (loglikelihood_seq_model - loglikelihood_default_model),
    p_value = pchisq(LRT_stat, df = 1, lower.tail = FALSE),
    significant = ifelse(p_value < 0.05, TRUE, FALSE)
  )

#End
