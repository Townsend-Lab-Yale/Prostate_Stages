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
Early_Late_groups <- cesa$samples[Gleason %in% c("Early", "Late") & coverage == "exome"]
Early_Metastasis_groups <- cesa$samples[Gleason %in% c("Early", "Metastasis") & coverage == "exome"]
Late_Metastasis_groups <- cesa$samples[Gleason %in% c("Late", "Metastasis") & coverage == "exome"]

cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa, covariates = "PRAD", samples = Early_groups, save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "PRAD", samples = Late_groups, save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "PRAD", samples = Metastasis_groups, save_all_dndscv_output = T)


selected_genes <- c("SPOP", "FOXA1", "AR", "PIK3CA", "PIK3CB", "TP53", "ROCK1", "RHOA", "AKT1", "ATM", "CUL3",
                    "APC", "CTNNB1", "PTEN", "KMT2C", "KMT2D")


### selection intensity of Late & Metastasis using mutation rate of Metastasis:

library(tidyverse)
### creating a data frame with mutation rate data for Metastasis_groups
mut_rate_df_Meta <- tibble(gene = cesa_samples_by_groups$gene_rates$gene,
                      Metastasis_gene_Rates = cesa_samples_by_groups$gene_rates$rate_grp_3)

set_cancer_rates_Meta <- mut_rate_df_Meta %>%
  select(gene, Metastasis_gene_Rates) %>%
  data.table::setDT()


### creating a data frame with mutation rate data for Late_groups
mut_rate_df_Late <- tibble(gene = cesa_samples_by_groups$gene_rates$gene,
                      Late_gene_Rates = cesa_samples_by_groups$gene_rates$rate_grp_2)

set_cancer_rates_Late <- mut_rate_df_Late %>%
  select(gene, Late_gene_Rates) %>%
  data.table::setDT()

### creating a data frame with mutation rate data for Early_groups
mut_rate_df_Early <- tibble(gene = cesa_samples_by_groups$gene_rates$gene,
                      Early_gene_Rates = cesa_samples_by_groups$gene_rates$rate_grp_1)

set_cancer_rates_Early <- mut_rate_df_Early %>%
  select(gene, Early_gene_Rates) %>%
  data.table::setDT()


# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# setting gene rates to Metastasis_rate
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = set_cancer_rates_Meta, missing_genes_take_nearest = T) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
cesa_samples_by_groups <- trinuc_mutation_rates(cesa = cesa_samples_by_groups, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)


# defining compound variants
compound <- define_compound_variants(cesa = cesa_samples_by_groups, 
                                     variant_table = cesa_samples_by_groups$variants |>
                                       filter(intergenic == F, gene %in% selected_genes),
                                     by = "gene", merge_distance = Inf)



for (comp_ind in 1:length(compound)) {
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  
  cesa_samples_by_groups <- ces_variant(
    cesa = cesa_samples_by_groups,
    variants = this_comp, samples = Late_Metastasis_groups,
    run_name = this_gene
  )
}


scientific <- function(x){ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))}

# selecting necessary data

selection_Late&Meta_rate_Meta <- rbindlist(cesa_samples_by_groups$selection)

# reformatting data set
selection_Late&Meta_rate_Meta <- selection_Late&Meta_rate_Meta |>
  select(variant_name, starts_with("selection"), starts_with("log"), starts_with("ci")) |>
  mutate(variant_name = stringr::str_remove(variant_name, "\\.1")) |>
  mutate(across(-variant_name, ~replace_na(., 0)))


data.table::fwrite(selection_Late&Meta_rate_Meta, file = "selection_Late&Meta_rate_Meta.txt", sep = "\t")


### selection intensity of Early & Metastasis using mutation rate of Metastasis:

#Clear variant effect output
cesa_samples_by_groups <- clear_effect_output(cesa_samples_by_groups)


for (comp_ind in 1:length(compound)) {
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  
  cesa_samples_by_groups <- ces_variant(
    cesa = cesa_samples_by_groups,
    variants = this_comp, samples = Early_Metastasis_groups,
    run_name = this_gene
  )
}


scientific <- function(x){ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))}

# selecting necessary data

selection_data_Early_Metastasis <- rbindlist(cesa_samples_by_groups$selection)

# reformatting data set
selection_data_Early_Metastasis <- selection_data_Early_Metastasis |>
  select(variant_name, starts_with("selection"), starts_with("log"), starts_with("ci")) |>
  mutate(variant_name = stringr::str_remove(variant_name, "\\.1")) |>
  mutate(across(-variant_name, ~replace_na(., 0)))


data.table::fwrite(selection_data_Early_Metastasis, file = "selection_data_Si_Early_Meta_rate_Meta.txt", sep = "\t")


### selection intensity of Late & Metastasis using mutation rate of Late:

#Clear variant effect output
cesa_samples_by_groups <- clear_effect_output(cesa_samples_by_groups)

# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# trinucleotide mutation rates 
cesa_samples_by_groups <- clear_trinuc_rates_and_signatures(cesa_samples_by_groups)

# setting gene rates to Late_rate
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = set_cancer_rates_Late, missing_genes_take_nearest = T) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
cesa_samples_by_groups <- trinuc_mutation_rates(cesa = cesa_samples_by_groups, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)


# defining compound variants
compound <- define_compound_variants(cesa = cesa_samples_by_groups, 
                                     variant_table = cesa_samples_by_groups$variants |>
                                       filter(intergenic == F, gene %in% selected_genes),
                                     by = "gene", merge_distance = Inf)



for (comp_ind in 1:length(compound)) {
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  
  cesa_samples_by_groups <- ces_variant(
    cesa = cesa_samples_by_groups,
    variants = this_comp, samples = Late_Metastasis_groups,
    run_name = this_gene
  )
}



# selecting necessary data

selection_data_Si_Late_Meta_rate_Late <- rbindlist(cesa_samples_by_groups$selection)

# reformatting data set
selection_data_Si_Late_Meta_rate_Late <- selection_data_Si_Late_Meta_rate_Late |>
  select(variant_name, starts_with("selection"), starts_with("log"), starts_with("ci")) |>
  mutate(variant_name = stringr::str_remove(variant_name, "\\.1")) |>
  mutate(across(-variant_name, ~replace_na(., 0)))


data.table::fwrite(selection_data_Si_Late_Meta_rate_Late, file = "selection_Late&Meta_rate_Late.txt", sep = "\t")


### selection intensity of Early & Metastasis using mutation rate of Early:

#Clear variant effect output
cesa_samples_by_groups <- clear_effect_output(cesa_samples_by_groups)

# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# trinucleotide mutation rates 
cesa_samples_by_groups <- clear_trinuc_rates_and_signatures(cesa_samples_by_groups)

# setting gene rates to Late_rate
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = set_cancer_rates_Early, missing_genes_take_nearest = T) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
cesa_samples_by_groups <- trinuc_mutation_rates(cesa = cesa_samples_by_groups, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)


# defining compound variants
compound <- define_compound_variants(cesa = cesa_samples_by_groups, 
                                     variant_table = cesa_samples_by_groups$variants |>
                                       filter(intergenic == F, gene %in% selected_genes),
                                     by = "gene", merge_distance = Inf)



for (comp_ind in 1:length(compound)) {
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  
  cesa_samples_by_groups <- ces_variant(
    cesa = cesa_samples_by_groups,
    variants = this_comp, samples = Early_Metastasis_groups,
    run_name = this_gene
  )
}


# selecting necessary data

selection_Early&Meta_rate_Early <- rbindlist(cesa_samples_by_groups$selection)

# reformatting data set
selection_Early&Meta_rate_Early <- selection_Early&Meta_rate_Early |>
  select(variant_name, starts_with("selection"), starts_with("log"), starts_with("ci")) |>
  mutate(variant_name = stringr::str_remove(variant_name, "\\.1")) |>
  mutate(across(-variant_name, ~replace_na(., 0)))


data.table::fwrite(selection_Early&Meta_rate_Early, file = "selection_Early&Meta_rate_Early.txt", sep = "\t")


#End
