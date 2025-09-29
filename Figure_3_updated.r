library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg19)
library(MutationalPatterns)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(scales)

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

#Defining group
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
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = Early_rate, missing_genes_take_nearest = T, samples = cesa$samples[Gleason=="Early"]) 
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = Late_rate, missing_genes_take_nearest = T, samples = cesa$samples[Gleason=="Late"]) 
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = Meta_rate, missing_genes_take_nearest = T, samples = cesa$samples[Gleason=="Metastasis"]) 


# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
cesa_samples_by_groups <- trinuc_mutation_rates(cesa = cesa_samples_by_groups, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)

threestage_final <- cesa_samples_by_groups
saveRDS(threestage_final, file="threestage_final.rds")

#groups:
mut_rate_stageless <- data.frame(gene=cesa_samples_by_groups@mutrates$gene,
                                 mutation_rate=cesa_samples_by_groups@mutrates$rate_grp_1)

highlight <- rep(FALSE, length(mut_rate_stageless$gene))
highlight[c(16063, 6373, 1221, 12567, 17456)] <- TRUE

#scatter plot of gene-level mutation rates, prim vs met
selected_colors <- brewer.pal(n = length(highlight[highlight==TRUE]), name = 'Dark2')
selected_colors[2] <- "#027cd9"

###Figure_early_late:
threestage_final <- load_cesa("threestage_final.rds")

mut_rate_earlylate <- data.frame(gene=threestage_final@mutrates$gene,
                       early_rate=threestage_final@mutrates$rate_grp_1,
                       late_rate=threestage_final@mutrates$rate_grp_2)

mutrates_early_late_plot <- ggplot()+
  geom_point(data = mut_rate_earlylate[!highlight,], aes(x=early_rate, y=late_rate), size=1, color="lightcoral", alpha=0.5, shape = 16) +
  geom_point(data = mut_rate_earlylate[highlight,], aes(x=early_rate, y=late_rate), size=2.5, color=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate_earlylate[highlight,], aes(x=early_rate, y=late_rate, label=gene)) +
  theme_bw() +
  theme(panel.border = element_blank(), plot.title = element_text(hjust=0.5)) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  xlab("Mutation rate in low-grade tumors") +
  ylab("Mutation rate in high-grade tumors") +
  geom_smooth(method="lm", color="navyblue") + 
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(5e-7, 1e-6), expand = c(0, 0)) +
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6, 2e-6, 3e-6, 4e-6), expand = c(0, 0))

###Figure_early_met:
mut_rate_earlymet <- data.frame(gene=threestage_final@mutrates$gene,
                       early_rate=threestage_final@mutrates$rate_grp_1,
                       metastasis_rate=threestage_final@mutrates$rate_grp_3)

mutrates_early_metastasis_plot <- ggplot()+
  geom_point(data = mut_rate_earlymet[!highlight,], aes(x=early_rate, y=metastasis_rate), size=1, color="lightcoral", alpha=0.5, shape = 16) +
  geom_point(data = mut_rate_earlymet[highlight,], aes(x=early_rate, y=metastasis_rate), size=2.5, color=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate_earlymet[highlight,], aes(x=early_rate, y=metastasis_rate, label=gene)) +
  theme_bw() +
  theme(panel.border = element_blank(), plot.title = element_text(hjust=0.5)) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  xlab("Mutation rate in low-grade tumors") +
  ylab("Mutation rate in mCRPCs") +
  geom_smooth(method="lm", color="navyblue") + 
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(5e-7, 1e-6), expand = c(0, 0)) +
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6, 2e-6, 3e-6, 4e-6), expand = c(0, 0))

###Figure_late_met:
mut_rate_latemet <- data.frame(gene=threestage_final@mutrates$gene,
                       late_rate=threestage_final@mutrates$rate_grp_2,
                       metastasis_rate=threestage_final@mutrates$rate_grp_3)

mutrates_late_metastasis_plot <- ggplot()+
  geom_point(data = mut_rate_latemet[!highlight,], aes(x=late_rate, y=metastasis_rate), size=1, color="lightcoral", alpha=0.5, shape = 16) +
  geom_point(data = mut_rate_latemet[highlight,], aes(x=late_rate, y=metastasis_rate), size=2.5, color=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate_latemet[highlight,], aes(x=late_rate, y=metastasis_rate, label=gene)) +
  theme_bw() +
  theme(panel.border = element_blank(), plot.title = element_text(hjust=0.5)) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  xlab("Mutation rate in high-grade tumors") +
  ylab("Mutation rate in mCRPCs") +
  geom_smooth(method="lm", color="navyblue") + 
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(5e-7, 1e-6), expand = c(0, 0)) +
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6, 2e-6, 3e-6, 4e-6), expand = c(0, 0))

### Figure_separate_tumors:
##Fig_mut_rate_early
mut_rate_early <- data.frame(gene=cesa_samples_by_groups@mutrates$gene,
                                 mutation_rate_early=cesa_samples_by_groups@mutrates$rate_grp_1)


mutrates_early_plot <- ggplot() +
  geom_density(data=mut_rate_early, size = 1, aes(color="lightcoral", x=mutation_rate_early, y=..scaled..))+
  xlab("Mutation rate in low-grade tumors") + ylab("Density") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none",
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  theme(plot.margin = margin(0,0,0,26, "pt")) +
  scale_x_continuous(labels=scientific, limits=c(0, 1.25e-6), breaks=c(5e-7, 1e-6), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 1.05), breaks=c(0,0.2,0.4,0.6,0.8,1), expand = c(0, 0))

##Fig_mut_rate_late
mut_rate_late <- data.frame(gene=cesa_samples_by_groups@mutrates$gene,
                                 mutation_rate_late=cesa_samples_by_groups@mutrates$rate_grp_2)


mutrates_late_plot <- ggplot() +
  geom_density(data=mut_rate_late, size = 1, aes(color="lightcoral", x=mutation_rate_late, y=..scaled..))+
  xlab("Mutation rate in high-grade tumors") + ylab("Density") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none",
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  theme(plot.margin = margin(0,0,0,26, "pt")) +
  scale_x_continuous(labels=scientific, limits=c(0, 1.25e-6), breaks=c(5e-7, 1e-6), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 1.05), breaks=c(0,0.2,0.4,0.6,0.8,1), expand = c(0, 0))
 
##Fig_mut_rate_metastasis
mut_rate_metastasis <- data.frame(gene=cesa_samples_by_groups@mutrates$gene,
                                 mutation_rate_metastasis=cesa_samples_by_groups@mutrates$rate_grp_3)


mutrates_metastasis_plot <- ggplot() +
  geom_density(data=mut_rate_metastasis, size = 1, aes(color="lightcoral", x=mutation_rate_metastasis, y=..scaled..))+
  xlab("Mutation rate in mCRPCs") + ylab("Density") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none",
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  theme(plot.margin = margin(0,0,0,26, "pt")) +
  scale_x_continuous(labels=scientific, limits=c(0, 1.25e-6), breaks=c(5e-7, 1e-6), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 1.05), breaks=c(0,0.2,0.4,0.6,0.8,1), expand = c(0, 0))
 

###Primary_Metastasis:
Primary_Met <- read.delim("C:/Moein/projects/prostate_stages/PRAD_files/Primary_Met.txt")

##Create CESAnalysis and load data:
Pri_Met <- CESAnalysis(refset = "ces.refset.hg19")
Pri_Met <- load_maf(cesa = Pri_Met, maf = MAF1, maf_name = "armenia")
Pri_Met <- load_maf(cesa = Pri_Met, maf = MAF2, maf_name = "boutros", coverage = "genome")
Pri_Met <- load_maf(cesa = Pri_Met, maf = MAF3, maf_name = "tcga", coverage = "genome")
Pri_Met <- load_maf(cesa = Pri_Met, maf = MAF4, maf_name = "SU2C",, coverage = "exome",
                     covered_regions = "SureSelect_All_Exon_covered_regions.bed",
                     covered_regions_name = "SureSelect_V4", covered_regions_padding = 100)
Pri_Met <- load_maf(cesa = Pri_Met, maf = MAF5, maf_name = "341", coverage = "targeted", 
                     covered_regions = "msk_341_exons.bed", 
                     covered_regions_name = "MSK_IMPACT_341", covered_regions_padding = 10)
Pri_Met <- load_maf(cesa = Pri_Met, maf = MAF6, maf_name = "410", coverage = "targeted", 
                     covered_regions = "msk_410_exons.bed", 
                     covered_regions_name = "MSK_IMPACT_410", covered_regions_padding = 10)
Pri_Met <- load_maf(cesa = Pri_Met, maf = MAF7, maf_name = "468", coverage = "targeted", 
                     covered_regions = "msk_468_exons.bed", 
                     covered_regions_name = "MSK_IMPACT_468", covered_regions_padding = 10)

Pri_Met <- load_sample_data(Pri_Met, Primary_Met)

#defining groups for gene mutation rate using exome:
Primary_groups <- Pri_Met$samples[Primary_Met == "Primary" & coverage == "exome", unique(Unique_Patient_Identifier)]
Metastasis_groups <- Pri_Met$samples[Primary_Met == "Metastasis" & coverage == "exome", unique(Unique_Patient_Identifier)]

#Defining group
Primary_groups_all <- Pri_Met$samples[Primary_Met == "Primary", unique(Unique_Patient_Identifier)]
Metastasis_groups_all <- Pri_Met$samples[Primary_Met == "Metastasis", unique(Unique_Patient_Identifier)]

Pri_Met <- gene_mutation_rates(cesa = Pri_Met, covariates = "PRAD", samples = Primary_groups, save_all_dndscv_output = T)
Pri_Met <- gene_mutation_rates(cesa = Pri_Met, covariates = "PRAD", samples = Metastasis_groups, save_all_dndscv_output = T)

RefCDS = ces.refset.hg19$RefCDS
dndscv_gene_names <- Pri_Met$gene_rates$gene
nsyn_sites = sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

# selecting mutation rate data for samples in Early_groups 
samples_in_Primary_groups <- length(unique(Pri_Met$dNdScv_results$rate_grp_1$annotmuts$sampleID ))

# selecting mutation rate data for samples in Late_groups
samples_in_Metastasis_groups <- length(unique(Pri_Met$dNdScv_results$rate_grp_2$annotmuts$sampleID ))

library(tidyverse)
### creating a data frame with mutation rate data for Primary_group and Metastasis_groups
mut_rate_df <- tibble(gene = Pri_Met$dNdScv_results$rate_grp_2$genemuts$gene_name,
                      exp_Primary_mu = Pri_Met$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_Metastasis_mu = Pri_Met$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv)

# Add n_syn_sites column to mut_rate_df:
mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df %>% 
  mutate(Primary_mu = (exp_Primary_mu / n_syn_sites) / samples_in_Primary_groups) %>%
  mutate(Metastasis_mu = (exp_Metastasis_mu / n_syn_sites) / samples_in_Metastasis_groups) %>%
  mutate(cancer_greater = (Metastasis_mu > Primary_mu)) -> 
  mut_rate_df

# saving gene mutation rates into separate data frame:
Primary_rate <- mut_rate_df %>%
  select(gene, rate = Primary_mu) %>%
  data.table::setDT()
Meta_rate <- mut_rate_df %>%
  select(gene, rate = Metastasis_mu) %>%
  data.table::setDT()

# clear the gene rates in the cesa object 
Pri_Met <- clear_gene_rates(cesa = Pri_Met)

# setting gene rates for Early, Late and Metastasis:
Pri_Met <- set_gene_rates(cesa = Pri_Met, rates = Primary_rate, missing_genes_take_nearest = T, samples = Pri_Met$samples[Primary_Met == "Primary"]) 
Pri_Met <- set_gene_rates(cesa = Pri_Met, rates = Meta_rate, missing_genes_take_nearest = T, samples = Pri_Met$samples[Primary_Met == "Metastasis"]) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
Pri_Met  <- trinuc_mutation_rates(cesa = Pri_Met, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)

twostage_final <- Pri_Met
saveRDS(twostage_final, file="twostage_final.rds")

###Figure_prim_met: 
###mutrates_prim_met_plo
mut_rate <- data.frame(gene=twostage_final@mutrates$gene,
                                 prim_rate=twostage_final@mutrates$rate_grp_1,
                                 met_rate=twostage_final@mutrates$rate_grp_2)

mutrates_prim_met_plot <- ggplot()+
  geom_point(data = mut_rate[!highlight,], aes(x=prim_rate, y=met_rate), size=1, color="lightcoral", alpha=0.5, shape = 16) +
  geom_point(data = mut_rate[highlight,], aes(x=prim_rate, y=met_rate), size=2.5, colour=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate[highlight,], aes(x=prim_rate, y=met_rate, label=gene)) +
  theme_bw() +
  theme(panel.border = element_blank(), plot.title = element_text(hjust=0.5)) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  xlab("Mutation rate in primary tumors") +
  ylab("Mutation rate in mCRPCs") +
  geom_smooth(method="lm", color="navyblue") +
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(5e-7, 1e-6), expand = c(0, 0)) + 
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6, 2e-6, 3e-6, 4e-6), expand = c(0, 0))

####Fig_mut_rate_primary
mut_rate_primary <- data.frame(gene=twostage_final@mutrates$gene,
                                 mutation_rate_primary=twostage_final@mutrates$rate_grp_1)

mutrates_primary_plot <- ggplot() +
  geom_density(data=mut_rate_primary, size = 1, aes(color="lightcoral", x=mutation_rate_primary, y=..scaled..))+
  xlab("Mutation rate in primary tumors") + ylab("Density") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none",
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  theme(plot.margin = margin(0,0,0,26, "pt")) +
  scale_x_continuous(labels=scientific, limits=c(0, 1.25e-6), breaks=c(5e-7, 1e-6), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 1.05), breaks=c(0,0.2,0.4,0.6,0.8,1), expand = c(0, 0))
  
###Combine_Figures    
library(cowplot)

combined_gene_mutrate <- plot_grid(
  mutrates_early_plot, mutrates_early_metastasis_plot,
  mutrates_late_plot, mutrates_late_metastasis_plot,
  mutrates_metastasis_plot, mutrates_prim_met_plot,
  labels = c("A", "D", "B", "E", "C", "F"),
  label_size = 12,
  align = "h", axis = "tb",
  nrow = 3, ncol = 2,
  rel_heights = c(1, 1, 1)
)

ggsave("Figure_3.png", width = 8, dpi=600, height = 10)


#End

