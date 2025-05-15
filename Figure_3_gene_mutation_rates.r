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

##Mutational processes and relative mutation rates:
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

cesa = trinuc_mutation_rates(cesa, ces.refset.hg19$signatures$COSMIC_v3.2,
                                  signature_exclusions = signature_exclusions)

library(ggplot2)
library(scales)

#defining groups:
Early_groups <- cesa$samples[Gleason == "Early", unique(Unique_Patient_Identifier)]
Late_groups <- cesa$samples[Gleason == "Late", unique(Unique_Patient_Identifier)]
Metastasis_groups <- cesa$samples[Gleason == "Metastasis", unique(Unique_Patient_Identifier)]

#gene_mutation_rates_analysis:
cesa <- gene_mutation_rates(cesa, covariates = "PRAD", samples = Early_groups)
cesa <- gene_mutation_rates(cesa, covariates = "PRAD", samples = Late_groups)
cesa <- gene_mutation_rates(cesa, covariates = "PRAD", samples = Metastasis_groups)

threestage_final <- cesa
saveRDS(threestage_final, file="threestage_final.rds")


#groups:
mut_rate_stageless <- data.frame(gene=cesa@mutrates$gene,
                                 mutation_rate=cesa@mutrates$rate_grp_1)

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
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(0, 5e-7, 1e-6)) +
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6))

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
  ylab("Mutation rate in metastases") +
  geom_smooth(method="lm", color="navyblue") + 
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(0, 5e-7, 1e-6)) +
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6))

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
  ylab("Mutation rate in metastases") +
  geom_smooth(method="lm", color="navyblue") + 
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(0, 5e-7, 1e-6)) +
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6))

### Figure_separate_tumors:
##Fig_mut_rate_early
mut_rate_early <- data.frame(gene=cesa@mutrates$gene,
                                 mutation_rate_early=cesa@mutrates$rate_grp_1)


mutrates_early_plot <- ggplot() +
  geom_density(data=mut_rate_early, size = 1, aes(color="lightcoral", x=mutation_rate_early, y=..scaled..))+
  xlab("Mutation rate in low-grade tumors") + ylab("Density") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none",
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  theme(plot.margin = margin(0,0,0,26, "pt")) +
  scale_x_continuous(labels=scientific, limits=c(0, 1.25e-6), breaks=c(0, 5e-7, 1e-6))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05), breaks=c(0,0.5,1))

##Fig_mut_rate_late
mut_rate_late <- data.frame(gene=cesa@mutrates$gene,
                                 mutation_rate_late=cesa@mutrates$rate_grp_2)


mutrates_late_plot <- ggplot() +
  geom_density(data=mut_rate_late, size = 1, aes(color="lightcoral", x=mutation_rate_late, y=..scaled..))+
  xlab("Mutation rate in high-grade tumors") + ylab("Density") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none",
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  theme(plot.margin = margin(0,0,0,26, "pt")) +
  scale_x_continuous(labels=scientific, limits=c(0, 1.25e-6), breaks=c(0, 5e-7, 1e-6))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05), breaks=c(0,0.5,1))
 
##Fig_mut_rate_metastasis
mut_rate_metastasis <- data.frame(gene=cesa@mutrates$gene,
                                 mutation_rate_metastasis=cesa@mutrates$rate_grp_3)


mutrates_metastasis_plot <- ggplot() +
  geom_density(data=mut_rate_metastasis, size = 1, aes(color="lightcoral", x=mutation_rate_metastasis, y=..scaled..))+
  xlab("Mutation rate in metastases") + ylab("Density") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none",
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  theme(plot.margin = margin(0,0,0,26, "pt")) +
  scale_x_continuous(labels=scientific, limits=c(0, 1.25e-6), breaks=c(0, 5e-7, 1e-6))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05), breaks=c(0,0.5,1))
 
##Fig_mut_rate_all
mut_rate_all <- clear_gene_rates(cesa)

mut_rate_all <- gene_mutation_rates(mut_rate_all, covariates = "PRAD")

mut_rate_all_stageless <- data.frame(gene=mut_rate_all@mutrates$gene,
                                 mutation_rate_all_stageless=mut_rate_all@mutrates$rate_grp_1)


mutrates_all_plot <- ggplot() +
  geom_density(data=mut_rate_all_stageless, size = 1, aes(color="lightcoral", x=mutation_rate_all_stageless, y=..scaled..))+
  xlab("Mutation rate in all tumors") + ylab("Density") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none",
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  theme(plot.margin = margin(0,0,0,26, "pt")) +
  scale_x_continuous(labels=scientific, limits=c(0, 1.25e-6), breaks=c(0, 5e-7, 1e-6))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05), breaks=c(0,0.5,1))

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

##Mutational processes and relative mutation rates:
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

Pri_Met = trinuc_mutation_rates(Pri_Met, ces.refset.hg19$signatures$COSMIC_v3.2,
                                  signature_exclusions = signature_exclusions)

#defining groups:
Primary_groups <- Pri_Met$samples[Primary_Met == "Primary", unique(Unique_Patient_Identifier)]
Metastasis_groups <- Pri_Met$samples[Primary_Met == "Metastasis", unique(Unique_Patient_Identifier)]

#gene_mutation_rates_analysis:
Pri_Met <- gene_mutation_rates(Pri_Met, covariates = "PRAD", samples = Primary_groups)
Pri_Met <- gene_mutation_rates(Pri_Met, covariates = "PRAD", samples = Metastasis_groups)

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
  ylab("Mutation rate in metastases") +
  geom_smooth(method="lm", color="navyblue") +
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(0, 5e-7, 1e-6)) + 
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6))

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
  scale_x_continuous(labels=scientific, limits=c(0, 1.25e-6), breaks=c(0, 5e-7, 1e-6))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.05), breaks=c(0,0.5,1))
  
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

ggsave("combined_gene_mutrate.png", width = 8, dpi=600, height = 10)


#End
