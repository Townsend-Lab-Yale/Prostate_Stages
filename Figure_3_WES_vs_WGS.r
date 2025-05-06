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
Early_WES <- cesa$samples[Gleason == "Early" & coverage == "exome", unique(Unique_Patient_Identifier)]
Early_WGS <- cesa$samples[Gleason == "Early" & coverage == "genome", unique(Unique_Patient_Identifier)]

Late_WES <- cesa$samples[Gleason == "Late" & coverage == "exome", unique(Unique_Patient_Identifier)]
Late_WGS <- cesa$samples[Gleason == "Late" & coverage == "genome", unique(Unique_Patient_Identifier)]

#gene_mutation_rates_analysis:
cesa <- gene_mutation_rates(cesa, covariates = "PRAD", samples = Early_WES)
cesa <- gene_mutation_rates(cesa, covariates = "PRAD", samples = Early_WGS)
cesa <- gene_mutation_rates(cesa, covariates = "PRAD", samples = Late_WES)
cesa <- gene_mutation_rates(cesa, covariates = "PRAD", samples = Late_WGS)

#count the number of Syn and nonSyn in each datagroup:
dndscv_Early_WES <- cesa$dNdScv_results[[1]]
dndscv_Early_WGS <- cesa$dNdScv_results[[2]]
dndscv_Late_WES <- cesa$dNdScv_results[[3]]
dndscv_Late_WGS <- cesa$dNdScv_results[[4]]

#Save the files:
write.csv(dndscv_Early_WES, "dndscv_Early_WES.csv", row.names = FALSE)
write.csv(dndscv_Early_WGS, "dndscv_Early_WGS.csv", row.names = FALSE)
write.csv(dndscv_Late_WES, "dndscv_Late_WES.csv", row.names = FALSE)
write.csv(dndscv_Late_WGS, "dndscv_Late_WGS.csv", row.names = FALSE)


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

###Figure_early_WES_WGS:
threestage_final <- load_cesa("threestage_final.rds")

mut_rate_early_WES_WGS <- data.frame(gene=threestage_final@mutrates$gene,
                       WES_rate=threestage_final@mutrates$rate_grp_1,
                       WGS_rate=threestage_final@mutrates$rate_grp_2)

mut_rate_early_WES_WGS_plot <- ggplot()+
  geom_point(data = mut_rate_early_WES_WGS[!highlight,], aes(x=WES_rate, y=WGS_rate), size=1, color="lightcoral", alpha=0.5, shape = 16) +
  geom_point(data = mut_rate_early_WES_WGS[highlight,], aes(x=WES_rate, y=WGS_rate), size=2.5, color=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate_early_WES_WGS[highlight,], aes(x=WES_rate, y=WGS_rate, label=gene)) +
  theme_bw() +
  theme(panel.border = element_blank(), plot.title = element_text(hjust=0.5)) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  xlab("Mutation rate in lower-grade WES tumors") +
  ylab("Mutation rate in lower-grade WGS tumors") +
  geom_smooth(method="lm", color="navyblue") + 
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(0, 5e-7, 1e-6)) +
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6))

###Figure_late_WES_WGS:
mut_rate_late_WES_WGS <- data.frame(gene=threestage_final@mutrates$gene,
                                     WES_rate=threestage_final@mutrates$rate_grp_3,
                                     WGS_rate=threestage_final@mutrates$rate_grp_4)

mut_rate_late_WES_WGS_plot <- ggplot()+
  geom_point(data = mut_rate_late_WES_WGS[!highlight,], aes(x=WES_rate, y=WGS_rate), size=1, color="lightcoral", alpha=0.5, shape = 16) +
  geom_point(data = mut_rate_late_WES_WGS[highlight,], aes(x=WES_rate, y=WGS_rate), size=2.5, color=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate_late_WES_WGS[highlight,], aes(x=WES_rate, y=WGS_rate, label=gene)) +
  theme_bw() +
  theme(panel.border = element_blank(), plot.title = element_text(hjust=0.5)) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  xlab("Mutation rate in higher-grade WES tumors") +
  ylab("Mutation rate in higher-grade WGS tumors") +
  geom_smooth(method="lm", color="navyblue") + 
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(0, 5e-7, 1e-6)) +
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6))


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
Primary_WES <- Pri_Met$samples[Primary_Met == "Primary" & coverage == "exome", unique(Unique_Patient_Identifier)]
Primary_WGS <- Pri_Met$samples[Primary_Met == "Primary" & coverage == "genome", unique(Unique_Patient_Identifier)]


#gene_mutation_rates_analysis:
Pri_Met <- gene_mutation_rates(Pri_Met, covariates = "PRAD", samples = Primary_WES)
Pri_Met <- gene_mutation_rates(Pri_Met, covariates = "PRAD", samples = Primary_WGS)

#count the number of Syn and nonSyn in each datagroup:
dndscv_Primary_WES <- Pri_Met$dNdScv_results[[1]]
dndscv_Primary_WGS <- Pri_Met$dNdScv_results[[2]]

#Save the files:
write.csv(dndscv_Primary_WES, "dndscv_Primary_WES.csv", row.names = FALSE)
write.csv(dndscv_Primary_WGS, "dndscv_Primary_WGS.csv", row.names = FALSE)


twostage_final <- Pri_Met
saveRDS(twostage_final, file="twostage_final.rds")

###Figure_prim_WES_WGS: 
###mutrates_prim_met_plo
mut_rate_primary_WES_WGS <- data.frame(gene=twostage_final@mutrates$gene,
                                 WES_rate=twostage_final@mutrates$rate_grp_1,
                                 WGS_rate=twostage_final@mutrates$rate_grp_2)

mut_rate_primary_WES_WGS_plot <- ggplot()+
  geom_point(data = mut_rate_primary_WES_WGS[!highlight,], aes(x=WES_rate, y=WGS_rate), size=1, color="lightcoral", alpha=0.5, shape = 16) +
  geom_point(data = mut_rate_primary_WES_WGS[highlight,], aes(x=WES_rate, y=WGS_rate), size=2.5, colour=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate_primary_WES_WGS[highlight,], aes(x=WES_rate, y=WGS_rate, label=gene)) +
  theme_bw() +
  theme(panel.border = element_blank(), plot.title = element_text(hjust=0.5)) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.title = element_text(size=13),
        axis.text = element_text(size=10)) +
  xlab("Mutation rate in primary WES tumors") +
  ylab("Mutation rate in primary WGS tumors") +
  geom_smooth(method="lm", color="navyblue") +
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(0, 5e-7, 1e-6)) + 
  scale_y_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6))


###Combine_Figures    
library(cowplot)

mut_rate_early_WES_WGS_plot <- mut_rate_early_WES_WGS_plot +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 10))

mut_rate_late_WES_WGS_plot <- mut_rate_late_WES_WGS_plot +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 10))

mut_rate_primary_WES_WGS_plot <- mut_rate_primary_WES_WGS_plot +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 10))



shared_limits <- c(0, 1.25e-6)
shared_breaks <- c(0, 5e-7, 1e-6)

mut_rate_early_WES_WGS_plot <- mut_rate_early_WES_WGS_plot +
  scale_y_continuous(labels = scientific, limits = shared_limits, breaks = shared_breaks)

mut_rate_late_WES_WGS_plot <- mut_rate_late_WES_WGS_plot +
  scale_y_continuous(labels = scientific, limits = shared_limits, breaks = shared_breaks)

mut_rate_primary_WES_WGS_plot <- mut_rate_primary_WES_WGS_plot +
  scale_y_continuous(labels = scientific, limits = shared_limits, breaks = shared_breaks)


combined_gene_mutrate_WES__WGS <- plot_grid(
  mut_rate_early_WES_WGS_plot, 
  mut_rate_late_WES_WGS_plot, 
  mut_rate_primary_WES_WGS_plot,
  labels = c("A", "B", "C"), 
  label_size = 16,
  align = "v",      
  axis = "l",      
  nrow = 1, 
  ncol = 3, 
  rel_heights = c(1, 1, 1)
)

ggsave("combined_gene_mutrate_WES__WGS.png", width = 16, dpi=600, height = 5)


#End
