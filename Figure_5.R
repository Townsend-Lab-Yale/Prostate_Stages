library(cancereffectsizeR)
library(MutationalPatterns)
library(ces.refset.hg19)
library(data.table)
library(scales)
library(stringr)
library(ggplot2)

###Figure_5_AR_Metastasis

###Preparing data
MAF1 <- preload_maf(maf = "prad_armenia_final.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF2 <- preload_maf(maf = "prad_boutros_wgs_final.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF3 <- preload_maf(maf = "tcga_wgs_final.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF4 <- preload_maf(maf = "SU2C_PCF_dedup_filtered_11-09-20.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF5 <- preload_maf(maf = "MSK_Eur Urol_2020-341.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF6 <- preload_maf(maf = "MSK_Eur Urol_2020-410.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF7 <- preload_maf(maf = "MSK_Eur Urol_2020-468.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")

### keep samples where the column "Gleason" is equal to "Metastasis" 
MAF1 <- MAF1[MAF1$Gleason == "Metastasis", ]
MAF2 <- MAF2[MAF2$Gleason == "Metastasis", ]
MAF3 <- MAF3[MAF3$Gleason == "Metastasis", ]
MAF4 <- MAF4[MAF4$Gleason == "Metastasis", ]
MAF5 <- MAF5[MAF5$Gleason == "Metastasis", ]
MAF6 <- MAF6[MAF6$Gleason == "Metastasis", ]
MAF7 <- MAF7[MAF7$Gleason == "Metastasis", ]

# In two of these files, MAF2 and MAF3, the column "Gleason" is equal to "Empty" ("Early" or "Late") (there are no Metastasis), so they would be removed from the rest of analysis.

# removing samples where column Problem is equal to NA
MAF1 <- MAF1[is.na(problem)]
MAF4 <- MAF4[is.na(problem)]
MAF5 <- MAF5[is.na(problem)]
MAF6 <- MAF6[is.na(problem)]
MAF7 <- MAF7[is.na(problem)]

# keeping only samples that do not occur at germline variant sites
MAF1 <- MAF1[germline_variant_site == F]
MAF4 <- MAF4[germline_variant_site == F]
MAF5 <- MAF5[germline_variant_site == F]
MAF6 <- MAF6[germline_variant_site == F]
MAF7 <- MAF7[germline_variant_site == F]

# keeping only samples that do not occur in repetitive regions 
MAF1 <- MAF1[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF4 <- MAF4[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF5 <- MAF5[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF6 <- MAF6[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF7 <- MAF7[(repetitive_region == F | cosmic_site_tier %in% 1:3)]

AR1 = CESAnalysis("ces.refset.hg19")
saveRDS(AR1, file = "AR1.rds")

AR2 = load_maf(cesa = AR1, maf = MAF1)
AR2 = load_maf(cesa = AR2, maf = MAF4, coverage = "exome",
               covered_regions = "SureSelect_All_Exon_covered_regions.bed", 
               covered_regions_name = "SureSelect_V4", covered_regions_padding = 100)
AR2 = load_maf(cesa = AR2, maf = MAF5, coverage = "exome",
               covered_regions = "msk_341_exons.bed", 
               covered_regions_name = "MSK_IMPACT_341", covered_regions_padding = 10)
AR2 = load_maf(cesa = AR2, maf = MAF6, coverage = "exome", 
               covered_regions = "msk_410_exons.bed",
               covered_regions_name = "MSK_IMPACT_410", covered_regions_padding = 10)
AR2 = load_maf(cesa = AR2, maf = MAF7, coverage = "targeted",
               covered_regions = "msk_468_exons.bed",
               covered_regions_name = "MSK_IMPACT_468", covered_regions_padding = 10)

saveRDS(AR2, file = "AR2.rds")

to_remove = suggest_cosmic_signatures_to_remove(cancer_type = "PRAD")

AR3 = trinuc_mutation_rates(AR2, signature_set = "COSMIC_v3.2",
                            signature_extractor = "deconstructSigs",
                            signatures_to_remove = to_remove)
saveRDS(AR3, file = "AR3.rds")

AR4 = gene_mutation_rates(AR3, covariates = "PRAD")
saveRDS(AR4, file = "AR4.rds")

AR_final = ces_variant(AR4, variants = select_variants(AR4, min_freq = 2), model = "sswm")
saveRDS(AR_final, file="AR_final.rds")

### making the figure:

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

common.text.size <- 4

PRAD_analysis <- load_cesa("AR_final.rds")
PRAD_results <- snv_results(PRAD_analysis)
PRAD_results <- PRAD_results$selection.1

# extract and add the gene names:
gene_name <- word(PRAD_results$variant_name, sep = "_")
PRAD_results$gene <- gene_name

# remove "AR_" from variant_name:
PRAD_results$variant_name <- str_replace(PRAD_results$variant_name, "AR_", "")

# keep only the rows where the variant_type == "aac".
aac <- PRAD_results$variant_type == "aac"
PRAD_results <- PRAD_results[aac,]

PRAD_results_recurrent <- PRAD_results[order(-selection_intensity),]
AR_true <- PRAD_results_recurrent$gene == "AR"
PRAD_results_recurrent <- PRAD_results_recurrent[AR_true,]

#########################################################################

bargraph_AR_SI <- ggplot(data=PRAD_results_recurrent, aes(x=reorder(variant_name, -selection_intensity), y=selection_intensity, fill=selection_intensity))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=ci_low_95, ymax=ci_high_95), width=0.333) +
  theme(axis.text.x = element_text (hjust = 1, angle = 45)) +
  xlab("AR amino acid substitution (and prevalence)") + ylab("Scaled selection coefficient") +
  scale_fill_gradient(low="gold", high="red2") +
  theme(legend.position = "none")+
  theme(panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  geom_text(aes(label=included_with_variant, y=-5000), size = common.text.size) +
  scale_y_continuous(labels=scientific, breaks = c(0, 1e4, 2e4, 3e4, 4e4, 5e4, 1e5, 1.5e5, 2e5))

ggsave("AR_recurrent_SI.png", width=8, height=5.25)

#End
