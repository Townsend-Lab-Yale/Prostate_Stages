
###Figure_4_SPOP_Primary

###Preparing data
MAF1 <- preload_maf(maf = "prad_armenia_final.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF2 <- preload_maf(maf = "prad_boutros_wgs_final.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF3 <- preload_maf(maf = "tcga_wgs_final.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF4 <- preload_maf(maf = "SU2C_PCF_dedup_filtered_11-09-20.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF5 <- preload_maf(maf = "MSK_Eur Urol_2020-341.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF6 <- preload_maf(maf = "MSK_Eur Urol_2020-410.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF7 <- preload_maf(maf = "MSK_Eur Urol_2020-468.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")

### keep samples where the column "Gleason" is equal to "Early" & "Late".  
MAF1 <- MAF1[MAF1$Gleason == "Early" | MAF1$Gleason == "Late", ]
MAF2 <- MAF2[MAF2$Gleason == "Early" | MAF2$Gleason == "Late", ]
MAF3 <- MAF3[MAF3$Gleason == "Early" | MAF3$Gleason == "Late", ]
MAF4 <- MAF4[MAF4$Gleason == "Early" | MAF4$Gleason == "Late", ]
MAF5 <- MAF5[MAF5$Gleason == "Early" | MAF5$Gleason == "Late", ]
MAF6 <- MAF6[MAF6$Gleason == "Early" | MAF6$Gleason == "Late", ]
MAF7 <- MAF7[MAF7$Gleason == "Early" | MAF7$Gleason == "Late", ]

# In four of these files , MAF4, MAF5, MAF6 and MAF5, the column "Gleason" is equal to either "Empty" or "Metastasis" (there are no information about Early and Late groups), so they would be removed from the rest of analysis.

# removing samples where column Problem is equal to NA
MAF1 <- MAF1[is.na(problem)]
MAF2 <- MAF2[is.na(problem)]
MAF3 <- MAF3[is.na(problem)]

# keeping only samples that do not occur at germline variant sites
MAF1 <- MAF1[germline_variant_site == F]
MAF2 <- MAF2[germline_variant_site == F]
MAF3 <- MAF3[germline_variant_site == F]


# keeping only samples that do not occur in repetitive regions 
MAF1 <- MAF1[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF2 <- MAF2[(repetitive_region == F | cosmic_site_tier %in% 1:3)]
MAF3 <- MAF3[(repetitive_region == F | cosmic_site_tier %in% 1:3)]


SPOP1 = CESAnalysis("ces.refset.hg19")
saveRDS(SPOP1, file = "SPOP1.rds")

SPOP2 = load_maf(cesa = SPOP1, maf = MAF1)
SPOP2 = load_maf(cesa = SPOP2, maf = MAF2, coverage = "genome")
SPOP2 = load_maf(cesa = SPOP2, maf = MAF3, coverage = "genome")

saveRDS(SPOP2, file = "SPOP2.rds")

to_remove = suggest_cosmic_signatures_to_remove(cancer_type = "PRAD")

SPOP3 = trinuc_mutation_rates(SPOP2, signature_set = "COSMIC_v3.2",
                              signature_extractor = "deconstructSigs",
                              signatures_to_remove = to_remove)

saveRDS(SPOP3, file = "SPOP3.rds")

SPOP4 = gene_mutation_rates(SPOP3, covariates = "PRAD")
saveRDS(SPOP4, file = "SPOP4.rds")

SPOP_final = ces_variant(SPOP4, variants = select_variants(SPOP4, min_freq = 2), model = "sswm")
saveRDS(SPOP_final, file="SPOP_final.rds")

### making the figure:

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

common.text.size <- 4

PRAD_analysis <- load_cesa("SPOP_final.rds")
PRAD_results <- snv_results(PRAD_analysis)
PRAD_results <- PRAD_results$selection.1

# extract  and add the gene names:
gene_name <- word(PRAD_results$variant_name, sep = "_")
PRAD_results$gene <- gene_name

# remove "SPOP_" from variant_name:
PRAD_results$variant_name <- str_replace(PRAD_results$variant_name, "SPOP_", "")

# keep only the rows where the variant_type == "aac".
aac <- PRAD_results$variant_type == "aac"
PRAD_results <- PRAD_results[aac,]

PRAD_results_recurrent <- PRAD_results[order(-selection_intensity),]
spop_true <- PRAD_results_recurrent$gene == "SPOP"
PRAD_results_recurrent <- PRAD_results_recurrent[spop_true,]

#########################################################################

bargraph_SPOP_SI <- ggplot(data=PRAD_results_recurrent, aes(x=reorder(variant_name, -selection_intensity), y=selection_intensity, fill=selection_intensity))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=ci_low_95, ymax=ci_high_95), width=0.333) +
  theme(axis.text.x = element_text (hjust = 1, angle = 45)) +
  xlab("SPOP amino acid substitution (and prevalence)") + ylab("Scaled selection coefficient") +
  scale_fill_gradient(low="gold", high="red2") +
  theme(legend.position = "none")+
  theme(panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  geom_text(aes(label=included_with_variant, y=-5000), size = common.text.size) +
  scale_y_continuous(labels=scientific, breaks = c(0, 5e4, 1e5, 1.5e5, 2e5))

ggsave("SPOP_recurrent_SI.png", width=8, height=5.25)

#End