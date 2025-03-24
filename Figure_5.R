library(ggplot2)
library(stringr)
library(scales)
library(gg.gap)
library(cowplot)
library(reshape2)
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


consensus_gr = Reduce(GenomicRanges::intersect, unlist(cesa$coverage_ranges[c('exome', 'targeted')]))



selected_genes <- c("SPOP", "FOXA1", "AR", "PIK3CA", "PIK3CB", "TP53", "ROCK1", "RHOA", "AKT1", "ATM", "CUL3",
                    "APC", "CTNNB1", "PTEN", "KMT2C", "KMT2D")

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
cesa <- trinuc_mutation_rates(cesa = cesa, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)

# calculate gene rates for all samples
cesa <- gene_mutation_rates(cesa, covariates = "PRAD", save_all_dndscv_output = T)

cesa <- ces_gene_epistasis(cesa = cesa, genes = selected_genes, variants = "recurrent", run_name = "gene_epistasis_example")

epistasiiiiiis <- cesa$epistasis$gene_epistasis_example

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

###SPOP###

#Select your gene pairs of interest
SPOP_list <- which(epistasis_results$variant_A == "SPOP" | epistasis_results$variant_B == "SPOP")

epistatic_change_SPOP <- c()

#Decoupling the gene pairs
for(x in SPOP_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(epistasis_results[x,5] - epistasis_results[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(epistasis_results[x,6] - epistasis_results[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(epistasis_results[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(epistasis_results[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(epistasis_results[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(epistasis_results[x,2]))
  epistatic_change_SPOP <- rbind(epistatic_change_SPOP, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_SPOP <- data.frame(gene = epistatic_change_SPOP[,1], change = as.numeric(epistatic_change_SPOP[,2]))

#Separating "Before" and "After"
epistatic_change_SPOP_before <- epistatic_change_SPOP[grep("SPOP_", epistatic_change_SPOP[,1]),]
epistatic_change_SPOP_before$time <- rep("Before", 14)
epistatic_change_SPOP_before <- epistatic_change_SPOP_before[order(epistatic_change_SPOP_before$change),]
epistatic_change_SPOP_after <- epistatic_change_SPOP[grep("_SPOP", epistatic_change_SPOP[,1]),]
epistatic_change_SPOP_after$time <- rep("After", 14)
epistatic_change_SPOP_after <- epistatic_change_SPOP_after[order(epistatic_change_SPOP_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_SPOP_before[,1] <- sub("SPOP_after_", "", epistatic_change_SPOP_before[,1])
epistatic_change_SPOP_after[,1] <- sub("after_SPOP", "", epistatic_change_SPOP_after[,1])


#Blank spot
blank <- data.frame(gene = "BLANK", change = -5000, time = "Before")

epistatic_change_SPOP <- rbind(epistatic_change_SPOP_before, epistatic_change_SPOP_after)

epistatic_change_SPOP$gene <- factor(epistatic_change_SPOP$gene,
                                     levels = c("KMT2C", "AKT1", "CUL3", "PIK3CB", "APC", "ATM", "KMT2D", "PTEN",
                                                "FOXA1", "TP53", "CTNNB1", "AR", "PIK3CA", "RHOA",
                                                "CUL3_", "PIK3CB_", "PIK3CA_", "AR_", "ATM_",
                                                "APC_", "TP53_", "KMT2D_", "CTNNB1_", "KMT2C_", "PTEN_", "FOXA1_",
                                                "AKT1_", "RHOA_"))

gene_labels_SPOP <- c(epistatic_change_SPOP_before$gene, "", epistatic_change_SPOP_after$gene)
gene_labels_SPOP <- sub("_", "", gene_labels_SPOP)


waterfall_SPOP <- ggplot(epistatic_change_SPOP, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.line.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("SPOP gene pairs") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels_SPOP) + scale_fill_discrete(breaks = c("Before", "After"))+
  scale_y_continuous(labels = scientific) +
  geom_hline(yintercept = 0)

waterfall_SPOP
ggsave("PRAD_figures/epistasis_16/waterfall_SPOP.png", width = 10, dpi=300, height = 7)

###PIK3CA###

#Select your gene pairs of interest
PIK3CA_list <- which(epistasis_results$variant_A == "PIK3CA" | epistasis_results$variant_B == "PIK3CA")

epistatic_change_PIK3CA <- c()

#Decoupling the gene pairs
for(x in PIK3CA_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(epistasis_results[x,5] - epistasis_results[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(epistasis_results[x,6] - epistasis_results[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(epistasis_results[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(epistasis_results[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(epistasis_results[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(epistasis_results[x,2]))
  epistatic_change_PIK3CA <- rbind(epistatic_change_PIK3CA, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_PIK3CA <- data.frame(gene = epistatic_change_PIK3CA[,1], change = as.numeric(epistatic_change_PIK3CA[,2]))

#Separating "Before" and "After"
epistatic_change_PIK3CA_before <- epistatic_change_PIK3CA[grep("PIK3CA_", epistatic_change_PIK3CA[,1]),]
epistatic_change_PIK3CA_before$time <- rep("Before", 14)
epistatic_change_PIK3CA_before <- epistatic_change_PIK3CA_before[order(epistatic_change_PIK3CA_before$change),]
epistatic_change_PIK3CA_after <- epistatic_change_PIK3CA[grep("_PIK3CA", epistatic_change_PIK3CA[,1]),]
epistatic_change_PIK3CA_after$time <- rep("After", 14)
epistatic_change_PIK3CA_after <- epistatic_change_PIK3CA_after[order(epistatic_change_PIK3CA_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_PIK3CA_before[,1] <- sub("PIK3CA_after_", "", epistatic_change_PIK3CA_before[,1])
epistatic_change_PIK3CA_after[,1] <- sub("after_PIK3CA", "", epistatic_change_PIK3CA_after[,1])

#Blank spot

epistatic_change_PIK3CA <- rbind(epistatic_change_PIK3CA_before, epistatic_change_PIK3CA_after)


epistatic_change_PIK3CA$gene <- factor(epistatic_change_PIK3CA$gene,
                                       levels = c("SPOP", "AR", "CTNNB1", "KMT2C", "PIK3CB", "AKT1", "CUL3", "APC",
                                                  "KMT2D", "ATM", "RHOA", "FOXA1", "TP53", "PTEN",
                                                  "CUL3_", "RHOA_", "AKT1_", "PIK3CB_", "ATM_",
                                                  "APC_", "TP53_", "PTEN_", "KMT2D_", "KMT2C_", "AR_", "FOXA1_",
                                                  "SPOP_", "CTNNB1_"))

gene_labels_PIK3CA <- c(epistatic_change_PIK3CA_before$gene, "", epistatic_change_PIK3CA_after$gene)
gene_labels_PIK3CA <- sub("_", "", gene_labels_PIK3CA)

waterfall_PIK3CA <- ggplot(epistatic_change_PIK3CA, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.line.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("PIK3CA gene pairs") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels_PIK3CA) + scale_fill_discrete(breaks = c("Before", "After"))+
  scale_y_continuous(labels = scientific) +
  geom_hline(yintercept = 0)

waterfall_PIK3CA
ggsave("PRAD_figures/epistasis_16/waterfall_PIK3CA.png", width = 10, dpi=300, height = 7)



###TP53###

#Select your gene pairs of interest
TP53_list <- which(epistasiiiiiis$variant_A == "TP53" | epistasiiiiiis$variant_B == "TP53")

epistatic_change_TP53 <- c()

#Decoupling the gene pairs
for(x in TP53_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(epistasiiiiiis[x,5] - epistasiiiiiis[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(epistasiiiiiis[x,6] - epistasiiiiiis[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(epistasiiiiiis[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(epistasiiiiiis[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(epistasiiiiiis[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(epistasiiiiiis[x,2]))
  epistatic_change_TP53 <- rbind(epistatic_change_TP53, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_TP53 <- data.frame(gene = epistatic_change_TP53[,1], change = as.numeric(epistatic_change_TP53[,2]))

#Separating "Before" and "After"
epistatic_change_TP53_before <- epistatic_change_TP53[grep("TP53_", epistatic_change_TP53[,1]),]
epistatic_change_TP53_before$time <- rep("Before", 15)
epistatic_change_TP53_before <- epistatic_change_TP53_before[order(epistatic_change_TP53_before$change),]
epistatic_change_TP53_after <- epistatic_change_TP53[grep("_TP53", epistatic_change_TP53[,1]),]
epistatic_change_TP53_after$time <- rep("After", 15)
epistatic_change_TP53_after <- epistatic_change_TP53_after[order(epistatic_change_TP53_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_TP53_before[,1] <- sub("TP53_after_", "", epistatic_change_TP53_before[,1])
epistatic_change_TP53_after[,1] <- sub("after_TP53", "", epistatic_change_TP53_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = -5000, time = "Before")

epistatic_change_TP53 <- rbind(epistatic_change_TP53_before, blank, epistatic_change_TP53_after)

epistatic_change_TP53$gene <- factor(epistatic_change_TP53$gene,
                                     levels = c("AKT1", "KMT2D", "RHOA", "AR", "PTEN", "ROCK1", "FOXA1", "SPOP",
                                                "ATM", "CTNNB1", "PIK3CA", "CUL3", "PIK3CB", "APC", "KMT2C", "BLANK",
                                                "CUL3_", "ROCK1_", "SPOP_", "RHOA_", "CTNNB1_", "PIK3CA_", "AKT1_",
                                                "APC_", "FOXA1_", "KMT2D_", "PIK3CB_", "KMT2C_", "ATM_", "AR_", "PTEN_"))


gene_labels_TP53 <- c(epistatic_change_TP53_before$gene, "", epistatic_change_TP53_after$gene)
gene_labels_TP53 <- sub("_", "", gene_labels_TP53)

waterfall_TP53 <- ggplot(epistatic_change_TP53, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.line.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("TP53 gene pairs") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels_TP53) + scale_fill_discrete(breaks = c("Before", "After"))+
  scale_y_continuous(labels = scientific, breaks = c(-2e4,-1.5e4, -1e4, 0, 1e4), limits = c(-2.5e4,1e4)) +
  geom_hline(yintercept = 0)

waterfall_TP53
ggsave("PRAD_figures/epistasis_16/waterfall_TP53.png", width = 10, dpi=300, height = 7)


###AR###

#Select your gene pairs of interest
AR_list <- which(epistasiiiiiis$variant_A == "AR" | epistasiiiiiis$variant_B == "AR")

epistatic_change_AR <- c()

#Decoupling the gene pairs
for(x in AR_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(epistasiiiiiis[x,5] - epistasiiiiiis[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(epistasiiiiiis[x,6] - epistasiiiiiis[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(epistasiiiiiis[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(epistasiiiiiis[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(epistasiiiiiis[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(epistasiiiiiis[x,2]))
  epistatic_change_AR <- rbind(epistatic_change_AR, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_AR <- data.frame(gene = epistatic_change_AR[,1], change = as.numeric(epistatic_change_AR[,2]))

#Separating "Before" and "After"
epistatic_change_AR_before <- epistatic_change_AR[grep("AR_", epistatic_change_AR[,1]),]
epistatic_change_AR_before$time <- rep("Before", 15)
epistatic_change_AR_before <- epistatic_change_AR_before[order(-epistatic_change_AR_before$change),]
epistatic_change_AR_after <- epistatic_change_AR[grep("_AR", epistatic_change_AR[,1]),]
epistatic_change_AR_after$time <- rep("After", 15)
epistatic_change_AR_after <- epistatic_change_AR_after[order(-epistatic_change_AR_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_AR_before[,1] <- sub("AR_after_", "", epistatic_change_AR_before[,1])
epistatic_change_AR_after[,1] <- sub("after_AR", "", epistatic_change_AR_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = 0, time = "Before")

epistatic_change_AR <- rbind(epistatic_change_AR_before, blank, epistatic_change_AR_after)

epistatic_change_AR$gene <- factor(epistatic_change_AR$gene,
                                   levels = c("PTEN", "ROCK1", "AKT1", "KMT2C", "CUL3", "PIK3CB", "APC", "RHOA", "TP53",
                                              "FOXA1", "CTNNB1", "PIK3CA", "SPOP", "ATM", "KMT2D", "BLANK", "CUL3_",
                                              "ROCK1_", "SPOP_", "RHOA_", "AKT1_", "PIK3CB_", "PIK3CA_", "APC_",
                                              "TP53_", "FOXA1_", "CTNNB1_", "KMT2D_", "KMT2C_", "PTEN_", "ATM_"))

gene_labels_AR <- sub("_", "", c("PTEN", "ROCK1", "AKT1", "KMT2C", "CUL3", "PIK3CB", "APC", "RHOA", "TP53",
                                 "FOXA1", "CTNNB1", "PIK3CA", "SPOP", "ATM", "KMT2D", "BLANK", "CUL3_",
                                 "ROCK1_", "SPOP_", "RHOA_", "AKT1_", "PIK3CB_", "PIK3CA_", "APC_",
                                 "TP53_", "FOXA1_", "CTNNB1_", "KMT2D_", "KMT2C_", "PTEN_", "ATM_"))

waterfall_AR <- ggplot(epistatic_change_AR, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.line.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("AR gene pairs") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels_AR) +
  scale_fill_discrete(breaks = c("Before", "After"))+
  scale_y_continuous(labels = scientific, limits = c(-3e4, 3e4), breaks = c(-2e4, -1.5e4, 0, 2e4, 2.8e4)) +
  geom_hline(yintercept = 0)

waterfall_AR

ggsave("PRAD_figures/epistasis_16/waterfall_AR.png", width = 10, dpi=300, height = 7)

###############################################

combined_waterfall <- plot_grid(waterfall_SPOP, waterfall_PIK3CA, waterfall_TP53, waterfall_AR,
                                labels = c("A", "B", "C", "D"), ncol = 2)

combined_waterfall

ggsave("PRAD_figures/epistasis_16/combined_waterfall.png", width = 10, dpi=300, height = 7.5)









# Clear gene rates and calculate gene rates for all samples (not separated by normal and tumor) for epistasis ----
cesa <- clear_gene_rates(cesa)
cesa <- gene_mutation_rates(cesa, covariates = "ESCA", save_all_dndscv_output = T)

dndscv_gene_names <- cesa$gene_rates$gene
nsyn_sites <- sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

samples_in_all <- length(unique(cesa$dNdScv_results$rate_grp_1$annotmuts$sampleID ))

mut_rate_df <- tibble(gene = cesa$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_mu = cesa$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df <- mut_rate_df %>% 
  mutate(total_mu = (exp_mu / n_syn_sites) / samples_in_all) %>%
  select(gene, total_mu) %>%
  data.table::setDT()

cesa <- clear_gene_rates(cesa = cesa)
cesa <- set_gene_rates(cesa = cesa, rates = mut_rate_df, missing_genes_take_nearest = T) 


cesa <- ces_epistasis(cesa, variants = compound, run_name = "epistasis_compound_variants_all_samples")


save_cesa(cesa = cesa, file = "analysis/eso_cesa_after_analysis.rds")
