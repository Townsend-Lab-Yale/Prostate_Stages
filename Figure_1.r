library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg19)
library(MutationalPatterns)
library(RColorBrewer)
library(ggrepel)


##Preparing data:

setwd("C:/Moein/projects/prostate_stages/PRAD_files")

gleason <- read.delim("C:/Moein/projects/prostate_stages/PRAD_files/gleason.txt")

MAF1 <- preload_maf(maf = "armenia_final.maf.txt", refset = ces.refset.hg19)
MAF2 <- preload_maf(maf = "boutros_final.maf.txt", refset = ces.refset.hg19)
MAF3 <- preload_maf(maf = "tcga_final.maf.txt", refset = ces.refset.hg19)
MAF4 <- preload_maf(maf = "SU2C_final.maf.txt", refset = ces.refset.hg19)
MAF5 <- preload_maf(maf = "MSK_341_final.maf.txt", refset = ces.refset.hg19)
MAF6 <- preload_maf(maf = "MSK_410_final.maf.txt", refset = ces.refset.hg19)
MAF7 <- preload_maf(maf = "MSK_468_final.maf.txt", refset = ces.refset.hg19)


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
                                  signature_extractor = "deconstructSigs",
                                  signature_exclusions = signature_exclusions)

##Figure_1:

snv_counts <- cesa$mutational_signatures$snv_counts

summed_snv_by_group <- data.table()
receptor_groups <- unique(na.omit(cesa$samples$Gleason))
samples_with_snvs <- cesa$samples[colnames(snv_counts), on = "Unique_Patient_Identifier"]
for (grp in receptor_groups) {
  curr_samples <- samples_with_snvs[grp, Unique_Patient_Identifier, on = "Gleason"]
  curr_snv_sum <- rowSums(snv_counts[, curr_samples])
  summed_snv_by_group[, (grp) := curr_snv_sum]
}
summed_snv_by_group <- as.matrix(summed_snv_by_group)
colnames(summed_snv_by_group)[c(1, 2, 3)] <- c("Metastases", "Higher-risk", "Lower-risk")
summed_snv_by_group <- summed_snv_by_group[, c("Lower-risk", "Higher-risk", "Metastases")]
rownames(summed_snv_by_group) <- rownames(snv_counts)
Figure_1 <- MutationalPatterns::plot_96_profile(summed_snv_by_group, ymax = 0.15)
ggsave("Figure_1.png", width = 8, height = 6, dpi = 600)

#End

