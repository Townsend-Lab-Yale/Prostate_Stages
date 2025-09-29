library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg19)
library(MutationalPatterns)
library(RColorBrewer)
library(ggrepel)

setwd("C:/Moein/projects/prostate_stages/Age_comparison_3bp_motifs")


###Preparing data

gleason <- read.delim("C:/Moein/projects/prostate_stages/Age_comparison_3bp_motifs/gleason_age_comparison.txt")

MAF1 <- preload_maf(maf = "armenia_final.maf.txt", refset = ces.refset.hg19)


# removing samples where column Problem is equal to NA
MAF1 <- MAF1[is.na(problem)]

# keeping only samples that do not occur at germline variant sites
MAF1 <- MAF1[germline_variant_site == F]


# keeping only samples that do not occur in repetitive regions 
MAF1 <- MAF1[(repetitive_region == F | cosmic_site_tier %in% 1:3)]


##Create CESAnalysis and load data:
cesa <- CESAnalysis(refset = "ces.refset.hg19")
cesa <- load_maf(cesa = cesa, maf = MAF1, maf_name = "armenia")


cesa <- load_sample_data(cesa, gleason)

##Mutational processes and relative mutation rates:
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

cesa = trinuc_mutation_rates(cesa, ces.refset.hg19$signatures$COSMIC_v3.2,
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
colnames(summed_snv_by_group)[c(1, 2)] <- c("High-grade", "Low-grade")
summed_snv_by_group <- summed_snv_by_group[, c("Low-grade", "High-grade")]
rownames(summed_snv_by_group) <- rownames(snv_counts)
Figure_1 <- MutationalPatterns::plot_96_profile(summed_snv_by_group, ymax = 0.15)
ggsave("Figure_Age_comparison.png", width = 8, height = 6, dpi = 600)

#End

