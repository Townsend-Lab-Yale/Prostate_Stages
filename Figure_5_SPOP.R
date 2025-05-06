library(cancereffectsizeR)
library(MutationalPatterns)
library(ces.refset.hg19)
library(data.table)
library(scales)
library(stringr)
library(ggplot2)

setwd("C:/Moein/projects/prostate_stages/PRAD_files")
###Figure_4_SPOP

###Preparing data
MAF1 <- preload_maf(maf = "prad_armenia_final.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF2 <- preload_maf(maf = "prad_boutros_wgs_final.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF3 <- preload_maf(maf = "tcga_wgs_final.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF4 <- preload_maf(maf = "SU2C_PCF_dedup_filtered_11-09-20.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF5 <- preload_maf(maf = "MSK_Eur Urol_2020-341.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF6 <- preload_maf(maf = "MSK_Eur Urol_2020-410.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")
MAF7 <- preload_maf(maf = "MSK_Eur Urol_2020-468.maf.txt", refset = ces.refset.hg19, keep_extra_columns = "Gleason")

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

# keeping snv:
MAF1 <- subset(MAF1, variant_type == "snv")
MAF2 <- subset(MAF2, variant_type == "snv")
MAF3 <- subset(MAF3, variant_type == "snv")
MAF4 <- subset(MAF4, variant_type == "snv")
MAF5 <- subset(MAF5, variant_type == "snv")
MAF6 <- subset(MAF6, variant_type == "snv")
MAF7 <- subset(MAF7, variant_type == "snv")

SPOP1 = CESAnalysis("ces.refset.hg19")
saveRDS(SPOP1, file = "SPOP1.rds")

SPOP2 = load_maf(cesa = SPOP1, maf = MAF1)
SPOP2 = load_maf(cesa = SPOP2, maf = MAF2, coverage = "genome")
SPOP2 = load_maf(cesa = SPOP2, maf = MAF3, coverage = "genome")
SPOP2 <- load_maf(cesa = SPOP2, maf = MAF4, maf_name = "SU2C",, coverage = "exome",
                 covered_regions = "SureSelect_All_Exon_covered_regions.bed",
                 covered_regions_name = "SureSelect_V4", covered_regions_padding = 100)
SPOP2 <- load_maf(cesa = SPOP2, maf = MAF5, maf_name = "341", coverage = "targeted", 
                 covered_regions = "msk_341_exons.bed", 
                 covered_regions_name = "MSK_IMPACT_341", covered_regions_padding = 10)
SPOP2 <- load_maf(cesa = SPOP2, maf = MAF6, maf_name = "410", coverage = "targeted", 
                 covered_regions = "msk_410_exons.bed", 
                 covered_regions_name = "MSK_IMPACT_410", covered_regions_padding = 10)
SPOP2 <- load_maf(cesa = SPOP2, maf = MAF7, maf_name = "468", coverage = "targeted", 
                 covered_regions = "msk_468_exons.bed", 
                 covered_regions_name = "MSK_IMPACT_468", covered_regions_padding = 10)

saveRDS(SPOP2, file = "SPOP2.rds")

##Mutational processes and relative mutation rates:
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

SPOP3 = trinuc_mutation_rates(SPOP2, ces.refset.hg19$signatures$COSMIC_v3.2,
                             signature_exclusions = signature_exclusions)


saveRDS(SPOP3, file = "SPOP3.rds")

SPOP4 = gene_mutation_rates(SPOP3, covariates = "PRAD")
saveRDS(SPOP4, file = "SPOP4.rds")

SPOP_final = ces_variant(SPOP4, variants = select_variants(SPOP4, min_freq = 2))
saveRDS(SPOP_final, file="SPOP_final.rds")

### making the figure:

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

common.text.size <- 4

PRAD_analysis <- load_cesa("SPOP_final.rds")
PRAD_results <- snv_results(PRAD_analysis)
PRAD_results <- PRAD_results$selection.1

# extract and add the gene names:
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
  geom_text(aes(label=included_with_variant, y= -3000), size = common.text.size) +  # adjusted y value
  scale_y_continuous(labels=scientific, breaks = c(0, 2e4, 4e4, 6e4, 8e4), limits = c(-3e3, 1e5))

setwd("C:/Moein/projects/prostate_stages/PRAD_files/PRAD_figures/Figures")

ggsave("SPOP_recurrent_SI.png", width=8, height=5.25)

#End

###Overlay SPOP_labeled_2.png onto SPOP_recurrent_SI.png:

library(magick)

# Set the file paths for the two PNG images you want to merge
file1 <- "SPOP_recurrent_SI.png"
file2 <- "SPOP_model_recurrent_resized.png"

# Read the images
image1 <- image_read(file1)
image2 <- image_read(file2)

# Get the dimensions of image1
width1 <- image_info(image1)$width
height1 <- image_info(image1)$height

# Calculate the offset to position image2 in the upper right corner of image1
offset_x <- width1 - image_info(image2)$width
offset_y <- 0

# Overlay image2 onto image1
merged_image <- image_composite(image1, image2, offset = paste0("+", offset_x, "+", offset_y))


# Write the merged image to a file
image_write(merged_image, "Figure_5_SPOP.png")

#End
