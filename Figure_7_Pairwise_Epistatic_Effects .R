library(cancereffectsizeR)
library(data.table)
library(ces.refset.hg19)
library(ggplot2)
library(stringr)
library(scales)
library(cowplot)
library(grid)


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


# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "PRAD")

# estimating trinucleotide mutation rates
cesa <- trinuc_mutation_rates(cesa = cesa, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)

#cesa <- clear_gene_rates(cesa_samples_by_groups)
cesa <- gene_mutation_rates(cesa, covariates = "PRAD", save_all_dndscv_output = T)

selected_genes <- c("SPOP", "FOXA1", "AR", "PIK3CA", "PIK3CB", "TP53", "ROCK1", "RHOA", "AKT1", "ATM", "CUL3",
                    "APC", "CTNNB1", "PTEN", "KMT2C", "KMT2D")


consensus_gr = Reduce(GenomicRanges::intersect, unlist(cesa$coverage_ranges[c('exome', 'targeted')]))
recurrent_cons_cvg = select_variants(cesa, gr = consensus_gr)[maf_prevalence > 1]

cesa <- ces_gene_epistasis(cesa = cesa, genes = selected_genes, variants = recurrent_cons_cvg, run_name = "gene_epistasis_PRAD")

epistasis_results <- cesa$epistasis$gene_epistasis_PRAD

library(dplyr)
library(purrr)

# Compute Fisher's exact test p-values:
fisher_results <- epistasis_results %>%
  mutate(
    p_value = pmap_dbl(
      list(nAB, nA0, nB0, n00),
      function(nAB, nA0, nB0, n00) {
        mat <- matrix(c(nAB, nA0, nB0, n00), nrow = 2, byrow = TRUE)
        if (all(dim(mat) == c(2, 2))) {
          fisher.test(mat)$p.value
        } else {
          NA_real_
        }
      }
    )
  ) %>%
  select(
    variant_A,  # Gene A
    variant_B,  # Gene B
    p_A_change,
    p_B_change,
    p_epistasis,
    nA0,        # Gene A mutated
    nB0,        # Gene B mutated
    nAB,        # Both mutated
    n00,        # Neither mutated
    n_total,    # Total samples
    p_value     # Fisher's exact test p-value
  )

genes_of_interest <- c("AR", "TP53", "PIK3CA", "SPOP")

significant_epistasis_filtered <- fisher_results %>%
  filter(
    (variant_A %in% genes_of_interest | variant_B %in% genes_of_interest) &
      (p_epistasis < 0.05 | p_A_change < 0.05 | p_B_change < 0.05 | p_value < 0.05)
  )


scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}


###SPOP###

# Select gene pairs involving SPOP
SPOP_list <- which(epistasis_results$variant_A == "SPOP" | epistasis_results$variant_B == "SPOP")

epistatic_change_SPOP <- c()

# Extract pairwise epistatic changes
for(x in SPOP_list){
  g1 <- as.character(epistasis_results[x,1])
  g2 <- as.character(epistasis_results[x,2])
  val_g1 <- as.numeric(epistasis_results[x,3])
  val_g2 <- as.numeric(epistasis_results[x,4])
  val_g1_after_g2 <- as.numeric(epistasis_results[x,5])
  val_g2_after_g1 <- as.numeric(epistasis_results[x,6])
  
  gene1_after_gene2 <- c(paste0(g1, "_after_", g2), val_g1_after_g2 - val_g1)
  gene2_after_gene1 <- c(paste0(g2, "_after_", g1), val_g2_after_g1 - val_g2)
  
  epistatic_change_SPOP <- rbind(epistatic_change_SPOP, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_SPOP <- data.frame(
  gene_pair = epistatic_change_SPOP[,1],
  change = as.numeric(epistatic_change_SPOP[,2])
)

# Label time and extract interacting gene
epistatic_change_SPOP$time <- ifelse(grepl("^SPOP_after_", epistatic_change_SPOP$gene_pair), "Before", "After")
epistatic_change_SPOP$gene <- ifelse(
  epistatic_change_SPOP$time == "Before",
  sub("SPOP_after_", "", epistatic_change_SPOP$gene_pair),
  sub("_after_SPOP", "", epistatic_change_SPOP$gene_pair)
)

# Add significance symbols
annotations_SPOP <- fisher_results %>%
  filter(variant_A == "SPOP" | variant_B == "SPOP") %>%
  mutate(
    gene1_after_gene2 = paste0(variant_A, "_after_", variant_B),
    gene2_after_gene1 = paste0(variant_B, "_after_", variant_A),
    symbol1 = ifelse(p_value < 0.05, "*", ""),
    symbol1 = ifelse(p_epistasis < 0.05, paste0(symbol1, "†"), symbol1),
    symbol2 = ifelse(p_value < 0.05, "*", ""),
    symbol2 = ifelse(p_epistasis < 0.05, paste0(symbol2, "†"), symbol2)
  ) %>%
  select(gene1_after_gene2, gene2_after_gene1, symbol1, symbol2)

epistatic_change_SPOP$symbol <- ""
epistatic_change_SPOP$symbol <- ifelse(
  epistatic_change_SPOP$gene_pair %in% annotations_SPOP$gene1_after_gene2,
  annotations_SPOP$symbol1[match(epistatic_change_SPOP$gene_pair, annotations_SPOP$gene1_after_gene2)],
  epistatic_change_SPOP$symbol
)
epistatic_change_SPOP$symbol <- ifelse(
  epistatic_change_SPOP$gene_pair %in% annotations_SPOP$gene2_after_gene1,
  annotations_SPOP$symbol2[match(epistatic_change_SPOP$gene_pair, annotations_SPOP$gene2_after_gene1)],
  epistatic_change_SPOP$symbol
)

# Apply your original factor ordering
epistatic_change_SPOP$gene <- ifelse(epistatic_change_SPOP$time == "Before",
                                     epistatic_change_SPOP$gene,
                                     paste0(epistatic_change_SPOP$gene, "_"))

epistatic_change_SPOP$gene <- factor(epistatic_change_SPOP$gene,
                                     levels = c("KMT2C", "AKT1", "CUL3", "PIK3CB", "APC", "ATM", "KMT2D", "PTEN",
                                                "FOXA1", "TP53", "CTNNB1", "AR", "PIK3CA", "RHOA",
                                                "CUL3_", "PIK3CB_", "PIK3CA_", "AR_", "ATM_",
                                                "APC_", "TP53_", "KMT2D_", "CTNNB1_", "KMT2C_", "PTEN_", "FOXA1_",
                                                "AKT1_", "RHOA_")
)

# Clean x-axis labels (remove _ for display only)
gene_labels_SPOP <- sub("_", "", levels(epistatic_change_SPOP$gene))

# Plot
waterfall_SPOP <- ggplot(epistatic_change_SPOP, aes(x = gene, y = change, fill = time)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_text(aes(label = symbol), position = position_dodge(width = 0.9),
            vjust = ifelse(epistatic_change_SPOP$change >= 0, -0.5, 1.2), size = 4) +
  theme_classic() +
  scale_fill_manual(values = c("Before" = "#F8766D", "After" = "#00BFC4")) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  xlab("Gene") + ylab("Epistatic change in selection") +
  scale_x_discrete(labels = gene_labels_SPOP) +
  scale_y_continuous(labels = scientific, limits = c(-4e4, 8e4),
                     breaks = c(-4e4, -2e4, 0, 2e4, 4e4, 6e4, 8e4)) +
  geom_hline(yintercept = 0)

# Save the plot
ggsave("PRAD_figures/epistasis_16/waterfall_SPOP.png", width = 10, dpi = 300, height = 7)

###PIK3CA###
# Select relevant gene pairs
PIK3CA_list <- which(epistasis_results$variant_A == "PIK3CA" | epistasis_results$variant_B == "PIK3CA")

epistatic_change_PIK3CA <- c()

# Extract and label directional changes
for(x in PIK3CA_list){
  g1 <- as.character(epistasis_results[x,1])
  g2 <- as.character(epistasis_results[x,2])
  val_g1 <- as.numeric(epistasis_results[x,3])
  val_g2 <- as.numeric(epistasis_results[x,4])
  val_g1_after_g2 <- as.numeric(epistasis_results[x,5])
  val_g2_after_g1 <- as.numeric(epistasis_results[x,6])
  
  gene1_after_gene2 <- c(paste0(g1, "_after_", g2), val_g1_after_g2 - val_g1)
  gene2_after_gene1 <- c(paste0(g2, "_after_", g1), val_g2_after_g1 - val_g2)
  
  epistatic_change_PIK3CA <- rbind(epistatic_change_PIK3CA, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_PIK3CA <- data.frame(
  gene_pair = epistatic_change_PIK3CA[,1],
  change = as.numeric(epistatic_change_PIK3CA[,2])
)

# Label Before/After and extract target gene
epistatic_change_PIK3CA$time <- ifelse(grepl("^PIK3CA_after_", epistatic_change_PIK3CA$gene_pair), "Before", "After")
epistatic_change_PIK3CA$gene <- ifelse(
  epistatic_change_PIK3CA$time == "Before",
  sub("PIK3CA_after_", "", epistatic_change_PIK3CA$gene_pair),
  sub("_after_PIK3CA", "", epistatic_change_PIK3CA$gene_pair)
)

# Add significance annotations
annotations_PIK3CA <- fisher_results %>%
  filter(variant_A == "PIK3CA" | variant_B == "PIK3CA") %>%
  mutate(
    gene1_after_gene2 = paste0(variant_A, "_after_", variant_B),
    gene2_after_gene1 = paste0(variant_B, "_after_", variant_A),
    symbol1 = ifelse(p_value < 0.05, "*", ""),
    symbol1 = ifelse(p_epistasis < 0.05, paste0(symbol1, "†"), symbol1),
    symbol2 = ifelse(p_value < 0.05, "*", ""),
    symbol2 = ifelse(p_epistasis < 0.05, paste0(symbol2, "†"), symbol2)
  ) %>%
  select(gene1_after_gene2, gene2_after_gene1, symbol1, symbol2)

epistatic_change_PIK3CA$symbol <- ""
epistatic_change_PIK3CA$symbol <- ifelse(
  epistatic_change_PIK3CA$gene_pair %in% annotations_PIK3CA$gene1_after_gene2,
  annotations_PIK3CA$symbol1[match(epistatic_change_PIK3CA$gene_pair, annotations_PIK3CA$gene1_after_gene2)],
  epistatic_change_PIK3CA$symbol
)
epistatic_change_PIK3CA$symbol <- ifelse(
  epistatic_change_PIK3CA$gene_pair %in% annotations_PIK3CA$gene2_after_gene1,
  annotations_PIK3CA$symbol2[match(epistatic_change_PIK3CA$gene_pair, annotations_PIK3CA$gene2_after_gene1)],
  epistatic_change_PIK3CA$symbol
)

# Apply original ordering and unique labels
epistatic_change_PIK3CA$gene <- ifelse(epistatic_change_PIK3CA$time == "Before",
                                       epistatic_change_PIK3CA$gene,
                                       paste0(epistatic_change_PIK3CA$gene, "_"))

epistatic_change_PIK3CA$gene <- factor(epistatic_change_PIK3CA$gene,
                                       levels = c("SPOP", "AR", "CTNNB1", "KMT2C", "PIK3CB", "AKT1", "CUL3", "APC",
                                                  "KMT2D", "ATM", "RHOA", "FOXA1", "TP53", "PTEN",
                                                  "CUL3_", "RHOA_", "AKT1_", "PIK3CB_", "ATM_",
                                                  "APC_", "TP53_", "PTEN_", "KMT2D_", "KMT2C_", "AR_", "FOXA1_",
                                                  "SPOP_", "CTNNB1_"))

gene_labels_PIK3CA <- sub("_", "", levels(epistatic_change_PIK3CA$gene))

# Plot
waterfall_PIK3CA <- ggplot(epistatic_change_PIK3CA, aes(x = gene, y = change, fill = time)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_text(aes(label = symbol), position = position_dodge(width = 0.9),
            vjust = ifelse(epistatic_change_PIK3CA$change >= 0, -0.5, 1.2), size = 4) +
  theme_classic() +
  scale_fill_manual(values = c("Before" = "#F8766D", "After" = "#00BFC4")) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  xlab("Gene") + ylab("Epistatic change in selection") +
  scale_x_discrete(labels = gene_labels_PIK3CA) +
  scale_y_continuous(labels = scientific, limits = c(-4e4, 5.5e4),
                     breaks = c(-4e4, -3e4, -2e4, -1e4, 0, 1e4, 2e4, 3e4, 4e4, 5e4)) +
  geom_hline(yintercept = 0)

# Save
ggsave("PRAD_figures/epistasis_16/waterfall_PIK3CA.png", width = 10, dpi=300, height = 7)

###TP53###
# Identify relevant gene pairs
TP53_list <- which(epistasis_results$variant_A == "TP53" | epistasis_results$variant_B == "TP53")

epistatic_change_TP53 <- c()

# Extract epistatic change values
for(x in TP53_list){
  g1 <- as.character(epistasis_results[x,1])
  g2 <- as.character(epistasis_results[x,2])
  val_g1 <- as.numeric(epistasis_results[x,3])
  val_g2 <- as.numeric(epistasis_results[x,4])
  val_g1_after_g2 <- as.numeric(epistasis_results[x,5])
  val_g2_after_g1 <- as.numeric(epistasis_results[x,6])
  
  gene1_after_gene2 <- c(paste0(g1, "_after_", g2), val_g1_after_g2 - val_g1)
  gene2_after_gene1 <- c(paste0(g2, "_after_", g1), val_g2_after_g1 - val_g2)
  
  epistatic_change_TP53 <- rbind(epistatic_change_TP53, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_TP53 <- data.frame(
  gene_pair = epistatic_change_TP53[,1],
  change = as.numeric(epistatic_change_TP53[,2])
)

# Label directionality and extract gene name
epistatic_change_TP53$time <- ifelse(grepl("^TP53_after_", epistatic_change_TP53$gene_pair), "Before", "After")
epistatic_change_TP53$gene <- ifelse(
  epistatic_change_TP53$time == "Before",
  sub("TP53_after_", "", epistatic_change_TP53$gene_pair),
  sub("_after_TP53", "", epistatic_change_TP53$gene_pair)
)

# Annotate significance symbols
annotations_TP53 <- fisher_results %>%
  filter(variant_A == "TP53" | variant_B == "TP53") %>%
  mutate(
    gene1_after_gene2 = paste0(variant_A, "_after_", variant_B),
    gene2_after_gene1 = paste0(variant_B, "_after_", variant_A),
    symbol1 = ifelse(p_value < 0.05, "*", ""),
    symbol1 = ifelse(p_epistasis < 0.05, paste0(symbol1, "†"), symbol1),
    symbol2 = ifelse(p_value < 0.05, "*", ""),
    symbol2 = ifelse(p_epistasis < 0.05, paste0(symbol2, "†"), symbol2)
  ) %>%
  select(gene1_after_gene2, gene2_after_gene1, symbol1, symbol2)

epistatic_change_TP53$symbol <- ""
epistatic_change_TP53$symbol <- ifelse(
  epistatic_change_TP53$gene_pair %in% annotations_TP53$gene1_after_gene2,
  annotations_TP53$symbol1[match(epistatic_change_TP53$gene_pair, annotations_TP53$gene1_after_gene2)],
  epistatic_change_TP53$symbol
)
epistatic_change_TP53$symbol <- ifelse(
  epistatic_change_TP53$gene_pair %in% annotations_TP53$gene2_after_gene1,
  annotations_TP53$symbol2[match(epistatic_change_TP53$gene_pair, annotations_TP53$gene2_after_gene1)],
  epistatic_change_TP53$symbol
)

# Restore original gene ordering with underscores for "After"
epistatic_change_TP53$gene <- ifelse(epistatic_change_TP53$time == "Before",
                                     epistatic_change_TP53$gene,
                                     paste0(epistatic_change_TP53$gene, "_"))

epistatic_change_TP53$gene <- factor(epistatic_change_TP53$gene,
                                     levels = c("PIK3CA", "AKT1", "PIK3CB", "KMT2D", "ATM", "RHOA", "SPOP", "AR",
                                                "PTEN", "FOXA1", "CTNNB1", "CUL3", "KMT2C", "APC",
                                                "CUL3_", "RHOA_", "CTNNB1_", "SPOP_", "APC_", "AR_", "FOXA1_",
                                                "PTEN_", "KMT2D_", "KMT2C_", "PIK3CA_", "AKT1_", "PIK3CB_", "ATM_"))

# Clean axis labels
gene_labels_TP53 <- sub("_", "", levels(epistatic_change_TP53$gene))

# Plot
waterfall_TP53 <- ggplot(epistatic_change_TP53, aes(x = gene, y = change, fill = time)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_text(aes(label = symbol), position = position_dodge(width = 0.9),
            vjust = ifelse(epistatic_change_TP53$change >= 0, -0.5, 1.2), size = 4) +
  theme_classic() +
  scale_fill_manual(values = c("Before" = "#F8766D", "After" = "#00BFC4")) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  xlab("Gene") + ylab("Epistatic change in selection") +
  scale_x_discrete(labels = gene_labels_TP53) +
  scale_y_continuous(labels = scientific, limits = c(-3.5e4,1e4),
                     breaks = c(-3e4, -2e4, -1e4, 0, 1e4)) +
  geom_hline(yintercept = 0)

# Save plot
ggsave("PRAD_figures/epistasis_16/waterfall_TP53.png", width = 10, dpi = 300, height = 7)

###AR###
# Select gene pairs involving AR
AR_list <- which(epistasis_results$variant_A == "AR" | epistasis_results$variant_B == "AR")

epistatic_change_AR <- c()

# Extract directional epistatic change values
for(x in AR_list){
  g1 <- as.character(epistasis_results[x,1])
  g2 <- as.character(epistasis_results[x,2])
  val_g1 <- as.numeric(epistasis_results[x,3])
  val_g2 <- as.numeric(epistasis_results[x,4])
  val_g1_after_g2 <- as.numeric(epistasis_results[x,5])
  val_g2_after_g1 <- as.numeric(epistasis_results[x,6])
  
  gene1_after_gene2 <- c(paste0(g1, "_after_", g2), val_g1_after_g2 - val_g1)
  gene2_after_gene1 <- c(paste0(g2, "_after_", g1), val_g2_after_g1 - val_g2)
  
  epistatic_change_AR <- rbind(epistatic_change_AR, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_AR <- data.frame(
  gene_pair = epistatic_change_AR[,1],
  change = as.numeric(epistatic_change_AR[,2])
)

# Label "Before" and "After", and extract partner gene
epistatic_change_AR$time <- ifelse(grepl("^AR_after_", epistatic_change_AR$gene_pair), "Before", "After")
epistatic_change_AR$gene <- ifelse(
  epistatic_change_AR$time == "Before",
  sub("AR_after_", "", epistatic_change_AR$gene_pair),
  sub("_after_AR", "", epistatic_change_AR$gene_pair)
)

# Add significance annotations
annotations_AR <- fisher_results %>%
  filter(variant_A == "AR" | variant_B == "AR") %>%
  mutate(
    gene1_after_gene2 = paste0(variant_A, "_after_", variant_B),
    gene2_after_gene1 = paste0(variant_B, "_after_", variant_A),
    symbol1 = ifelse(p_value < 0.05, "*", ""),
    symbol1 = ifelse(p_epistasis < 0.05, paste0(symbol1, "†"), symbol1),
    symbol2 = ifelse(p_value < 0.05, "*", ""),
    symbol2 = ifelse(p_epistasis < 0.05, paste0(symbol2, "†"), symbol2)
  ) %>%
  select(gene1_after_gene2, gene2_after_gene1, symbol1, symbol2)

epistatic_change_AR$symbol <- ""
epistatic_change_AR$symbol <- ifelse(
  epistatic_change_AR$gene_pair %in% annotations_AR$gene1_after_gene2,
  annotations_AR$symbol1[match(epistatic_change_AR$gene_pair, annotations_AR$gene1_after_gene2)],
  epistatic_change_AR$symbol
)
epistatic_change_AR$symbol <- ifelse(
  epistatic_change_AR$gene_pair %in% annotations_AR$gene2_after_gene1,
  annotations_AR$symbol2[match(epistatic_change_AR$gene_pair, annotations_AR$gene2_after_gene1)],
  epistatic_change_AR$symbol
)

# Use original factor ordering and encode After with "_"
epistatic_change_AR$gene <- ifelse(epistatic_change_AR$time == "Before",
                                   epistatic_change_AR$gene,
                                   paste0(epistatic_change_AR$gene, "_"))

epistatic_change_AR$gene <- factor(epistatic_change_AR$gene,
                                   levels = c("SPOP", "CTNNB1", "KMT2C", "AKT1", "PIK3CB", "CUL3", "APC", "RHOA",
                                              "ATM", "TP53", "FOXA1", "PIK3CA", "PTEN", "KMT2D",
                                              "CUL3_", "RHOA_", "AKT1_", "PIK3CB_", "PIK3CA_", "APC_", "FOXA1_",
                                              "PTEN_", "KMT2D_", "TP53_", "KMT2C_", "CTNNB1_", "SPOP_", "ATM_"))

# Clean gene labels
gene_labels_AR <- sub("_", "", levels(epistatic_change_AR$gene))

# Plot
waterfall_AR <- ggplot(epistatic_change_AR, aes(x = gene, y = change, fill = time)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_text(aes(label = symbol), position = position_dodge(width = 0.9),
            vjust = ifelse(epistatic_change_AR$change >= 0, -0.5, 1.2), size = 4) +
  theme_classic() +
  scale_fill_manual(values = c("Before" = "#F8766D", "After" = "#00BFC4")) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size = 12)) +
  xlab("Gene") + ylab("Epistatic change in selection") +
  scale_x_discrete(labels = gene_labels_AR) +
  scale_y_continuous(labels = scientific, limits = c(-5e4, 8.5e4),
                     breaks = c(-4e4, -2e4, 0, 2e4, 4e4, 6e4, 8e4)) +
  geom_hline(yintercept = 0)

# Save plot
ggsave("PRAD_figures/epistasis_16/waterfall_AR.png", width = 10, dpi = 300, height = 7)

###############################################

combined_waterfall <- plot_grid(waterfall_SPOP, waterfall_PIK3CA, waterfall_TP53, waterfall_AR,
                                labels = c("A", "B", "C", "D"), ncol = 2)


pdf("combined_waterfall_with_legends.pdf", width = 12, height = 8)
grid.newpage()
print(combined_waterfall, vp = viewport())


# SPOP
grid.rect(x = unit(0.125, "npc"), y = unit(0.96, "npc"), width = unit(0.018, "npc"), height = unit(0.028, "npc"), gp = gpar(fill = "#F8766D", col = NA))
grid.text("Change in selection for", x = unit(0.14, "npc"), y = unit(0.98, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression("mutated "*italic("SPOP")*" after"), x = unit(0.14, "npc"), y = unit(0.96, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text("the gene below is mutated", x = unit(0.14, "npc"), y = unit(0.94, "npc"), just = "left", gp = gpar(fontsize = 8))

grid.rect(x = unit(0.335, "npc"), y = unit(0.96, "npc"), width = unit(0.018, "npc"), height = unit(0.028, "npc"), gp = gpar(fill = "#00BFC4", col = NA))
grid.text("Change in selection for", x = unit(0.35, "npc"), y = unit(0.98, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression("mutated gene below after"), x = unit(0.35, "npc"), y = unit(0.96, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression(italic("SPOP")*" is mutated"), x = unit(0.35, "npc"), y = unit(0.94, "npc"), just = "left", gp = gpar(fontsize = 8))

# TP53
grid.rect(x = unit(0.125, "npc"), y = unit(0.24, "npc"), width = unit(0.018, "npc"), height = unit(0.028, "npc"), gp = gpar(fill = "#F8766D", col = NA))
grid.text("Change in selection for", x = unit(0.14, "npc"), y = unit(0.26, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression("mutated "*italic("TP53")*" after"), x = unit(0.14, "npc"), y = unit(0.24, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text("the gene below is mutated", x = unit(0.14, "npc"), y = unit(0.22, "npc"), just = "left", gp = gpar(fontsize = 8))

grid.rect(x = unit(0.335, "npc"), y = unit(0.24, "npc"), width = unit(0.018, "npc"), height = unit(0.028, "npc"), gp = gpar(fill = "#00BFC4", col = NA))
grid.text("Change in selection for", x = unit(0.35, "npc"), y = unit(0.26, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression("mutated gene below after"), x = unit(0.35, "npc"), y = unit(0.24, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression(italic("TP53")*" is mutated"), x = unit(0.35, "npc"), y = unit(0.22, "npc"), just = "left", gp = gpar(fontsize = 8))

# PIK3CA
grid.rect(x = unit(0.615, "npc"), y = unit(0.96, "npc"), width = unit(0.018, "npc"), height = unit(0.028, "npc"), gp = gpar(fill = "#F8766D", col = NA))
grid.text("Change in selection for", x = unit(0.63, "npc"), y = unit(0.98, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression("mutated "*italic("PIK3CA")*" after"), x = unit(0.63, "npc"), y = unit(0.96, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text("the gene below is mutated", x = unit(0.63, "npc"), y = unit(0.94, "npc"), just = "left", gp = gpar(fontsize = 8))

grid.rect(x = unit(0.825, "npc"), y = unit(0.96, "npc"), width = unit(0.018, "npc"), height = unit(0.028, "npc"), gp = gpar(fill = "#00BFC4", col = NA))
grid.text("Change in selection for", x = unit(0.84, "npc"), y = unit(0.98, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression("mutated gene below after"), x = unit(0.84, "npc"), y = unit(0.96, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression(italic("PIK3CA")*" is mutated"), x = unit(0.84, "npc"), y = unit(0.94, "npc"), just = "left", gp = gpar(fontsize = 8))

# AR
grid.rect(x = unit(0.615, "npc"), y = unit(0.46, "npc"), width = unit(0.02, "npc"), height = unit(0.03, "npc"), gp = gpar(fill = "#F8766D", col = NA))
grid.text("Change in selection for", x = unit(0.63, "npc"), y = unit(0.48, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression("mutated "*italic("AR")*" after"), x = unit(0.63, "npc"), y = unit(0.46, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text("the gene below is mutated", x = unit(0.63, "npc"), y = unit(0.44, "npc"), just = "left", gp = gpar(fontsize = 8))

grid.rect(x = unit(0.825, "npc"), y = unit(0.46, "npc"), width = unit(0.02, "npc"), height = unit(0.03, "npc"), gp = gpar(fill = "#00BFC4", col = NA))
grid.text("Change in selection for", x = unit(0.84, "npc"), y = unit(0.48, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression("mutated gene below after"), x = unit(0.84, "npc"), y = unit(0.46, "npc"), just = "left", gp = gpar(fontsize = 8))
grid.text(expression(italic("AR")*" is mutated"), x = unit(0.84, "npc"), y = unit(0.44, "npc"), just = "left", gp = gpar(fontsize = 8))

dev.off()

#End
