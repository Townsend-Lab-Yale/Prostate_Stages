This repository provides the data and code used in the project "Somatic evolution of prostate cancer: mutation, selection, and epistasis across disease stages".
Large-Scale Prostate Cancer Evolution

Computational analysis of somatic evolution across primary and metastatic prostate cancer using large-scale, multi-cohort genomic datasets.

This repository contains the R analysis and visualization workflows supporting our study:

Rajaei M., Yang A., Cross C.N., Glasmacher K., Fisk J.N., Perry E.B., Mandell J.D., Gaffney S.G., Yamaguchi T.N., Livingstone J., Costa J., Humphrey P., Cannataro V.L., Boutros P.C., Townsend J.P.
Somatic evolution of prostate cancer: mutation, selection, and epistasis across disease stages.
bioRxiv (2025).

Project Overview

Prostate cancer progresses through multiple stages, from the emergence of a primary tumor to advanced and metastatic disease. The evolutionary forces driving these transitions can differ substantially across disease stages.

In this project, we integrated genomic data from 2,618 prostate tumors, including 1,566 primary tumors and 1,052 metastatic castration-resistant prostate cancers (mCRPC), to investigate how mutation and selection shape prostate cancer evolution.

The study combines tumors profiled using:

Whole-genome sequencing (WGS)
Whole-exome sequencing (WES)
Targeted sequencing

The computational framework separates the contribution of mutation rate from positive selection, allowing us to investigate how cancer driver effects change during disease progression.

We additionally evaluate genetic interactions (epistasis) among prostate cancer drivers to determine how earlier genomic alterations modify selection on subsequent driver mutations.

Research Questions

The analyses address several major questions:

How do somatic mutation profiles change across stages of prostate cancer?
Which driver genes experience the strongest positive selection during primary tumor development and metastatic progression?
How do mutation rates and cancer effect sizes differ across disease stages?
Which driver mutations are preferentially selected early versus late during tumor evolution?
How does the presence of one driver mutation alter selection on subsequent mutations?
Which genetic interactions may contribute to progression toward metastatic and treatment-resistant disease?
Dataset

The integrated dataset contains 2,618 prostate tumors:

Disease group	Number of tumors
Primary prostate cancer	1,566
Metastatic castration-resistant prostate cancer (mCRPC)	1,052
Total	2,618

Primary tumors were additionally stratified by Gleason Grade Group to investigate evolutionary differences across disease severity.

The combined dataset incorporates multiple independent prostate cancer cohorts, including data from:

TCGA
MSK cohorts
SU2C/PCF
Armenia cohort
Boutros cohort
Additional published prostate cancer sequencing studies

Because these studies were generated using different sequencing technologies and capture designs, cohort-specific genomic coverage was incorporated where necessary.

Sequencing Platforms

The combined dataset includes:

Sequencing strategy	Approximate number of tumors
Whole-genome sequencing	1,036
Whole-exome sequencing	1,092
Targeted sequencing	570

Integrating these datasets required harmonization of genomic coordinates, mutation calls, sample identifiers, disease-stage annotations, and genomic regions accessible to each sequencing strategy.

Analytical Framework

The overall workflow can be summarized as:

Multiple prostate cancer cohorts
              │
              ▼
    Mutation data harmonization
              │
              ▼
   Disease-stage classification
              │
              ▼
  Mutation prevalence analysis
              │
              ▼
Trinucleotide mutation profiles
              │
              ▼
 Gene-specific mutation rates
              │
              ▼
   Cancer effect size analysis
              │
        ┌─────┴─────┐
        ▼           ▼
 Stage-specific   Driver-specific
   selection       selection
        │           │
        └─────┬─────┘
              ▼
    Pairwise epistasis
              │
              ▼
 Evolutionary interpretation
Major Analyses
Mutation Prevalence

Figure_1_prevalence.R

Compares the prevalence of somatic alterations across prostate cancer disease stages and tumor groups.

These analyses provide an initial description of how frequently major genomic alterations occur during prostate cancer progression.

Trinucleotide Mutation Profiles

Figure_2_trinucleotide_mutation_profiles.R

Characterizes trinucleotide mutation profiles across disease stages.

Mutation-spectrum analysis helps determine whether observed differences in driver frequencies arise from changes in underlying mutational processes, changes in selection, or both.

Gene-Specific Mutation Rates

Figure_3_gene_mutrate.R

Estimates gene-specific mutation rates while accounting for the mutational opportunities available within each dataset.

Separating mutation rate from selection is important because high mutation prevalence alone does not necessarily imply strong positive selection.

Cancer Effect Sizes

Figure_4_CES.R

Quantifies the selective effects of recurrent cancer driver mutations across prostate cancer stages.

Cancer effect sizes provide an estimate of the strength of positive selection acting on specific driver alterations after accounting for their expected mutation rates.

SPOP Evolution

Figure_5_SPOP.R

Examines stage-specific selection and mutation patterns involving SPOP, an important early prostate cancer driver.

The analysis evaluates how selection on SPOP-associated alterations changes across prostate cancer development and progression.

Related visualization:

SPOP_model_recurrent_resized.png
Androgen Receptor Evolution

Figure_6_AR.R

Examines evolutionary patterns involving the androgen receptor (AR), a major driver of advanced prostate cancer and treatment resistance.

The analysis evaluates changes in selection on AR during progression toward metastatic castration-resistant disease.

Related visualization:

AR_labeled_2.png
Pairwise Epistatic Effects

Figure_7_Pairwise_Epistatic_Effects.R

Quantifies pairwise genetic interactions among recurrent prostate cancer drivers.

These analyses evaluate whether the presence of one driver alteration changes the selective advantage associated with acquiring another driver alteration.

This framework allows the study of historical contingency in cancer evolution, where the fitness effect of a mutation depends on the genomic background in which it occurs.

Key Biological Findings

The analyses reveal substantial changes in somatic selection during prostate cancer progression.

Major observations include:

SPOP experiences strong selection during early prostate tumorigenesis.
Selection on AR becomes substantially stronger during progression toward metastatic castration-resistant disease.
Driver genes including APC, ROCK1, RHOA, and PIK3CB show stage-dependent patterns of selection associated with progression of primary disease.
Genetic background can substantially modify the selective advantage of subsequent driver mutations.
Alterations in PTEN, for example, modify selection on additional signaling and prostate cancer driver genes.
Pairwise epistasis analyses identify strong interactions among recurrent prostate cancer drivers, demonstrating that the evolutionary effect of a mutation can depend on mutations already present in the tumor.

Together, these results demonstrate that prostate cancer evolution reflects both changing mutation processes and changing selective pressures across disease stages.

Repository Contents
Main Analysis Scripts
File	Analysis
Figure_1_prevalence.R	Somatic mutation prevalence across disease groups
Figure_2_trinucleotide_mutation_profiles.R	Trinucleotide mutation-spectrum analysis
Figure_3_gene_mutrate.R	Gene-specific mutation-rate analysis
Figure_4_CES.R	Cancer effect size and selection analysis
Figure_5_SPOP.R	SPOP-specific evolutionary analysis
Figure_6_AR.R	Androgen receptor evolutionary analysis
Figure_7_Pairwise_Epistatic_Effects.R	Pairwise epistasis among prostate cancer drivers
Figure_S1.R	Supplementary analysis
Figure_S2.R	Supplementary analysis
Figure_S3.R	Supplementary analysis
new_sequential_lik.R	Sequential likelihood/statistical modeling functions
Mutation Datasets

The repository contains harmonized mutation data from multiple prostate cancer cohorts, including:

MSK_341_final.maf.txt
MSK_410_final.maf.txt
MSK_468_final.maf.txt
SU2C_final.maf.txt
tcga_final.maf.txt
tcga_wgs_final.maf.txt
armenia_final.maf.txt
boutros_final.maf.txt
prad_armenia_final.maf.txt
prad_boutros_wgs_final.maf.txt

Mutation data are primarily represented using Mutation Annotation Format (MAF) files.

Genomic Coverage Files
SureSelect_All_Exon_covered_regions.bed
msk_341_exons.bed
msk_410_exons.bed
msk_468_exons.bed

These BED files describe genomic regions accessible to specific sequencing or capture platforms and are used when accounting for differences in genomic coverage across cohorts.

Clinical and Disease-Stage Information
gleason.txt
gleason_age_comparison.txt

These files contain information used for stratification and comparison of primary prostate tumors.

Figures and Visualizations
AR_labeled_2.png
SPOP_model_recurrent_resized.png

These files contain selected visualizations generated from the analyses.

Computational Environment

The analyses were primarily developed using:

R
Linux
High-performance computing environments
Git/GitHub for version control

R packages and computational tools used across the project include packages for genomic data manipulation, statistical modeling, cancer evolutionary analysis, and visualization.

Reproducibility

The main figure scripts correspond directly to analyses presented in the associated manuscript.

A typical analysis workflow is:

1. Obtain and harmonize cohort-specific mutation data
2. Define genomic coverage for each sequencing platform
3. Assign samples to disease-stage groups
4. Analyze mutation prevalence
5. Calculate trinucleotide mutation profiles
6. Estimate gene-specific mutation rates
7. Quantify cancer effect sizes
8. Analyze individual prostate cancer drivers
9. Estimate pairwise epistatic effects
10. Generate manuscript figures

Some analyses depend on cohort-specific source data, clinical annotations, or intermediate files that may be subject to data-access restrictions and therefore may not be distributed directly through this repository.

Related Publication

Rajaei M, Yang A, Cross CN, Glasmacher K, Fisk JN, Perry EB, Mandell JD, Gaffney SG, Yamaguchi TN, Livingstone J, Costa J, Humphrey P, Cannataro VL, Boutros PC, Townsend JP.
Somatic evolution of prostate cancer: mutation, selection, and epistasis across disease stages.
bioRxiv. 2025.

Author

Moein Rajaei, Ph.D.
Computational Biologist | Cancer Genomics Data Scientist

Research areas:

Cancer genomics
Tumor evolution
Somatic mutation
Statistical genomics
Cancer driver selection
Genetic interactions
Reproducible computational biology

GitHub | ORCID | Google Scholar | LinkedIn

Citation

If you use code or analytical methods from this repository, please cite:

Rajaei M, et al. Somatic evolution of prostate cancer: mutation, selection, and epistasis across disease stages. bioRxiv. 2025.

License

Please see the repository license for information regarding reuse of code and materials.
