# Somatic evolution of prostate cancer: mutation, selection, and epistasis across disease stages.

Computational analysis of somatic evolution across primary and metastatic castration-resistant prostate cancer (mCRPC) using large-scale, multi-cohort genomic datasets.

This repository contains the **R analysis and visualization workflows** supporting our study:

**Rajaei M., Yang A., Cross C.N., Glasmacher K., Fisk J.N., Perry E.B., Mandell J.D., Gaffney S.G., Yamaguchi T.N., Livingstone J., Costa J., Humphrey P., Cannataro V.L., Boutros P.C., Townsend J.P.**
*Somatic evolution of prostate cancer: mutation, selection, and epistasis across disease stages.*
**bioRxiv (2025).**

## Project Overview

We analyzed 2,704 primary and mCRPC tumors from multiple cohorts and sequencing platforms to characterize how mutation, selection, and genetic interactions shape prostate cancer evolution.

The study quantifies mutation rates and somatic selection across disease stages and identifies synergistic and antagonistic epistatic interactions among key cancer drivers, revealing a dynamic evolutionary landscape from tumor initiation to mCRPC.

The study combines tumors profiled using:

* Whole-genome sequencing (WGS)
* Whole-exome sequencing (WES)
* Targeted sequencing

## Research Questions

The analyses address several major questions:

1. How do somatic mutation profiles change across stages of prostate cancer?
2. Which driver genes experience the strongest  selection during primary tumor development and metastatic progression?
3. How do mutation rates and cancer effect sizes differ across disease stages?
4. Which driver mutations are selected early versus late during tumor evolution?
5. How does the presence of one driver mutation alter selection on subsequent mutations?
6. Which genetic interactions may contribute to progression toward metastatic and treatment-resistant disease?

## Dataset

The integrated dataset contains **2,704 prostate tumors**. Samples were excluded if they contained only nucleotide base substitutions at known germline variant sites, within repetitive regions, or no nucleotide base substitutions, yielding a final set of 2,618 samples for downstream analysis

| Disease group                                           | Number of tumors |
| ------------------------------------------------------- | ---------------: |
| Primary prostate cancer                                 |            1,593 |
| Metastatic castration-resistant prostate cancer (mCRPC) |            1,025 |
| **Total**                                               |        **2,618** |

Primary tumors were additionally stratified by Gleason Grade Group to investigate evolutionary differences across disease severity.

The combined dataset incorporates multiple independent prostate cancer cohorts, including data from:

* TCGA
* MSK cohorts
* SU2C/PCF
* Armenia cohort
* Boutros cohort
* Additional published prostate cancer sequencing studies

## Sequencing Platforms

The combined dataset includes:

| Sequencing strategy     | Approximate number of tumors |
| ----------------------- | ---------------------------: |
| Whole-genome sequencing |                          293 |
| Whole-exome sequencing  |                        1,122 |
| Targeted sequencing     |                        1,289 |

## Analytical Framework

The overall workflow can be summarized as:

```text
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
```

## Major Analyses

### Mutation Prevalence

`Figure_1_prevalence.R`

Prevalence of variants in 16 selected driver genes, in low-grade primary tumors, high-grade primary tumors, and mCRPC.

### Trinucleotide Mutation Profiles

`Figure_2_trinucleotide_mutation_profiles.R`

Percent of single-nucleotide somatic variants within each trinucleotide context in low-grade primary tumors, high-grade primary tumors, and mCRPC.

### Gene-Specific Mutation Rates

`Figure_3_gene_mutrate.R`

The gene-level mutation rates spanning from organogenesis to low-grade primary tumors, organogenesis to high-grade primary tumors, organogenesis to mCRPC, and in  organogenesis to mCRPC versus organogenesis to low-grade primary tumorigenesis, organogenesis to mCRPC versus organogenesis to high-grade primary tumorigenesis, organogenesis to mCRPC versus primary tumorigenesis. 

### Cancer Effect Sizes

`Figure_4_CES.R`

Gene-level estimates and 95% confidence intervals for scaled selection coefficients on somatic variants in oncogenic sites of 16 genes that are known to act as drivers in prostate cancer tumorigenesis and metastasis.

### SPOP Evolution

`Figure_5_SPOP.R`

Scaled selection coefficients for recurrent single-nucleotide variant amino-acid substitutions in SPOP during the evolutionary trajectory from prostate organogenesis to primary and mCRPC tumors.

Related visualization:

```text
SPOP_model_recurrent_resized.png
```

### AR Evolution

`Figure_6_AR.R`

Scaled selection coefficients of recurrent single-nucleotide variant amino-acid substitutions in AR along the step from primary tumors to mCRPC.


Related visualization:

```text
AR_labeled_2.png
```

### Pairwise Epistatic Effects

`Figure_7_Pairwise_Epistatic_Effects.R`

Pairwise epistatic effect trends of 15 genes, focusing on SPOP, PIK3CA, TP53, and AR.

## Key Biological Findings

* Mutation load and mutation rates increase during prostate cancer progression, while **trinucleotide mutational patterns remain relatively stable**.
* **SPOP** mutations in the BRD3-binding domain experience strong early positive selection and increase subsequent selection for **RHOA** while decreasing selection for **TP53**.
* **CUL3** shows antagonistic selective epistasis with both **SPOP** and **PIK3CA**.
* **KMT2C** mutations increase selection for subsequent **TP53** mutations.
* **PTEN** mutations increase selection for both **PIK3CA** and **AR**, revealing strong synergistic epistatic interactions.

Together, these results reveal a dynamic landscape of **mutation, selection, and epistasis** across prostate cancer progression.


## Repository Contents

### Main Analysis Scripts

| File                                         | Analysis                                             |
| -------------------------------------------- | ---------------------------------------------------- |
| `Figure_1_prevalence.R`                      | Somatic mutation prevalence across disease groups    |
| `Figure_2_trinucleotide_mutation_profiles.R` | Trinucleotide mutation-spectrum analysis             |
| `Figure_3_gene_mutrate.R`                    | Gene-specific mutation-rate analysis                 |
| `Figure_4_CES.R`                             | Cancer effect size analysis                          |
| `Figure_5_SPOP.R`                            | SPOP-specific evolutionary analysis                  |
| `Figure_6_AR.R`                              | AR-specific evolutionary analysis                    |
| `Figure_7_Pairwise_Epistatic_Effects.R`      | Pairwise epistasis among prostate cancer drivers     |
| `Figure_S1.R`                                | Supplementary analysis                               |
| `Figure_S2.R`                                | Supplementary analysis                               |
| `Figure_S3.R`                                | Supplementary analysis                               |
| `new_sequential_lik.R`                       | Sequential likelihood/statistical modeling functions |

### Mutation Datasets

The repository contains harmonized mutation data from multiple prostate cancer cohorts, including:

```text
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
```

Mutation data are primarily represented using **Mutation Annotation Format (MAF)** files.

### Genomic Coverage Files

```text
SureSelect_All_Exon_covered_regions.bed
msk_341_exons.bed
msk_410_exons.bed
msk_468_exons.bed
```

### Clinical and Disease-Stage Information

```text
gleason.txt
gleason_age_comparison.txt
```

## Computational Environment

The analyses were primarily developed using:

* **R**
* Linux
* High-performance computing environments
* Git/GitHub for version control

## Reproducibility

The main figure scripts correspond directly to analyses presented in the associated manuscript.

A typical analysis workflow is:

```text
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
```

## Related Publication

**Rajaei M, Yang A, Cross CN, Glasmacher K, Fisk JN, Perry EB, Mandell JD, Gaffney SG, Yamaguchi TN, Livingstone J, Costa J, Humphrey P, Cannataro VL, Boutros PC, Townsend JP.**
*Somatic evolution of prostate cancer: mutation, selection, and epistasis across disease stages.*
**bioRxiv. 2025.**

## Author

**Moein Rajaei, Ph.D.**
Computational Biologist | Cancer Genomics & Bioinformatics

* GitHub: https://github.com/moeinrajaei
* ORCID: https://orcid.org/my-orcid?orcid=0009-0006-7079-2033
* Google Scholar: https://scholar.google.com/citations?user=LbgtrokAAAAJ&hl=en&oi=ao
* LinkedIn: https://www.linkedin.com/in/moein-rajaei-ph-d-3a5a136b/

## Citation

If you use code or analytical methods from this repository, please cite:

> Rajaei M, et al. *Somatic evolution of prostate cancer: mutation, selection, and epistasis across disease stages.* bioRxiv. 2025.
