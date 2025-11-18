# Reproducible ChIP-seq Analysis of YAP1, TAZ, and TEAD4  
### Reproducing Figures from *Zanconato et al., Cell* (2018) via CrazyHotTommy’s Tutorial

---

## 1. Project Overview

This repository documents my attempt to **reproduce Figure 1 and related analyses** from:

**Original article (paper)**  
Zanconato et al., 2018  
*YAP/TAZ–TEAD and AP-1 Cooperate to Control Oncogenic Enhancers in Cancer Cells*  
PMC link: https://pmc.ncbi.nlm.nih.gov/articles/PMC6186417/

The analysis is guided by the excellent reproducible genomics tutorial:

**Tutorial (workflow guide)**  
CrazyHotTommy – *Reproduce genomics paper figures*  
https://crazyhottommy.github.io/reproduce_genomics_paper_figures/index.html

The tutorial is the **step-by-step pipeline**, while the PMC article is the **original biological study** whose figures I reproduced.

This repo contains:

- All the **R / shell scripts** used for processing and plotting  
- A **Jupyter notebook** with narrative explanations  
- Exported **figures** closely matching those in the original paper

---

## 2. Data

I used the ChIP-seq datasets described in the paper and tutorial, focusing on:

- **YAP1 ChIP-seq**
- **TAZ ChIP-seq**
- **TEAD4 ChIP-seq**
- **IgG control**

Key derived files:

- `*.sorted.bam` – aligned, sorted reads (Bowtie2 + samtools)  
- `*.narrowPeak` – peak calls (MACS3)  
- `*.bw` – normalized signal tracks (deepTools `bamCoverage`)

Example MACS3 call:

```bash
macs3 callpeak -t <IP>.sorted.bam \
  -c IgG.sorted.bam \
  -f BAM \
  -n <sample_name>
