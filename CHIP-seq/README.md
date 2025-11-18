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
```

##3. Pipeline Summary
3.1 Alignment
```bash
bowtie2 -x hg38 \
  -1 sample_R1.fq.gz -2 sample_R2.fq.gz \
  | samtools sort -o sample.sorted.bam
samtools index sample.sorted.bam
```

##3.2 Peak Calling
```bash
macs3 callpeak -t sample.sorted.bam \
  -c IgG.sorted.bam \
  -n sample \
  -f BAM
```

##3.3 BigWig Generation
```bash
bamCoverage -b sample.sorted.bam \
  -o sample.bw \
  --normalizeUsing RPKM
```

##3.4 Heatmaps (deepTools)
```bash
computeMatrix reference-point \
  -S YAP.bw TAZ.bw TEAD4.bw \
  -R peaks.bed \
  --referencePoint center \
  -b 3000 -a 3000 \
  -o matrix.gz

plotHeatmap -m matrix.gz -o figures/heatmap.png
```

##3.5 R / ggplot2 Analysis

All downstream analyses (TSS distance classification, summit distance distribution, scatter plots, and Venn diagrams) were performed in R and are available in:

scripts/04_R_analysis.R

notebooks/YAP_TAZ_TEAD4_reproduction.ipynb
