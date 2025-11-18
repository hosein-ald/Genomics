Reproducible ChIP-seq Analysis of YAP1, TAZ, and TEAD4
Reproducing Figures from Zanconato et al., Cell (2018) — Using CrazyHotTommy’s Tutorial
1. Project Overview

This repository documents my full attempt to reproduce Figure 1 and related analyses from the paper:

Zanconato et al. (2018)
YAP/TAZ–TEAD and AP-1 Cooperate to Control Oncogenic Enhancers in Cancer Cells
Link (Original Article): https://pmc.ncbi.nlm.nih.gov/articles/PMC6186417/

I followed (and adapted) the excellent reproducible-analysis workflow provided by:

CrazyHotTommy Reproducible Genomics Tutorial
Link (Tutorial): https://crazyhottommy.github.io/reproduce_genomics_paper_figures/index.html

This tutorial is the guide, and the PMC article is the original study whose figures are being reproduced.

My reproduced versions of the figures—including heatmaps, enhancer/promoter partitions, distance-to-summit distributions, peak overlap Venn diagrams, IGV snapshots, and ChIP-seq signal correlations—are all included in this repository.

2. Data Used

I used the ChIP-seq datasets described in the tutorial and paper, focusing on:

YAP1 ChIP-seq

TAZ ChIP-seq

TEAD4 ChIP-seq

IgG control

All datasets were aligned, filtered, and processed into the following key files:

.bam (aligned reads)

.bw (normalized tracks for IGV visualization)

.narrowPeak (MACS2 peak calls)
