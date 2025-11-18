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

Example MACS3 call:

```bash
macs3 callpeak -t <IP>.sorted.bam \
  -c IgG.sorted.bam \
  -f BAM \
  -n <sample_name>
```

3. Pipeline Summary
3.1 Alignment
```bash
bowtie2 -x hg38 \
  -1 sample_R1.fq.gz -2 sample_R2.fq.gz \
  | samtools sort -o sample.sorted.bam
samtools index sample.sorted.bam
```

3.2 Peak Calling
```bash
macs3 callpeak -t sample.sorted.bam \
  -c IgG.sorted.bam \
  -n sample \
  -f BAM
```

3.3 BigWig Generation
```bash
bamCoverage -b sample.sorted.bam \
  -o sample.bw \
  --normalizeUsing RPKM
```

3.4 Heatmaps (deepTools)
```bash
computeMatrix reference-point \
  -S YAP.bw TAZ.bw TEAD4.bw \
  -R peaks.bed \
  --referencePoint center \
  -b 3000 -a 3000 \
  -o matrix.gz

plotHeatmap -m matrix.gz -o figures/heatmap.png
```

3.5 R / ggplot2 Analysis

All downstream analyses (TSS distance classification, summit distance distribution, scatter plots, and Venn diagrams) were performed in R and are available in:

scripts/04_R_analysis.R

notebooks/YAP_TAZ_TEAD4_reproduction.ipynb

4. Reproduced Figures (Overview)

All figure PNGs are stored in figures/.

Heatmap of ChIP-seq signal around peak summits

File: heatmap.png

YAP1, TAZ, and TEAD4 signal over ±1–3 kb around peak centers

Separate panels for promoters vs enhancers

Genomic partition of shared YAP/TAZ/TEAD4 peaks

File: pichart_with_percentage.png

Pie chart of promoters, active enhancers, inactive enhancers, unclassified

IGV snapshot of a representative locus (e.g. CCN2/CTGF)

File: Screenshot_2025-11-17_at_13.09.53.png

BigWig tracks for YAP, TAZ, TEAD4, IgG; narrowPeak beds shown below

Distance to summit distribution

File: distance_to_the_summit_of_TAZ_peaks.png

Histogram/density of read positions relative to TAZ peak summits

Distance to TSS classification

File: Distance_to_TSS.png

Stacked bar plot with categories <1 kb, 1–10 kb, 10–100 kb, >100 kb

Signal correlations

Files:

TAZ_TEAD4_Splot.png

YAP1_TEAD4_Splot.png

YAP-TEAD4_plot1.png

Scatter plots of ChIP-seq signal at shared peaks with regression lines and R²

Peak overlap Venn diagrams

Files:

TEAD4_YAPTAZ.png (YAP/TAZ combined vs TEAD4)

YAP_TAZ.png (YAP vs TAZ)

5. Biological Conclusions

From these reproduced analyses, I confirm key conclusions from Zanconato et al.:

YAP1, TAZ, and TEAD4 show sharp, centered enrichment at shared peak summits.

The majority of shared binding sites fall within active enhancer regions, not promoters.

Most peaks are distal to TSSs (1–10 kb and 10–100 kb), consistent with enhancer targeting.

There is strong correlation of signal intensity between TEAD4 and YAP/TAZ at shared peaks (R² ~ 0.66–0.67).

Venn diagrams show substantial overlap between YAP, TAZ, and TEAD4 peaks, supporting cooperative enhancer regulation.

6. Repository Structure

Suggested layout (what I used):
├── data/
│   ├── bam/
│   ├── bigwig/
│   ├── peaks/
│   └── annotations/
├── figures/
│   ├── heatmap.png
│   ├── pichart_with_percentage.png
│   ├── ...
│   └── YAP1_TEAD4_Splot.png
├── scripts/
│   ├── 01_alignment.sh
│   ├── 02_peakcalling.sh
│   ├── 03_make_bigwig.sh
│   └── 04_R_analysis.R
└── README.md
