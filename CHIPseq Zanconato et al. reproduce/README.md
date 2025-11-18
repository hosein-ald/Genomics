# Reproducible ChIP-seq Analysis of YAP1, TAZ, and TEAD4

### Reproducing Figures from *Zanconato et al., Cell (2018)* Using CrazyHotTommy's Pipeline

## 📌 Project Overview

This repository documents my complete attempt to reproduce **Figure 1
and related analyses** from:

**Zanconato et al. (2018)**\
*YAP/TAZ--TEAD and AP-1 Cooperate to Control Oncogenic Enhancers in
Cancer Cells*\
https://pmc.ncbi.nlm.nih.gov/articles/PMC6186417/

Workflow based on:\
**CrazyHotTommy --- Reproduce Genomics Figures**\
https://crazyhottommy.github.io/reproduce_genomics_paper_figures/index.html

------------------------------------------------------------------------

## 📂 Data Used

-   YAP1 ChIP-seq\
-   TAZ ChIP-seq\
-   TEAD4 ChIP-seq\
-   IgG control

Processed into `.bam`, `.bw`, `.narrowPeak` files.

**MACS3 peak call example:**

``` bash
macs3 callpeak -t <IP>.sorted.bam   -c IgG.sorted.bam   -f BAM   -n <sample_name>
```

------------------------------------------------------------------------

## ⚙️ Pipeline Summary

### Alignment

``` bash
bowtie2 -x hg38   -1 sample_R1.fq.gz -2 sample_R2.fq.gz   | samtools sort -o sample.sorted.bam
samtools index sample.sorted.bam
```

### Peak Calling

``` bash
macs3 callpeak -t sample.sorted.bam   -c IgG.sorted.bam   -n sample   -f BAM
```

### BigWig Generation

``` bash
bamCoverage -b sample.sorted.bam   -o sample.bw   --normalizeUsing RPKM
```

### Heatmaps (deepTools)

``` bash
computeMatrix reference-point   -S YAP.bw TAZ.bw TEAD4.bw   -R peaks.bed   --referencePoint center   -b 3000 -a 3000   -o matrix.gz

plotHeatmap -m matrix.gz -o figures/heatmap.png
```

------------------------------------------------------------------------

## 📊 Reproduced Figures

(Stored in `figures/`)

-   heatmap.png\
-   pichart_with_percentage.png\
-   Screenshot_2025-11-17.png\
-   distance_to_the_summit_of_TAZ_peaks.png\
-   Distance_to_TSS.png\
-   TAZ_TEAD4_Splot.png\
-   TEAD4_YAPTAZ.png\
-   YAP_TAZ.png\
-   YAP-TEAD4_plot1.png\
-   YAP1_TEAD4_Splot.png

------------------------------------------------------------------------

## 📁 Repository Structure

``` text
├── data/
├── figures/
├── notebooks/
│   └── YAP_TAZ_TEAD4_reproduction.ipynb
├── scripts/
│   ├── 01_alignment.sh
│   ├── 02_peakcalling.sh
│   ├── 03_make_bigwig.sh
│   └── 04_R_analysis.R
└── README.md
```

------------------------------------------------------------------------

## 📝 How to Run

``` bash
git clone <repo_url>
cd <repo>
```

Create environment:

``` bash
conda create -n chipseq_env python=3.10 bowtie2 samtools macs3 deeptools
conda activate chipseq_env
```

Run workflow:

``` bash
bash scripts/01_alignment.sh
bash scripts/02_peakcalling.sh
bash scripts/03_make_bigwig.sh
Rscript scripts/04_R_analysis.R
```

------------------------------------------------------------------------

## 📚 References

-   Zanconato et al. (2018). *YAP/TAZ--TEAD and AP-1 Cooperate to
    Control Oncogenic Enhancers in Cancer Cells.*\
-   CrazyHotTommy --- *Reproduce Genomics Paper Figures*.
