# G-Quadruplex × DNA Methylation × Gene Expression in Colorectal Cancer

## Project overview

This project investigates the relationship between **G-quadruplex (G4) DNA structures**
and **differential DNA methylation** in colorectal cancer (CRC), and integrates these
findings with **RNA-seq gene expression data** to identify G4-regulated genes that are
both differentially methylated and differentially expressed in tumour tissue.

### Central biological question

G-quadruplexes are non-canonical DNA secondary structures that form in guanine-rich
sequences. They are enriched at gene promoters and are implicated in transcriptional
regulation. This project asks:

1. Are CpG sites inside G4-forming regions (BG4 ChIP-seq peaks) **more or less likely**
   to be differentially methylated in CRC compared to random expectation?
2. Does G4 methylation status predict changes in **gene expression** at nearby genes?
3. Is the pattern consistent across **two independent methylation datasets** — an
   RRBS sequencing cohort (GSE95656) and an Illumina 450k array cohort (GSE75550/GSE77955)?

---

## Datasets

| ID | Type | Samples | Genome | Description |
|----|------|---------|--------|-------------|
| GSE95656 | RRBS methylation | 20 (10T / 10N paired) | hg19 | Primary CRC vs adjacent normal colon |
| GSE75550 / GSE77955 | Illumina 450k array | 12 (6T / 6N paired) | hg19 | CRC vs adjacent normal |
| GSE95132 | RNA-seq (counts) | 20 (10T / 10N paired) | hg38 | CRC vs adjacent normal — same patients as GSE95656 |
| GSE213402 | RRBS + RNA-seq | CRC vs liver metastasis | hg38→hg19 | Used for comparison only |
| BG4 ChIP-seq (e3 PE, e5 PE) | ChIP-seq peaks | — | hg19 / hg38 | G4 structures captured by BG4 antibody |
| Quadron predictions | Computational | — | hg19 | In silico G4 predictions from the Quadron algorithm |

**T = Tumour, N = adjacent Normal.**

---

## Analysis pipeline — three scripts

The project is split into three R scripts that must be run **in order**. Each script
depends on outputs from the previous one.

```
methylation_pipeline_clean.R
        ↓  produces: myobj_norm.rds, meth_united_strict.rds,
        ↓            meth_diff_strict.rds, hyper/hypo BED files,
        ↓            total_cpg_universe_strict.rds,
        ↓            hyper_bed_new, hypo_bed_new,
        ↓            BG4_hypo/hyper_genes_inPromoter
        ↓
chipseq_array_pipeline_clean.R
        ↓  produces: BG4_e3_pe_peaks_hg19_pure,
        ↓            cpg_BG4_e3_pe_peaks_hg19_pure.bed (→ fed back to Script 1 Fisher test),
        ↓            cpg_quadron_df_fixed.bed,
        ↓            anno_bg4_df (ChIPseeker annotation),
        ↓            bed_ann_GSE77955_hyper/hypo (450k array DMPs)
        ↓
g4_methylation_rnaseq_integration_clean.R
           produces: all figures, enrichment results, final DEG sets
```

> **Important cross-dependency:** The BG4 CpG file
> (`cpg_BG4_e3_pe_peaks_hg19_pure.bed`) is generated in **Script 2 Section 3** but
> consumed in **Script 1 Section 7** (Fisher test). If running from raw data for the
> first time, run Script 2 Section 3 first, then continue with Script 1 from Section 7.
> If the BED file already exists on disk, run the scripts in order as written.

---

## Script 1 — `methylation_pipeline_clean.R`

### Purpose
Full RRBS methylation analysis for GSE95656: from raw Bismark coverage files through
to differentially methylated CpGs (DMCs), statistical enrichment tests, and gene
annotation.

### Section-by-section guide

#### Section 0 — Libraries
Loads all Bioconductor and CRAN packages in a strict order that avoids namespace
conflicts. The load order is: Bioc base → annotation → methylKit → file I/O →
tidyverse → valr → parallel → visualisation.

**Critical conflict fixes applied after loading:**
```r
seqinfo    <- GenomeInfoDb::seqinfo
seqnames   <- GenomicRanges::seqnames
seqlengths <- GenomeInfoDb::seqlengths
start      <- BiocGenerics::start
end        <- BiocGenerics::end
strand     <- BiocGenerics::strand
```
These are necessary because `dplyr`, `data.table`, and `GenomicRanges` all export
functions with the same names. Without explicit assignment, later-loaded packages
silently overwrite the Bioconductor versions and cause cryptic errors.

#### Section 1 — Paths & global parameters
All hardcoded file paths are collected here. When moving the project to a new machine,
only this section needs to be updated.

Key globals:
- `N_PERM = 1000` — number of permutations for the permutation test
- `N_CORES = 7` — parallel cores (adjust to your machine)

#### Section 2 — Load preprocessed methylKit objects
Loads `myobj_norm.rds` — the filtered, normalised methylKit object produced by the
methylKit pipeline described at the bottom of this section. This avoids re-running
`methRead()` (slow) every time the script is opened.

**If running from raw `.cov` files for the first time**, the full pipeline is in the
archived Section at the bottom of the script under the `#### methylKit pipeline` heading.
The key steps are:
1. `methRead()` with `pipeline = "bismarkCoverage"` and `dbtype = "tabix"`
2. `filterByCoverage(lo.count = 10, hi.perc = 99.9)` — removes low and abnormally
   high coverage CpGs (likely PCR duplicates)
3. `normalizeCoverage(method = "median")` — corrects for sequencing depth differences
4. Save as `myobj_norm.rds`

#### Section 3 — Unite samples & QC
`methylKit::unite()` is called with `min.per.group = NULL` which means **all 20 samples**
must have coverage at a CpG for it to be included. This is the **strict intersection**
approach. The alternative (union, i.e. filling missing values) inflates CpG counts and
introduces imputed data, which is avoided here.

Expected output: 500,000–1,500,000 CpGs for 20-sample RRBS.

QC plots produced:
- `PCASamples()` — checks for outlier samples and batch effects
- `clusterSamples(dist = "correlation", method = "ward")` — hierarchical clustering

#### Section 4 — Differential methylation
`calculateDiffMeth()` fits a **logistic regression** at each CpG comparing tumour vs
normal. The model is:
```
methylation ~ treatment (0 = normal, 1 = tumour)
```
Output is a `methylDiff` object saved as `meth_diff_strict.rds`, and also extracted
as a plain data frame `meth_diff` for downstream use with `valr` and `dplyr`.

Columns in `meth_diff`:
| Column | Description |
|--------|-------------|
| chr | Chromosome |
| start | 1-based start position (methylKit convention) |
| end | 1-based end position |
| strand | Strand |
| pvalue | Unadjusted p-value from logistic regression |
| qvalue | BH-adjusted q-value |
| meth.diff | Methylation difference (%) — positive = hypermethylated in tumour |

#### Section 5 — DMC extraction & BED export
DMCs are extracted at four thresholds (5%, 10%, 15%, 25% methylation difference)
all at q < 0.05.

**Coordinate system conversion:** methylKit outputs 1-based coordinates. BED format
requires 0-based start. The `to_bed0based()` function subtracts 1 from `start`:
```r
start = as.integer(df$start) - 1L   # 1-based → 0-based BED
end   = as.integer(df$end)           # end is inclusive in both systems
```

`hyper_bed_new` and `hypo_bed_new` are additionally end-corrected (`end = end + 1`)
so every CpG has width = 1. This is required for `valr::bed_intersect()` to work
correctly with single-base features.

The `total_cpg_df` object (all tested CpGs, 0-based) serves as the **CpG universe**
for the Fisher and permutation tests — it is the denominator.

#### Section 6 — Helper functions
Three functions used throughout Sections 7–9:

**`to_valr(df)`** — converts any data frame with (chrom, start, end) as first three
columns into a valr-compatible tibble. Filters to standard chromosomes only
(`chr1–22, X, Y`) and sorts by coordinate.

**`to_gr(df, zero_based)`** — converts a data frame to a `GRanges` object.
Set `zero_based = TRUE` for BED/NarrowPeak files (subtracts 1 from start back to
1-based for GRanges). Set `zero_based = FALSE` for already-1-based data.

**`genome_df`** — chromosome size table built from BSgenome, needed by
`valr::bed_shuffle()` in the permutation test.

#### Section 7 — Fisher test (covered BG4 CpGs only)
This is the primary enrichment test. It asks: **"Are CpG sites inside G4 peaks more
or less likely to be differentially methylated than CpG sites outside G4 peaks?"**

The key design decision is to use **only BG4 CpGs covered by the RRBS data** as the
background, not all BG4 CpGs. This matters because RRBS has incomplete CpG coverage —
BG4 peaks outside RRBS-covered regions cannot possibly show DMCs, so including them
would deflate the rate inside peaks and bias the OR toward depletion.

Contingency table:

```
                DM       not DM
BG4 covered:    A           B
non-BG4:        C           D
```

- A = BG4-covered CpGs that are DMCs
- B = BG4-covered CpGs that are not DMCs
- C = non-BG4 CpGs that are DMCs
- D = non-BG4 CpGs that are not DMCs

`fisher.test(matrix(c(A,C,B,D), nrow=2), alternative = "two.sided")`

OR < 1 = G4 CpGs are **depleted** for DMCs (protected from methylation change)
OR > 1 = G4 CpGs are **enriched** for DMCs

#### Section 8 — Permutation test (parallel, CpG-constrained)
A permutation test that shuffles BG4 peaks randomly within the RRBS-covered CpG
universe and counts overlaps with DMCs at each shuffle. This tests whether the
**spatial distribution** of G4 peaks relative to DMCs is non-random.

Key design choices:
- **`incl = cpg_universe`** in `valr::bed_shuffle()` restricts shuffled peaks to
  land only within RRBS-covered CpG regions. This is more biologically correct than
  genome-wide shuffling because it matches the null hypothesis to the actual data
  coverage.
- **`parLapply` instead of `mclapply`** — `mclapply` uses forking which is unstable
  on macOS with Bioconductor packages. `parLapply` with `makeCluster()` uses socket
  connections which are safer.
- Workers explicitly load `valr` and `dplyr` via `clusterEvalQ()`.

Output for each comparison:
- `observed` — actual overlap count
- `perm_mean` / `perm_sd` — null distribution summary
- `pval` — two-tailed empirical p-value
- `zscore` — effect size in standard deviation units
- `fold_enrichment` — observed / perm_mean
- `perm_dist` — full vector of permuted counts (for plotting)

#### Section 9 — Visualisation
All plots use a consistent colour scheme across the project:
- Hypermethylated: `#D85A30` (orange-red)
- Hypomethylated: `#378ADD` (blue)
- Expected/background: `#AAAAAA` (grey)

Five plots are produced and saved as PDFs in `output_dir`.

#### Section 10 — Genome annotations
Defines promoter regions (±2kb around TSS from TxDb hg19), TSS points, and
SCREEN enhancers (mock-subtracted to remove promoter overlap).

#### Section 11 — File format conversions
Utilities for BigWig → BED conversion and hg38 → hg19 liftover using
`rtracklayer::liftOver()` and a downloaded chain file.

#### Section 12 — RNA-seq (DESeq2, GSE95132 + GSE213402)
Full DESeq2 analysis for the GSE95132 RNA-seq dataset (CRC vs normal, 10 paired
samples). Also includes a limma analysis on the merged GSE95132 + GSE213402 log2-FPKM
matrix. The DESeq2 object `dds_gse95132` saved here is loaded by Script 3.

---

## Script 2 — `chipseq_array_pipeline_clean.R`

### Purpose
Process BG4 and Quadron G4 peak files, extract CpG dinucleotides within peaks,
run the 450k array methylation analysis (GSE75550/GSE77955), and compute all
pairwise intersections between G4 sets and DMC sets across both datasets.

### Section-by-section guide

#### Section 0 — Libraries
Includes `minfi`, `IlluminaHumanMethylation450k*` packages for 450k array processing,
plus `ChIPseeker` for peak annotation.

#### Section 1 — Paths & global parameters
Sample IDs for the GSE77955 450k dataset are defined as constants `TUMOR_SAMPLES`
and `NORMAL_SAMPLES` to avoid re-typing them throughout the mean methylation
calculations.

#### Section 2 — Load & filter ChIP-seq peak files
Loads four peak sets:
- `BG4_e3_pe_peaks_hg19` — BG4 ChIP-seq, experiment 3, paired-end, hg19
- `BG4_e5_pe_peaks_hg19` — BG4 ChIP-seq, experiment 5, paired-end, hg19
- `mock_e3_pe_peaks_hg19` / `mock_e5_pe_peaks_hg19` — IgG mock controls

Mock subtraction is performed using `valr::bed_subtract()` to remove background
peaks from each BG4 experiment, producing `_pure` variants. This is the correct
approach — the original script erroneously used `bed_intersect()` for this step
(which returns the overlap rather than removing it).

The `_2kbUP` variants extend each peak's start by 2kb upstream to capture promoter
regions proximal to G4 structures. `pmax(0, start - 2000)` prevents negative
coordinates at chromosome starts.

#### Section 3 — CpG dinucleotide extraction from peak intervals
For each peak set, the genomic sequence is retrieved from BSgenome and all "CG"
dinucleotides are located using `Biostrings::matchPattern()`. This converts
broad ChIP-seq peaks into single-base CpG resolution, which is necessary because
methylation data is measured per CpG.

**Coordinate arithmetic:**
```
genomic_start (1-based) = peak_start + match_start - 1
→ convert to 0-based BED: subtract 1 from both start and end
```

The `extract_cpgs_from_gr()` helper encapsulates this logic and is applied to three
peak sets (BG4 e5 raw, BG4 e3 pure, Quadron hg19).

**The output BED file `cpg_BG4_e3_pe_peaks_hg19_pure.bed` is the primary input for
the Fisher and permutation tests in Script 1, Section 7.**

#### Section 4 — Load methylation results (GSE95656)
Reads the hyper and hypo DMC BED files produced by Script 1 Section 5. These are
0-based BED files with columns: chrom, start, end, name, score (= meth.diff), strand.

#### Section 5 — Array methylation (GSE75550 / GSE77955, minfi + limma)

**5A: Load IDAT files** using `minfi::read.metharray.exp()`.

**5B: QC** — Detection p-values are computed for each probe. Samples with > 1%
failed probes (detP > 0.01) are removed, then probes failing in any remaining
sample are filtered out. This is the standard minfi QC protocol.

**5C: Annotation** — `getAnnotation()` retrieves probe coordinates and genomic
context. These are exported as a 0-based BED file for use with valr.

**5D: Normalisation** using `preprocessNoob()` — background correction and dye-bias
normalisation recommended for 450k data.

**5E: Differential methylation** with `limma` on M-values (logit-transformed beta
values). M-values are used for the linear model because they are more statistically
well-behaved (heteroskedasticity) than raw beta values. Delta-beta (Δβ) is added
afterwards as an interpretable effect size on the 0–1 scale.

DMPs are labelled: hyper (Δβ > 0.10, FDR < 5%), hypo (Δβ < −0.10, FDR < 5%).

#### Section 6 — Genomic intersections
All combinations of (G4 set) × (methylation direction) × (dataset) are computed
using `valr::bed_intersect()`. This produces 16 intersection objects covering:
- BG4 e5 / BG4 e5 +2kb / Quadron / Quadron +2kb
- × hyper / hypo
- × GSE95656 (RRBS) / GSE77955 (450k)

#### Section 7 — Fisher test: hyper vs hypo enrichment rate
This Fisher test asks a **different question** from the one in Script 1 Section 7.
Here the question is: **"Does the proportion of peaks overlapping DMCs differ between
hypermethylated and hypomethylated CpGs?"**

Contingency table:
```
          overlap   no_overlap
hyper:       a          b
hypo:        c          d
```
OR > 1 = hyper CpGs are more likely than hypo CpGs to overlap G4 peaks.

#### Section 8 — Mean methylation at G4 CpGs (450k array)
Identifies which 450k probes map to positions inside BG4 / Quadron peaks by
intersecting probe coordinates with CpG-level BED files from Section 3. For each
matching probe, mean beta values are computed separately for tumour and normal groups.

#### Section 9 — ChIPseeker peak annotation
Annotates BG4 (hg38) and Quadron (hg19) peaks relative to genomic features using
the appropriate TxDb for each genome version. Note: BG4 e5 hg38 uses
`TxDb.Hsapiens.UCSC.hg38.knownGene`; Quadron hg19 uses
`TxDb.Hsapiens.UCSC.hg19.knownGene`. Using the wrong TxDb is a common mistake
that produces incorrect annotations without any error message.

---

## Script 3 — `g4_methylation_rnaseq_integration_clean.R`

### Purpose
Integrate all upstream results into figures and final biological conclusions:
overlap visualisations, enrichment statistics with plots, RNA-seq DEG analysis
with paired design, pathway enrichment, and G4-gene expression integration.

### Section-by-section guide

#### Section 0 — Libraries & shared functions
Defines three functions used throughout:

**`to_gr(df, zero_based)`** — as in Script 1. Defined here again because this script
can be loaded independently in a new R session.

**`count_hits(peaks_gr, label, hyper_gr, hypo_gr)`** — uses `countOverlaps()` to
count how many hyper and hypo DMC CpGs fall inside each peak. Returns only peaks
with at least one DM CpG.

**`add_category(df)`** — classifies each peak:
- "Hyper only" — peak contains hyper DMCs but no hypo
- "Hypo only" — peak contains hypo DMCs but no hyper
- "Mixed" — peak contains both hyper and hypo DMCs (mechanistically interesting —
  suggests within-peak heterogeneity)

#### Section 1 — Paths & data setup
Creates aliases so that objects from Scripts 1 and 2 (which use descriptive long
names) map to the shorter names used in this script:
```r
BG4e5PE_hg19_pure      <- BG4_e5_pe_peaks_hg19_pure
quadron_e5PE_peak_hg19 <- quadron_BG4_e5_hg19
```
Then builds `GRanges` objects for all four sets (`bg4`, `quadron`, `hyper`, `hypo`).

#### Section 2 — Per-peak DM CpG counting & categorisation
Produces `bg4_dm` and `quadron_dm` — data frames where each row is a peak that
contains at least one DMC CpG, with columns: `peak_id`, `set`, `hyper` (count),
`hypo` (count), `category`.

These are the source data for all downstream visualisations and the G4 × RNA-seq
integration in Section 15.

#### Section 3 — Euler diagram
Uses the `eulerr` package to draw proportional circles showing set sizes and overlaps.
The key numbers hardcoded here (3867 BG4 peaks, 815 Quadron, 174 hypo, 65 hyper, etc.)
come from the `bg4_dm` / `quadron_dm` analysis in Section 2.

`error_plot(fit)` should be checked before finalising the figure — it shows how well
the circle areas match the actual set sizes. For very skewed size ratios, the fit may
be poor and an UpSet plot (Section 4) is more accurate.

#### Section 4 — UpSet plot
Handles higher-order intersections that Venn/Euler diagrams cannot represent cleanly.
The "Mixed" category (11 peaks with both hyper and hypo CpGs) is only representable
in this format. The `fromExpression()` input counts are:
- All categories are **subsets of BG4** — Quadron-only, HyperDMC-only, HypoDMC-only
  counts are zero by construction and omitted automatically.

#### Section 5 — DM CpGs per peak: violin & stacked bar
Two complementary views of the same data:
- **Violin** (log scale): shows the distribution of DM CpG counts per peak —
  useful for seeing whether peaks tend to have one CpG or many.
- **Stacked bar**: shows the proportion of Hyper only / Hypo only / Mixed peaks —
  useful for comparing BG4 vs Quadron patterns.

#### Section 6 — ChIPseeker annotation × methylation category
Joins the ChIPseeker annotation (from Script 2 Section 9) with the per-peak
methylation category (from Section 2 of this script) using `peak_id` as the key.
Produces a stacked bar showing what proportion of Hyper only / Hypo only / Mixed
peaks fall in promoters, exons, introns, etc.

#### Section 7 — CpG island overlap by methylation category
Loads CpG islands from the UCSC `cpgIslandExt.txt.gz` file (the same file used for
permutation shuffling in Script 1 Section 8) and tests whether hyper/hypo/mixed
peaks differ in their CpG island overlap rate.

#### Section 8 — Fisher test (summary counts)
A third Fisher test, distinct from those in Scripts 1 and 2. This one uses
**pre-summarised counts** (hardcoded from the full analysis) and tests enrichment
at the level of raw CpG counts genome-wide. It is the most straightforward version
to report in a paper methods section.

The three Fisher tests in this project ask:
1. **Script 1 Section 7** — Are G4 CpGs enriched for DMCs vs. non-G4 CpGs?
   (background = RRBS-covered CpGs, split by G4 membership)
2. **Script 2 Section 7** — Is the G4 overlap rate different between hyper and
   hypo DMCs? (tests directionality)
3. **Script 3 Section 8** — Are DMCs enriched inside G4 peaks vs. genome-wide
   expected? (based on peak CpG fraction)

All three ask related but distinct questions and should all be reported together.

#### Section 9 — Permutation test (genome-wide GRanges shuffle)
A faster but less conservative version of the permutation test compared to Script 1
Section 8. Peaks are shuffled anywhere in the genome (not restricted to CpG regions).
Use this for a quick sanity check; use Script 1 Section 8 for the result to report.

The `shuffle_peaks()` function samples chromosomes weighted by their size, then
picks a random start position. This preserves the original peak widths.

#### Sections 10A–C — Enrichment visualisation
Three complementary figures:
- **Permutation distribution** (10A): histogram of null distribution with observed
  value marked — shows how far the observed overlap is from random expectation.
- **Forest plot** (10B): OR with 95% CI on log scale — the standard format for
  reporting Fisher test results.
- **Observed vs expected** (10C): bar chart comparing the raw counts — most
  intuitive for non-statistical audiences.

#### Section 11 — RNA-seq: DESeq2 paired analysis & robust DEG set
Loads the DESeq2 object from Script 1 and **refits with a paired design**:
```r
design(dds) <- ~ patient + condition
```
This is important — the original analysis in Script 1 used `~ condition` only.
Adding `patient` as a blocking factor removes inter-patient variability and
substantially increases power for detecting tumour-specific effects.

The **robust DEG set** strategy:
1. Run DESeq2 on all 20 samples → DEG set A
2. Remove 3 outlier patients (41, 42, 51) identified by PCA → DEG set B
3. Robust DEGs = intersection of A and B

This approach is conservative but reduces false positives driven by outlier samples.
The outlier patients were identified visually from the PCA plot (Section 12) — they
cluster with the wrong condition, suggesting sample labelling uncertainty or unusual
biology.

LFC shrinkage is applied using `lfcShrink(type = "apeglm")` which reduces noise in
log2 fold-change estimates for low-count genes. This is the recommended approach for
visualisation and ranking, though it does not affect the p-values.

#### Section 12 — RNA-seq QC
PCA with patient labels to identify outliers. `blind = TRUE` is used for
visualisation (correct — ignores the design when computing variance). `blind = FALSE`
would use the design during VST, which is appropriate for downstream analysis but
distorts the QC view.

#### Section 13 — Volcano, MA, heatmap
Standard DEG visualisation. Three layers on the volcano plot:
- Grey points: non-significant genes
- Faded coloured points: significant in full dataset but not robust
- Bright coloured points + labels: robust DEGs (significant in both analyses)

This layered approach makes the uncertainty in the DEG set visually explicit.

#### Section 14 — GO / KEGG enrichment
Uses `clusterProfiler::enrichGO()` and `enrichKEGG()` on the robust DEG set only.
A `universe` parameter is passed (all tested genes, not just DEGs) to avoid
enrichment inflation from small universe size. Results are visualised as dotplots
showing enrichment ratio and adjusted p-value.

#### Section 15 — G4 × RNA-seq integration
Joins the per-peak methylation category (Section 2) with the ChIPseeker gene
annotation (Script 2 Section 9), then intersects with robust DEGs to identify
**genes that are both differentially expressed AND located near a G4 peak that
is differentially methylated**.

Two figures:
- **Violin plot**: log2FC distribution for genes near Hyper only / Hypo only /
  Mixed G4 peaks — tests whether G4 methylation direction predicts expression
  direction.
- **G4-coloured volcano**: all genes shown, but those near DM G4 peaks are
  coloured and labelled. This is the key integration figure.

#### Section 16 — BG4 promoter gene heatmap
Shows expression of genes whose promoters contain BG4-covered DMC CpGs.
`BG4_hypo_genes_inPromoter` and `BG4_hyper_genes_inPromoter` are produced at the
end of Script 1 using `findOverlaps()` + `mapIds()`.

The heatmap uses manual ordering: rows are clustered within each G4 group (hypo /
both / hyper) separately, then concatenated. Columns are clustered within
Normal and Tumour groups separately. `gaps_row` and `gaps_col` draw visual
separators between the groups.

---

## Colour palette (used consistently across all scripts)

| Colour | Hex | Usage |
|--------|-----|-------|
| Orange-red | `#D85A30` | Hypermethylated, Tumour, Upregulated, Observed |
| Blue | `#378ADD` | Hypomethylated, Normal, Downregulated, Expected |
| Purple | `#7F77DD` | Hypomethylated (alternative), Hypo only peaks |
| Teal | `#1D9E75` | Quadron peaks, CpG islands |
| Maroon | `#B4506A` | Mixed methylation category |
| Grey | `#AAAAAA` | Background, non-significant, expected |
| Dark grey | `#444441` | UpSet bar/matrix elements |
| Off-white | `#F1EFE8` | UpSet row shading |

---

## Key statistical decisions and rationale

### Why strict intersection in methylKit::unite()?
Using `min.per.group = NULL` requires all 20 samples to have RRBS coverage at a
CpG. This produces fewer CpGs (~500K–1.5M) but avoids imputing missing values,
which would inflate DMC counts and reduce the stringency of the analysis.

### Why use covered BG4 CpGs as the Fisher denominator?
The standard approach would compare DMC rate inside vs outside BG4 peaks using
all BG4 CpGs as the denominator. However, RRBS only covers ~10–20% of all CpGs
genome-wide. BG4 CpGs not covered by RRBS cannot show up as DMCs by definition.
Including them in the denominator makes the inside-BG4 DMC rate appear artificially
low, biasing the odds ratio toward depletion. Using only the RRBS-covered BG4 CpGs
as the denominator gives a fair comparison.

### Why shuffle within the CpG universe for permutations?
When shuffling BG4 peaks to build a null distribution, shuffling genome-wide allows
peaks to land in heterochromatin, repeat elements, or other regions that RRBS never
covers. A shuffled peak in those regions will trivially have zero DMC overlaps,
making the null distribution unrealistically low and inflating significance.
Restricting shuffling to RRBS-covered CpG regions (`incl = cpg_universe`) means
the null distribution reflects what we would expect if G4 peaks were distributed
randomly among the regions we actually measured.

### Why the paired DESeq2 design matters for RNA-seq?
The 10 tumour/normal pairs in GSE95132 come from the same patient. Without blocking
on patient, the model conflates within-patient (tumour vs normal) effects with
between-patient variability. This reduces power and can produce false positives if
some patients have unusually high or low baseline expression. The paired design
`~ patient + condition` removes between-patient variability before testing
the condition effect.

### Why define a robust DEG set?
Three patients (41, 42, 51) were identified as outliers by PCA — their samples
cluster with the wrong condition. This could reflect sample labelling errors,
unusual tumour biology, or batch effects. Rather than simply removing them (which
would be a post-hoc decision that inflates the DEG count), the robust DEG strategy
keeps them in the primary analysis but requires each DEG to be significant in both
the full and the cleaned dataset. This is conservative and transparent.

---

## Output files

| File | Produced by | Description |
|------|-------------|-------------|
| `myobj_norm.rds` | Script 1 Section 2 | Filtered, normalised methylKit object |
| `meth_united_strict.rds` | Script 1 Section 3 | United methylBase object (strict intersection) |
| `meth_diff_strict.rds` | Script 1 Section 4 | methylDiff object |
| `hyper_GSE95656_strict_5.bed` | Script 1 Section 5 | Hyper DMCs at 5% threshold (0-based BED) |
| `hypo_GSE95656_strict_5.bed` | Script 1 Section 5 | Hypo DMCs at 5% threshold (0-based BED) |
| `total_cpg_universe_strict.rds` | Script 1 Section 5 | All RRBS-covered CpGs (CpG universe) |
| `cpg_BG4_e3_pe_peaks_hg19_pure.bed` | Script 2 Section 3 | CpG positions inside BG4 e3 peaks |
| `cpg_quadron_df_fixed.bed` | Script 2 Section 3 | CpG positions inside Quadron peaks |
| `ResLimma_GSE77955_meth.csv` | Script 2 Section 5 | 450k array limma results |
| `dds_results_gse95132.rds` | Script 1 Section 12 | DESeq2 object for RNA-seq |
| `plot_permutation_distribution.pdf` | Script 1 Section 9 | Permutation null distributions |
| `plot_DM_rate_inside_outside_G4.pdf` | Script 1 Section 9 | DM rate inside/outside G4 |
| `dm_cpgs_per_peak.pdf` | Script 3 Section 5 | DM CpGs per peak violin |
| `permutation_test_BG4.pdf` | Script 3 Section 10 | Permutation distribution (genome-wide) |
| `fisher_test_plot.pdf` | Script 3 Section 10 | Forest plot + obs vs expected |
| `RNAseq_G4_integration.pdf` | Script 3 Section 15 | All RNA-seq + G4 figures |

---

## R package requirements

```r
# Bioconductor
BiocManager::install(c(
  "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges",
  "Biostrings", "XVector", "AnnotationDbi", "org.Hs.eg.db", "GenomicFeatures",
  "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "BSgenome", "BSgenome.Hsapiens.UCSC.hg19",
  "methylKit", "Biobase", "rtracklayer", "BiocIO", "Rsamtools",
  "minfi", "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "ChIPseeker", "clusterProfiler", "enrichplot",
  "DESeq2", "apeglm", "limma"
))

# CRAN
install.packages(c(
  "data.table", "readr", "dplyr", "tidyr",
  "valr", "parallel",
  "ggplot2", "ggrepel", "pheatmap", "patchwork", "scales",
  "UpSetR", "eulerr",
  "tibble"
))
```

---

## File paths that need updating

All hardcoded paths are collected in Section 1 of each script. Search for
`/Users/hossein.allahdadi/` and replace with your own base directory.
The paths you need to update are:

| Variable | Script | Points to |
|----------|--------|-----------|
| `rds_dir` | Script 1 | Folder containing `.rds` intermediate files |
| `output_dir` | Script 1, 3 | Folder for output BED files and plots |
| `bgz_dir` | Script 1 | methylKit tabix database folder |
| `bg4_cpg_path` | Script 1 | `cpg_BG4_e3_pe_peaks_hg19_pure.bed` |
| `chip_dir` | Script 2 | Folder with NarrowPeak files |
| `meth_dir` | Script 2 | GSE75550 IDAT folder |
| `rds_dds` | Script 3 | Path to `dds_results_gse95132.rds` |
| Chain file path | Script 1 Section 11 | `hg38ToHg19.over.chain` |
| CpG islands path | Scripts 1, 3 | `cpgIslandExt.txt.gz` |
