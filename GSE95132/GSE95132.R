# ══════════════════════════════════════════════════════════════════════════════
# METHYLATION ANALYSIS PIPELINE
# Study:   G-quadruplex (BG4) enrichment at differentially methylated CpGs
# Dataset: GSE95656 — RRBS colorectal cancer (20 samples, hg19)
# Author:  Hossein Allahdadi
#
# Pipeline overview:
#   0.  Libraries
#   1.  Paths & global parameters
#   2.  Load preprocessed methylKit objects
#   3.  Unite samples & QC
#   4.  Differential methylation (calculateDiffMeth)
#   5.  DMC extraction at multiple thresholds & BED export
#   6.  Helper functions  (shared across sections 7–9)
#   7.  Fisher test       (covered BG4 CpGs only — honest denominator)
#   8.  Permutation test  (parallel, shuffled within CpG universe)
#   9.  Visualisation
#  10.  Genome annotations (TSS, Enhancers, CpG Islands)
#  11.  File format conversions (BigWig → BED, LiftOver hg38 → hg19)
#  12.  RNA-seq analysis (GSE95132 + GSE213402 — DESeq2 / limma)
#  ──  ARCHIVE: older / exploratory code (commented out)
# ══════════════════════════════════════════════════════════════════════════════


# ── SECTION 0: LIBRARIES ──────────────────────────────────────────────────────
# Load order matters: Bioc base → annotation → analysis → tidyverse → viz
# Conflict fixes are applied at the end of this block.

# Block 1 — Bioconductor base infrastructure
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)      # seqinfo(), seqnames(), seqlengths()
library(GenomicRanges)     # GRanges(), IRanges()
library(Biostrings)        # matchPattern(), getSeq()
library(XVector)

# Block 2 — Genome annotation
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

# Block 3 — Methylation analysis
library(methylKit)
library(Biobase)

# Block 4 — File I/O and genomic formats
library(rtracklayer)       # import(), export(), import.chain()
library(BiocIO)
library(Rsamtools)

# Block 5 — Data wrangling (load AFTER Bioc to avoid masking)
library(data.table)        # fread()
library(readr)             # read_table(), read_tsv()
library(dplyr)             # filter(), mutate(), distinct()
library(tidyr)             # pivot_longer()

# Block 6 — Genomic interval operations (load AFTER dplyr — pipe conflict)
library(valr)              # bed_intersect(), bed_shuffle()

# Block 7 — Parallel computing
library(parallel)

# Block 8 — Visualisation (always last)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# Conflict fixes — explicitly assign the correct namespace versions
seqinfo    <- GenomeInfoDb::seqinfo
seqnames   <- GenomicRanges::seqnames
seqlengths <- GenomeInfoDb::seqlengths
start      <- BiocGenerics::start
end        <- BiocGenerics::end
strand     <- BiocGenerics::strand

# Quick sanity check
test_gr <- GRanges("chr1", IRanges(1, 100))
cat("seqnames test:", as.character(GenomicRanges::seqnames(test_gr)), "\n")
cat("seqinfo  test:", class(GenomicRanges::seqinfo(test_gr)), "\n")
cat("All packages loaded successfully\n")

# DESeq2 — only load when needed (heavy namespace conflicts with methylKit)
# library(DESeq2)
# library(limma)


# ── SECTION 1: PATHS & GLOBAL PARAMETERS ─────────────────────────────────────
# Change these paths before running on a new machine.

rds_dir    <- "/Users/hossein.allahdadi/Downloads/ChIP/GSE95656/meth_files/GSE95656_RAW/GSE95565 results"
output_dir <- "/Users/hossein.allahdadi/Downloads/ChIP/GSE95656/meth_files/GSE95656_RAW/GSE95565 results"
bed_dir    <- "/Users/hossein.allahdadi/Downloads/ChIP/GSE95656/meth_files/GSE95656_RAW/result bed files"
bgz_dir    <- "/Users/hossein.allahdadi/methylDB_2026-03-12_oVU"

# BG4 ChIP-seq CpG file (CpG dinucleotides inside BG4 peaks, hg19)
bg4_cpg_path <- "/Users/hossein.allahdadi/Downloads/ChIP/files/cpg_BG4_e3_pe_peaks_hg19_pure.bed"

# Permutation settings
N_PERM  <- 1000   # increase from 100 once results look stable
N_CORES <- 7


# ── SECTION 2: LOAD PREPROCESSED METHYLKIT OBJECTS ───────────────────────────
# These RDS files are produced by the methylKit pipeline (see Section 12 for
# the full methRead → filter → normalize steps if you need to re-run from raw).

myobj_norm <- readRDS(file.path(rds_dir, "myobj_norm.rds"))

cat("Samples loaded:", length(myobj_norm), "\n")
cat("Sample IDs:\n")
print(sapply(seq_along(myobj_norm), function(i) myobj_norm[[i]]@sample.id))
cat("CpGs per sample:\n")
print(sapply(seq_along(myobj_norm), function(i) nrow(methylKit::getData(myobj_norm[[i]]))))


# ── SECTION 3: UNITE SAMPLES & QC ────────────────────────────────────────────
# min.per.group = NULL → ALL samples must have coverage (strict intersection).
# This avoids inflating CpG counts with union-based merging.

meth_strict <- methylKit::unite(myobj_norm,
                                destrand      = FALSE,
                                min.per.group = NULL)

cat("CpGs after strict intersection:", nrow(methylKit::getData(meth_strict)), "\n")
# Expected for RRBS with 20 samples: 500K – 1.5M

saveRDS(meth_strict, file.path(rds_dir, "meth_united_strict.rds"))
cat("Saved: meth_united_strict.rds\n")

# QC plots
PCASamples(meth_strict)
clusterSamples(meth_strict, dist = "correlation", method = "ward", plot = TRUE)


# ── SECTION 4: DIFFERENTIAL METHYLATION ──────────────────────────────────────
# calculateDiffMeth uses logistic regression by default.
# Results are saved as both RDS (methylDiff object) and a plain data frame.

meth_diff_obj <- calculateDiffMeth(meth_strict)
cat("Total CpGs tested:", nrow(methylKit::getData(meth_diff_obj)), "\n")
saveRDS(meth_diff_obj, file.path(rds_dir, "meth_diff_strict.rds"))

# Extract as a plain data frame for downstream use with valr/dplyr
meth_diff <- methylKit::getData(meth_diff_obj)
colnames(meth_diff) <- c("chr", "start", "end", "strand",
                         "pvalue", "qvalue", "meth.diff")


# ── SECTION 5: DMC EXTRACTION & BED EXPORT ───────────────────────────────────
# Extract differentially methylated CpGs (DMCs) at several thresholds.
# Positive meth.diff = hypermethylated in tumour; negative = hypomethylated.

meth_hyper_5  <- meth_diff[meth_diff$meth.diff >  5  & meth_diff$qvalue < 0.05, ]
meth_hypo_5   <- meth_diff[meth_diff$meth.diff < -5  & meth_diff$qvalue < 0.05, ]
meth_hyper_10 <- meth_diff[meth_diff$meth.diff >  10 & meth_diff$qvalue < 0.05, ]
meth_hypo_10  <- meth_diff[meth_diff$meth.diff < -10 & meth_diff$qvalue < 0.05, ]
meth_hyper_15 <- meth_diff[meth_diff$meth.diff >  15 & meth_diff$qvalue < 0.05, ]
meth_hypo_15  <- meth_diff[meth_diff$meth.diff < -15 & meth_diff$qvalue < 0.05, ]
meth_hyper_25 <- meth_diff[meth_diff$meth.diff >  25 & meth_diff$qvalue < 0.05, ]
meth_hypo_25  <- meth_diff[meth_diff$meth.diff < -25 & meth_diff$qvalue < 0.05, ]

cat("─── DMC counts ───────────────────────────────\n")
cat("Hyper_5: ",  nrow(meth_hyper_5),  " | Hypo_5: ",  nrow(meth_hypo_5),  "\n")
cat("Hyper_10:",  nrow(meth_hyper_10), " | Hypo_10:",  nrow(meth_hypo_10), "\n")
cat("Hyper_15:",  nrow(meth_hyper_15), " | Hypo_15:",  nrow(meth_hypo_15), "\n")
cat("Hyper_25:",  nrow(meth_hyper_25), " | Hypo_25:",  nrow(meth_hypo_25), "\n")

# Convert methylKit 1-based coordinates to 0-based BED format
to_bed0based <- function(df, label) {
  data.frame(
    chrom     = as.character(df$chr),
    start     = as.integer(df$start) - 1L,   # 1-based → 0-based
    end       = as.integer(df$end),
    name      = paste0(label, "_", seq_len(nrow(df))),
    score     = round(df$meth.diff, 3),
    strand    = as.character(df$strand),
    pvalue    = df$pvalue,
    qvalue    = df$qvalue,
    meth.diff = df$meth.diff,
    stringsAsFactors = FALSE
  )
}

# Save 5% threshold BED files (primary output)
write.table(to_bed0based(meth_hyper_5, "hyper"),
            file.path(bed_dir, "hyper_GSE95656_strict_5.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(to_bed0based(meth_hypo_5, "hypo"),
            file.path(bed_dir, "hypo_GSE95656_strict_5.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Build the total CpG universe (all tested CpGs, 0-based, standard chrs only)
# Used as the background for Fisher and permutation tests.
total_cpg_df <- data.frame(
  chrom = as.character(meth_diff$chr),
  start = as.integer(meth_diff$start) - 1L,
  end   = as.integer(meth_diff$end),
  stringsAsFactors = FALSE
)
total_cpg_df <- total_cpg_df[grep("^chr([0-9]+|X|Y)$", total_cpg_df$chrom), ]
total_cpg_df <- total_cpg_df[order(total_cpg_df$chrom, total_cpg_df$start), ]
cat("Total CpG universe:", nrow(total_cpg_df), "\n")

saveRDS(total_cpg_df, file.path(rds_dir, "total_cpg_universe_strict.rds"))


# ── SECTION 6: HELPER FUNCTIONS ───────────────────────────────────────────────
# Shared utilities used in Sections 7, 8, and 9.

# to_valr: convert any data frame with (chrom, start, end) as first 3 columns
# into a valr-compatible table (standard chromosomes only, sorted).
to_valr <- function(df) {
  df <- data.frame(
    chrom = as.character(df[[1]]),
    start = as.integer(df[[2]]),
    end   = as.integer(df[[3]]),
    stringsAsFactors = FALSE
  )
  df <- df[grep("^chr([0-9]+|X|Y)$", df$chrom), ]
  df[order(df$chrom, df$start), ]
}

# to_gr: convert a data frame to a GRanges object.
# Set zero_based = TRUE if the input uses 0-based BED start coordinates.
to_gr <- function(df, zero_based = FALSE) {
  GRanges(
    seqnames = df[[1]],
    ranges   = IRanges(start = df[[2]] + as.integer(zero_based),
                       end   = df[[3]])
  )
}

# Chromosome sizes for valr::bed_shuffle()
chrom_sizes <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19)
chrom_sizes <- chrom_sizes[grep("^chr([0-9]+|X|Y)$", names(chrom_sizes))]
genome_df   <- data.frame(
  chrom = names(chrom_sizes),
  size  = as.integer(chrom_sizes),
  stringsAsFactors = FALSE
)


# ── SECTION 7: FISHER TEST (covered BG4 CpGs only) ───────────────────────────
# We restrict the analysis to BG4 CpGs that are actually covered in the RRBS
# data (i.e., present in total_cpg_df). This gives an honest denominator and
# avoids inflating significance with uncovered sites.

# Load BG4 CpG file
cpg_df_e3PE_pure <- read.table(bg4_cpg_path)
colnames(cpg_df_e3PE_pure) <- c("chrom", "start", "end")

bg4_valr   <- to_valr(cpg_df_e3PE_pure)
total_valr <- to_valr(total_cpg_df)

# Choose which DMC threshold to use for Fisher (10% used here — adjust as needed)
hyper_valr <- to_valr(meth_hyper_10)
hypo_valr  <- to_valr(meth_hypo_10)

# Subset BG4 to CpGs covered in RRBS
bg4_covered <- bed_intersect(bg4_valr, total_valr, min_overlap = 0L) %>%
  dplyr::distinct(chrom, start.x, end.x) %>%
  dplyr::transmute(chrom = as.character(chrom),
                   start = as.integer(start.x),
                   end   = as.integer(end.x))

non_bg4_covered <- valr::bed_subtract(total_valr, bg4_valr)

cat("BG4 CpGs total:          ", nrow(bg4_valr),         "\n")
cat("BG4 CpGs covered in RRBS:", nrow(bg4_covered),      "\n")
cat("Coverage rate:           ",
    round(nrow(bg4_covered) / nrow(bg4_valr) * 100, 2),  "%\n")
cat("Non-BG4 covered CpGs:   ", nrow(non_bg4_covered),   "\n")

# Fisher test function — uses covered BG4 CpGs as the true background
run_fisher_covered <- function(dmc_valr, label) {

  A <- bed_intersect(bg4_covered, dmc_valr, min_overlap = 0L) %>%
    dplyr::distinct(chrom, start.x, end.x) %>% nrow()
  B <- nrow(bg4_covered) - A

  C <- bed_intersect(non_bg4_covered, dmc_valr, min_overlap = 0L) %>%
    dplyr::distinct(chrom, start.x, end.x) %>% nrow()
  D <- nrow(non_bg4_covered) - C

  cat("\n===", label, "— covered BG4 only ===\n")
  cat("  A BG4 + DM:         ", A, "\n")
  cat("  B BG4 + not DM:     ", B, "\n")
  cat("  C non-BG4 + DM:     ", C, "\n")
  cat("  D non-BG4 + not DM: ", D, "\n")
  cat("  DM rate inside  BG4:", round(A / nrow(bg4_covered)      * 100, 3), "%\n")
  cat("  DM rate outside BG4:", round(C / nrow(non_bg4_covered)  * 100, 3), "%\n")

  ft <- fisher.test(matrix(c(A, C, B, D), nrow = 2), alternative = "two.sided")

  cat("  OR:", round(ft$estimate, 4),
      "| CI [", round(ft$conf.int[1], 4), "-", round(ft$conf.int[2], 4), "]",
      "| p:", ft$p.value,
      "| direction:", ifelse(ft$estimate > 1, "ENRICHED", "DEPLETED"), "\n")

  data.frame(
    label           = label,
    A_bg4_DM        = A,
    B_bg4_noDM      = B,
    C_nonbg4_DM     = C,
    D_nonbg4_noDM   = D,
    bg4_covered     = nrow(bg4_covered),
    DM_rate_inside  = round(A / nrow(bg4_covered)     * 100, 3),
    DM_rate_outside = round(C / nrow(non_bg4_covered) * 100, 3),
    OR              = round(ft$estimate,    4),
    CI_low          = round(ft$conf.int[1], 4),
    CI_high         = round(ft$conf.int[2], 4),
    pval            = ft$p.value,
    direction       = ifelse(ft$estimate > 1, "ENRICHED", "DEPLETED"),
    stringsAsFactors = FALSE
  )
}

results_covered <- rbind(
  run_fisher_covered(hyper_valr, "Hyper"),
  run_fisher_covered(hypo_valr,  "Hypo")
)
print(results_covered)


# ── SECTION 8: PERMUTATION TEST ──────────────────────────────────────────────
# BG4 peaks are shuffled N_PERM times within the RRBS-covered CpG universe
# (total_cpg_df) to obtain a null distribution of overlaps.
# Uses parLapply (more stable than mclapply on macOS with fork-safety issues).

# CpG universe for shuffling = all RRBS-covered CpGs (0-based)
cpg_universe <- total_cpg_df   # already 0-based from Section 5

cat("CpG universe for permutation:", nrow(cpg_universe), "sites\n")

# Core permutation function
permutation_test <- function(query_valr, subject_valr,
                             cpg_universe, genome_df,
                             n_perm = N_PERM, n_cores = N_CORES,
                             seed = 42) {
  set.seed(seed)
  cat("  Query rows:  ", nrow(query_valr),   "\n")
  cat("  Subject rows:", nrow(subject_valr), "\n")

  obs_hits <- valr::bed_intersect(query_valr, subject_valr, min_overlap = 0L) %>%
    dplyr::distinct(chrom, start.x, end.x) %>%
    nrow()
  cat("  Observed hits:", obs_hits, "\n")

  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl,
                          varlist = c("query_valr", "subject_valr",
                                      "cpg_universe", "genome_df"),
                          envir = environment())
  parallel::clusterEvalQ(cl, { library(valr); library(dplyr) })

  perm_hits <- unlist(parallel::parLapply(cl, seq_len(n_perm), function(i) {
    shuffled <- valr::bed_shuffle(query_valr, genome_df,
                                  incl = cpg_universe, seed = i)
    valr::bed_intersect(shuffled, subject_valr, min_overlap = 0L) %>%
      dplyr::distinct(chrom, start.x, end.x) %>%
      nrow()
  }))
  parallel::stopCluster(cl)

  cat("  Permutations done:", length(perm_hits), "\n")
  cat("  Perm mean:", round(mean(perm_hits), 1), "\n")

  pval_enr  <- sum(perm_hits >= obs_hits) / n_perm
  pval_dep  <- sum(perm_hits <= obs_hits) / n_perm
  pval      <- min(2 * min(pval_enr, pval_dep), 1)
  zscore    <- (obs_hits - mean(perm_hits)) / sd(perm_hits)
  fold      <- obs_hits / mean(perm_hits)

  list(
    observed        = obs_hits,
    perm_mean       = round(mean(perm_hits), 1),
    perm_sd         = round(sd(perm_hits),   1),
    pval            = pval,
    zscore          = round(zscore, 3),
    fold_enrichment = round(fold,   3),
    direction       = ifelse(fold >= 1, "ENRICHED", "DEPLETED"),
    perm_dist       = perm_hits
  )
}

# Run BG4 vs Hyper / Hypo
cat("\n=== Hyper DMCs vs BG4 ===\n")
res_hyper <- permutation_test(bg4_valr, hyper_valr, cpg_universe, genome_df)

cat("\n=== Hypo DMCs vs BG4 ===\n")
res_hypo  <- permutation_test(bg4_valr, hypo_valr,  cpg_universe, genome_df)

# Summarise results
results_perm <- data.frame(
  comparison      = c("Hyper vs BG4", "Hypo vs BG4"),
  observed        = c(res_hyper$observed,        res_hypo$observed),
  perm_mean       = c(res_hyper$perm_mean,        res_hypo$perm_mean),
  fold_enrichment = c(res_hyper$fold_enrichment,  res_hypo$fold_enrichment),
  zscore          = c(res_hyper$zscore,           res_hypo$zscore),
  pval            = c(res_hyper$pval,             res_hypo$pval),
  direction       = c(res_hyper$direction,        res_hypo$direction),
  significant     = c(ifelse(res_hyper$pval < 0.05, "YES", "NO"),
                      ifelse(res_hypo$pval  < 0.05, "YES", "NO")),
  stringsAsFactors = FALSE
)
print(results_perm)


# ── SECTION 9: VISUALISATION ──────────────────────────────────────────────────
# All plots follow the same colour scheme:
#   Hypermethylated = #D85A30 (orange-red)
#   Hypomethylated  = #378ADD (blue)

# ── Plot 1: Permutation distribution ─────────────────────────────────────────
perm_df <- data.frame(
  hits = c(res_hyper$perm_dist, res_hypo$perm_dist),
  type = rep(c("Hypermethylated", "Hypomethylated"),
             each = length(res_hyper$perm_dist))
)
obs_df <- data.frame(
  hits  = c(res_hyper$observed, res_hypo$observed),
  type  = c("Hypermethylated", "Hypomethylated"),
  label = c(paste0("Observed = ", res_hyper$observed),
            paste0("Observed = ", res_hypo$observed))
)

ggplot(perm_df, aes(x = hits, fill = type)) +
  geom_histogram(bins = 50, alpha = 0.8, color = NA) +
  geom_vline(data = obs_df, aes(xintercept = hits),
             linetype = "dashed", linewidth = 1, color = "black") +
  geom_text(data = obs_df,
            aes(x = hits, y = Inf, label = label),
            hjust = 1.1, vjust = 2, size = 3.5,
            color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Hypermethylated" = "#D85A30",
                               "Hypomethylated"  = "#378ADD")) +
  facet_wrap(~type, scales = "free") +
  labs(x = "Permuted overlap count", y = "Frequency") +
  theme_minimal(base_size = 13) +
  theme(legend.position  = "none",
        panel.grid.minor = element_blank(),
        strip.text       = element_text(size = 12))

ggsave(file.path(output_dir, "plot_permutation_distribution.pdf"),
       width = 10, height = 5)

# ── Plot 2: Observed vs Expected bar chart ────────────────────────────────────
bar_df <- data.frame(
  type     = rep(c("Hypermethylated", "Hypomethylated"), each = 2),
  category = rep(c("Observed", "Expected"), 2),
  value    = c(res_hyper$observed, res_hyper$perm_mean,
               res_hypo$observed,  res_hypo$perm_mean)
)

ggplot(bar_df, aes(x = type, y = value, fill = category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("Observed" = "#D85A30", "Expected" = "#AAAAAA")) +
  labs(x = NULL, y = "BG4 sites overlapping DMCs", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(output_dir, "plot_observed_vs_expected.pdf"),
       width = 6, height = 5)

# ── Plot 3: Fold enrichment bar chart ────────────────────────────────────────
ggplot(results_perm, aes(x = comparison, y = fold_enrichment, fill = direction)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_text(aes(label = round(fold_enrichment, 3)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("DEPLETED" = "#378ADD", "ENRICHED" = "#D85A30")) +
  labs(x = NULL, y = "Fold enrichment over random", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(output_dir, "plot_fold_enrichment.pdf"), width = 6, height = 5)

# ── Plot 4: DM rate inside vs outside BG4 (from Fisher results) ───────────────
dm_rate_df <- data.frame(
  type     = rep(c("Hypermethylated", "Hypomethylated"), each = 2),
  location = rep(c("Inside G4", "Outside G4"), 2),
  rate     = c(results_covered$DM_rate_inside[1],  results_covered$DM_rate_outside[1],
               results_covered$DM_rate_inside[2],  results_covered$DM_rate_outside[2])
)

label_df <- data.frame(
  type     = c("Hypermethylated", "Hypomethylated"),
  location = "Outside G4",
  rate     = c(results_covered$DM_rate_outside[1], results_covered$DM_rate_outside[2]),
  label    = paste0("OR = ", results_covered$OR, "\np < 0.0001")
)

ggplot(dm_rate_df, aes(x = location, y = rate, fill = location)) +
  geom_bar(stat = "identity", width = 0.55) +
  geom_text(aes(label = paste0(rate, "%")), vjust = -0.5, size = 4, fontface = "bold") +
  geom_text(data = label_df,
            aes(x = location, y = rate, label = label),
            vjust = -1.8, size = 3.2, color = "gray30", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Inside G4" = "#378ADD", "Outside G4" = "#AAAAAA")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  facet_wrap(~type) +
  labs(x = NULL, y = "Differential methylation rate (%)", fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "top",
        strip.text       = element_text(size = 13, face = "bold"))

ggsave(file.path(output_dir, "plot_DM_rate_inside_outside_G4.pdf"),
       width = 7, height = 5)

# ── Plot 5: Methylation level inside vs outside G4 (per sample) ───────────────
# Note: meth_strict must be loaded; uses percMethylation() which needs the
# full united methylBase object.

pct_matrix <- methylKit::percMethylation(meth_strict)
colnames(pct_matrix) <- meth_strict@sample.ids

# Label each CpG as inside / outside BG4
data_strict <- methylKit::getData(meth_strict)
cpg_coords  <- data.frame(
  chrom = as.character(data_strict$chr),
  start = as.integer(data_strict$start) - 1L,
  end   = as.integer(data_strict$end),
  idx   = seq_len(nrow(data_strict)),
  stringsAsFactors = FALSE
) %>% dplyr::filter(grepl("^chr([0-9]+|X|Y)$", chrom))

in_bg4 <- valr::bed_intersect(cpg_coords, bg4_valr, min_overlap = 0L) %>%
  dplyr::distinct(chrom, start.x, end.x)

cpg_coords$location <- ifelse(
  paste(cpg_coords$chrom, cpg_coords$start) %in%
    paste(in_bg4$chrom, in_bg4$start.x),
  "Inside G4", "Outside G4"
)

cat("Inside G4: ", sum(cpg_coords$location == "Inside G4"),  "\n")
cat("Outside G4:", sum(cpg_coords$location == "Outside G4"), "\n")

treatment  <- meth_strict@treatment
sample_ids <- meth_strict@sample.ids

sample_means <- do.call(rbind, lapply(seq_along(sample_ids), function(i) {
  meth_vals <- pct_matrix[cpg_coords$idx, i]
  data.frame(
    sample    = sample_ids[i],
    condition = ifelse(treatment[i] == 1, "Tumor", "Normal"),
    location  = cpg_coords$location,
    meth_pct  = meth_vals,
    stringsAsFactors = FALSE
  )
})) %>%
  dplyr::filter(!is.na(meth_pct)) %>%
  dplyr::group_by(sample, condition, location) %>%
  dplyr::summarise(mean_meth = mean(meth_pct, na.rm = TRUE), .groups = "drop")

# Violin plot
ggplot(sample_means, aes(x = location, y = mean_meth, fill = condition)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA,
               position = position_dodge(0.9)) +
  geom_point(aes(color = condition),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
             size = 2.5, alpha = 0.9) +
  scale_fill_manual(values  = c("Normal" = "#85B7EB", "Tumor" = "#D85A30")) +
  scale_color_manual(values = c("Normal" = "#0C447C", "Tumor" = "#712B13")) +
  labs(x = NULL, y = "Mean methylation per sample (%)",
       fill = "Group", color = "Group") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave(file.path(output_dir, "plot_meth_inside_outside_G4_violin.pdf"),
       width = 7, height = 5)


# ── SECTION 10: GENOME ANNOTATIONS ───────────────────────────────────────────

# ── Promoters & TSS (hg19 / TxDb) ────────────────────────────────────────────
txdb19    <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoters <- promoters(transcripts(txdb19), upstream = 2000, downstream = 2000)
promoters <- GenomicRanges::reduce(promoters)
promoters_df <- as.data.frame(promoters)
colnames(promoters_df) <- c("chrom", "start", "end", "width", "strand")

TSS_19 <- promoters(transcripts(txdb19), upstream = 0, downstream = 1)
TSS_19 <- GenomeInfoDb::keepStandardChromosomes(TSS_19, pruning.mode = "coarse")
TSS_19 <- GenomicRanges::reduce(TSS_19)
cat("TSS sites:", length(TSS_19), "\n")

# ── SCREEN enhancers (hg19, no promoter overlap) ──────────────────────────────
# Requires the SCREEN bed file from ENCODE / your local copy
enhancer_REF_hg19 <- as.data.frame(
  read.table("/Users/hossein.allahdadi/Downloads/ChIP/rnaSeq/enhancer bed file/enhancer_peaks_hg19.bed")
)
dim(enhancer_REF_hg19)
colnames(enhancer_REF_hg19) <- c("chrom", "start", "end", "ID1", "ID2", "dist")

# Remove promoter-overlapping enhancers
SCREEN_pure_tss <- valr::bed_subtract(enhancer_REF_hg19, promoters_df)
cat("SCREEN enhancers (no promoter overlap):", nrow(SCREEN_pure_tss), "\n")

# ── CpG Islands (UCSC, downloaded) ───────────────────────────────────────────
# Only needed if you want to restrict permutation shuffling to CpG islands
# instead of all RRBS-covered CpGs (see Section 8).
# download.file(
#   url      = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz",
#   destfile = "/Users/hossein.allahdadi/Downloads/ChIP/files/cpgIslandExt.txt.gz"
# )
cpg_islands_raw <- read.table(
  "/Users/hossein.allahdadi/Downloads/ChIP/files/cpgIslandExt.txt.gz",
  header = FALSE, sep = "\t"
)
cpg_islands_df <- data.frame(
  chrom = as.character(cpg_islands_raw$V2),
  start = as.integer(cpg_islands_raw$V3),
  end   = as.integer(cpg_islands_raw$V4),
  stringsAsFactors = FALSE
)
cpg_islands_df <- cpg_islands_df[grep("^chr([0-9]+|X|Y)$", cpg_islands_df$chrom), ]
cat("CpG islands loaded:", nrow(cpg_islands_df), "\n")   # ~27,000


# ── SECTION 11: FILE FORMAT CONVERSIONS ──────────────────────────────────────

# ── BigWig → BED ─────────────────────────────────────────────────────────────
bigwig_path <- "/Users/hossein.allahdadi/Downloads/ChIP/GSE266366(active seq methylation_raw_agustin)/GSM9170395_Active-Seq_SW48_Rep-1_coverage_RPGC.bw"
bigwig      <- rtracklayer::import(bigwig_path, format = "BigWig")
bw_df       <- as.data.frame(bigwig)[, c("seqnames", "start", "end")]
colnames(bw_df) <- c("chrom", "start", "end")

write.table(bw_df,
            file      = sub("\\.bw$", ".bed", bigwig_path),
            sep       = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# ── LiftOver hg38 → hg19 ─────────────────────────────────────────────────────
chain_file <- "/Users/hossein.allahdadi/Downloads/ChIP/GSE266366(active seq methylation_raw_agustin)/hg38ToHg19.over.chain"

gr_hg38 <- GRanges(
  seqnames = bw_df$chrom,
  ranges   = IRanges(start = bw_df$start, end = bw_df$end)
)

ch      <- rtracklayer::import.chain(chain_file)
gr_hg19 <- unlist(liftOver(gr_hg38, ch))
genome(gr_hg19) <- "hg19"

agustin_active1 <- data.frame(
  chrom = as.character(seqnames(gr_hg19)),
  start = start(gr_hg19),
  end   = end(gr_hg19)
)

rtracklayer::export(gr_hg19,
                    con    = sub("\\.bw$", "_hg19.bed", bigwig_path),
                    format = "BED")


# ── SECTION 12: RNA-SEQ ANALYSIS ─────────────────────────────────────────────
# Two datasets merged: GSE95132 (CRC vs adjacent normal) + GSE213402 (CRC vs LM)
# Workflow A: DESeq2 for GSE95132 alone
# Workflow B: limma on log2-FPKM for the merged dataset

library(DESeq2)
library(limma)

# ── 12A: GSE95132 — DESeq2 ───────────────────────────────────────────────────
counts_raw_GSE95132 <- readr::read_table(
  "/Users/hossein.allahdadi/Downloads/ChIP/GSE95656/GSE95132_RNAseq_read_counts.txt"
)
counts_raw_GSE95132 <- as.data.frame(counts_raw_GSE95132)
rownames(counts_raw_GSE95132) <- counts_raw_GSE95132$`"Gene_symbol"`
counts_raw_GSE95132 <- counts_raw_GSE95132[, -1]

# Aggregate duplicate gene symbols by summing counts
aggr_data <- aggregate(. ~ `"Gene_symbol"`,
                       data = counts_raw_GSE95132, FUN = sum)

# Fix HGNC-renamed gene symbols (MARCH → MARCHF, SEPT → SEPTIN)
genes_fixed <- c("DELEC1","MARCHF1","SEPTIN1","MARCHF10","SEPTIN10","MARCHF11",
                 "SEPTIN11","SEPTIN12","SEPTIN14","SEPTIN15","MARCHF2","SEPTIN2",
                 "MARCHF3","SEPTIN3","MARCHF4","SEPTIN4","MARCHF5","SEPTIN5",
                 "MARCHF6","SEPTIN6","MARCHF7","SEPTIN7","MARCHF8","SEPTIN8",
                 "MARCHF9","SEPTIN9")
aggr_data[seq_along(genes_fixed), 1] <- genes_fixed
aggr_data$`"Gene_symbol"` <- gsub('^"|"$', "", aggr_data$`"Gene_symbol"`)
rownames(aggr_data) <- aggr_data$`"Gene_symbol"`
aggr_data <- aggr_data[, -1]
aggr_data <- aggr_data[, -c(21:31)]   # remove extra columns

# Sample metadata
coldata_gse95132 <- data.frame(
  row.names = colnames(aggr_data),
  GSM_names = c("GSM2496685","GSM2496686","GSM2496687","GSM2496688","GSM2496689",
                "GSM2496690","GSM2496691","GSM2496692","GSM2496693","GSM2496694",
                "GSM2496695","GSM2496696","GSM2496697","GSM2496698","GSM2496699",
                "GSM2496700","GSM2496701","GSM2496702","GSM2496703","GSM2496704"),
  patient   = c("32","32","41","41","42","42","15","15","16","16","8","8",
                "50","50","51","51","54","54","60","60"),
  condition = c("Normal","Tumor","Normal","Tumor","Normal","Tumor","Normal","Tumor",
                "Normal","Tumor","Normal","Tumor","Tumor","Normal","Tumor","Normal",
                "Tumor","Normal","Tumor","Normal"),
  stringsAsFactors = FALSE
)
coldata_gse95132$patient   <- factor(coldata_gse95132$patient)
coldata_gse95132$condition <- factor(coldata_gse95132$condition, levels = c("Normal","Tumor"))

# DESeq2
dds_gse95132 <- DESeqDataSetFromMatrix(
  countData = aggr_data,
  colData   = coldata_gse95132,
  design    = ~ condition
)
dds_gse95132 <- DESeq(dds_gse95132)
saveRDS(dds_gse95132, file.path(output_dir, "dds_results_gse95132.rds"))

ddsResult_NormalvsTumor <- results(dds_gse95132,
                                   contrast = c("condition", "Normal", "Tumor"))
ddsresult_gse95132 <- na.omit(as.data.frame(ddsResult_NormalvsTumor))

up_gse95132   <- ddsresult_gse95132[ddsresult_gse95132$log2FoldChange >  1 &
                                      ddsresult_gse95132$padj < 0.05, ]
down_gse95132 <- ddsresult_gse95132[ddsresult_gse95132$log2FoldChange < -1 &
                                      ddsresult_gse95132$padj < 0.05, ]
cat("Up-regulated genes:  ", nrow(up_gse95132),   "\n")
cat("Down-regulated genes:", nrow(down_gse95132), "\n")

# QC plots for GSE95132
vsd_gse95132 <- vst(dds_gse95132, blind = TRUE)
DESeq2::plotPCA(vsd_gse95132, intgroup = "condition")

dists <- dist(t(assay(vsd_gse95132)))
pheatmap::pheatmap(as.matrix(dists),
                   clustering_distance_rows = dists,
                   clustering_distance_cols = dists,
                   main = "Sample-to-sample distances — GSE95132")

# ── 12B: Merged GSE95132 + GSE213402 — limma on log2-FPKM ───────────────────
gse95132_fpkm <- read.table(
  "/Users/hossein.allahdadi/Downloads/ChIP/merge rnaseq 95656_213402/GSE95132_norm_counts_FPKM_GRCh38.p13_NCBI.tsv"
)
colnames(gse95132_fpkm)  <- gse95132_fpkm[1, ]
gse95132_fpkm            <- gse95132_fpkm[-1, ]
rownames(gse95132_fpkm)  <- gse95132_fpkm$GeneID
gse95132_fpkm            <- gse95132_fpkm[, -1]

gse213402_fpkm <- readr::read_table(
  "/Users/hossein.allahdadi/Downloads/ChIP/merge rnaseq 95656_213402/GSE213402_norm_counts_FPKM_GRCh38.p13_NCBI.tsv"
)
gse213402_fpkm           <- as.data.frame(gse213402_fpkm)
rownames(gse213402_fpkm) <- gse213402_fpkm$GeneID
gse213402_fpkm           <- gse213402_fpkm[, -1]

cat("Shared genes:", length(intersect(rownames(gse213402_fpkm),
                                      rownames(gse95132_fpkm))), "\n")

# Merge on common genes
merged_fpkm <- merge(gse213402_fpkm, gse95132_fpkm, by = 0)
rownames(merged_fpkm) <- merged_fpkm$Row.names
merged_fpkm <- merged_fpkm[, -1]

# Remove 3 outlier normal samples identified by PCA
merged_fpkm <- merged_fpkm[, !colnames(merged_fpkm) %in%
                              c("GSM2496685", "GSM2496698", "GSM2496701")]

# Log2 transform (FPKM + 1 pseudocount)
merged_log2 <- log2(sapply(merged_fpkm, as.numeric) + 1)
rownames(merged_log2) <- rownames(merged_fpkm)

# Remove LM samples (keep only CRC vs Normal)
lm_samples  <- paste0("GSM6585", 715:724)
merged_log2 <- merged_log2[, !colnames(merged_log2) %in% lm_samples]

# Sample metadata for merged dataset
samples_merged <- colnames(merged_log2)
condition_vec  <- c(rep("CRC", 8),
                    "CRC","Normal","CRC","Normal","CRC","Normal","CRC","Normal",
                    "CRC","Normal","CRC","CRC","CRC","Normal","Normal","CRC","Normal")
batch_vec      <- c(rep("GSE213402", 8), rep("GSE95132", 17))

coldata_merged <- data.frame(
  row.names = samples_merged,
  condition = factor(condition_vec, levels = c("Normal", "CRC")),
  batch     = factor(batch_vec),
  stringsAsFactors = FALSE
)

# Limma — adjust for batch
design_merged  <- model.matrix(~ 0 + condition, data = coldata_merged)
fit_merged     <- lmFit(merged_log2, design_merged)
fit_merged     <- contrasts.fit(fit_merged,
                                makeContrasts(CRCvsNormal = conditionCRC - conditionNormal,
                                              levels = design_merged))
fit_merged     <- eBayes(fit_merged)
res_merged     <- topTable(fit_merged, coef = "CRCvsNormal", n = Inf, adjust = "BH")

cat("Limma DEGs (FDR < 0.05):", sum(res_merged$adj.P.Val < 0.05, na.rm = TRUE), "\n")


# ══════════════════════════════════════════════════════════════════════════════
# ARCHIVE — Older / exploratory code kept for reference
# Everything below was replaced by the cleaner versions above.
# ══════════════════════════════════════════════════════════════════════════════

# ── A1: Earlier Fisher test (all BG4 CpGs as denominator, not just covered) ──
# Replaced by run_fisher_covered() in Section 7 which uses only RRBS-covered
# BG4 CpGs, giving a more honest denominator.
#
# run_fisher_old <- function(dmc_valr, label) {
#   A <- bed_intersect(bg4_valr, dmc_valr, min_overlap = 1L) %>%
#     distinct(chrom, start.x, end.x) %>% nrow()
#   B <- nrow(bg4_valr) - A
#   C <- nrow(dmc_valr) - (bed_intersect(dmc_valr, bg4_valr, min_overlap = 1L) %>%
#                            distinct(chrom, start.x, end.x) %>% nrow())
#   D <- (nrow(total_valr) - nrow(bg4_valr)) - C
#   ft <- fisher.test(matrix(c(A,C,B,D), nrow=2), alternative="two.sided")
#   cat(label, "— OR:", round(ft$estimate, 3), "| p:", ft$p.value, "\n")
# }
# run_fisher_old(hyper_valr, "Hyper")
# run_fisher_old(hypo_valr,  "Hypo")

# ── A2: BG4-side Fisher (alternative matrix orientation) ─────────────────────
# Same result as run_fisher_covered, different bookkeeping.
# Kept here for cross-checking OR direction.
#
# run_fisher_bg4side <- function(dmc_valr, dmc_label) {
#   A <- bed_intersect(bg4_valr, dmc_valr) %>%
#     dplyr::distinct(chrom, start.x, end.x) %>% nrow()
#   B <- nrow(bg4_valr) - A
#   dmc_in_bg4 <- bed_intersect(dmc_valr, bg4_valr) %>%
#     dplyr::distinct(chrom, start.x, end.x) %>% nrow()
#   C <- nrow(dmc_valr) - dmc_in_bg4
#   D <- (nrow(total_valr) - nrow(bg4_valr)) - C
#   ft <- fisher.test(matrix(c(A,C,B,D), nrow=2), alternative = "two.sided")
#   cat(dmc_label, "OR:", round(ft$estimate, 4), "| p:", ft$p.value, "\n")
# }
# results_fisher_bg4side <- rbind(
#   run_fisher_bg4side(hyper_valr, "Hyper"),
#   run_fisher_bg4side(hypo_valr,  "Hypo")
# )

# ── A3: Permutation test with mclapply (replaced by parLapply version) ────────
# mclapply can have fork-safety issues on macOS; parLapply is more stable.
#
# permutation_test_mclapply <- function(query_valr, subject_valr, ...) {
#   perm_hits <- unlist(parallel::mclapply(1:n_perm, function(i) {
#     shuffled <- valr::bed_shuffle(query_valr, genome_df,
#                                   incl = cpg_universe, seed = i)
#     bed_intersect(shuffled, subject_valr) %>%
#       dplyr::distinct(chrom, start.x, end.x) %>% nrow()
#   }, mc.cores = n_cores))
# }

# ── A4: GRanges-based permutation (shuffle_peaks, no CpG constraint) ─────────
# Simple but biologically incorrect — shuffles peaks anywhere in the genome,
# not restricted to CpG-accessible regions.
#
# shuffle_peaks <- function(peaks, chr_sizes) {
#   chr_pick <- sample(names(chr_sizes), length(peaks), replace = TRUE,
#                      prob = chr_sizes / sum(chr_sizes))
#   peak_w   <- width(peaks)
#   starts   <- mapply(function(chr, w) {
#     max_start <- chr_sizes[chr] - w
#     if (max_start < 1) return(1L)
#     sample.int(max_start, 1L)
#   }, chr_pick, peak_w)
#   GRanges(seqnames = chr_pick,
#           ranges   = IRanges(start = as.integer(starts), width = peak_w))
# }
# observed_hyper <- sum(countOverlaps(bg4_gr, hyper_gr))
# perm_hyper     <- replicate(n_perm,
#   sum(countOverlaps(shuffle_peaks(bg4_gr, chr_sizes), hyper_gr)))
# p_hyper <- mean(perm_hyper >= observed_hyper)

# ── A5: Manual methylDiff object reconstruction ───────────────────────────────
# new("methylDiff") approach needed when loading from BGZ outside methylKit.
# Replaced by the simpler read_bgz() data-frame approach.
#
# read_bgz <- function(bgz_path) {
#   tbx       <- Rsamtools::TabixFile(bgz_path)
#   open(tbx)
#   raw_lines <- scanTabix(tbx)[[1]]
#   close(tbx)
#   df <- data.table::fread(text = paste(raw_lines, collapse = "\n"))
#   colnames(df) <- c("chr","start","end","strand","pvalue","qvalue","meth.diff")
#   as.data.frame(df)
# }
# myDiff_all   <- read_bgz(file.path(bgz_dir, "methylDiff_985977353f49_all.txt.bgz"))
# myDiff_hyper <- read_bgz(file.path(bgz_dir, "methylDiff_985977353f49_hyper.txt.bgz"))
# myDiff_hypo  <- read_bgz(file.path(bgz_dir, "methylDiff_985977353f49_hypo.txt.bgz"))

# ── A6: Raw per-sample bedGraph reads (exploratory; superseded by methylKit) ──
# These were read individually before the methylKit pipeline was set up.
# Kept for reference — individual GSM file paths are recorded here.
#
# GSM2520721_6065_8T_meth <- read.table(
#   ".../GSM2520721_6065_8T_R1_val_1.fq_bismark_pe.bedGraph",
#   sep="\t", header=FALSE)
# colnames(GSM2520721_6065_8T_meth) <- c("chrom","start","end","pct","mC","uC")
# ... (similar for GSM2520722 – GSM2520726)

# ── A7: GSE213402 methylation file (CRC vs LM — separate study) ──────────────
# Loaded separately; not used in the main BG4 analysis pipeline.
#
# meth_gse13402 <- readr::read_table(
#   "/Users/hossein.allahdadi/Downloads/ChIP/GSE213402/GSE213402_CRCvsLM_RRBS_processed.txt"
# )
# colnames(meth_gse13402) <- c("chrom","start","end","mean.x","mean.y","mean.diff",
#   "var.x","var.y","stderr","df","statistic","pvalue","Enddist","Overlap",
#   "CpGIs_range","CpGIs_dist","CpGI_relation","Strand","GeneID")
# meth_gse13402_bed <- meth_gse13402 %>%
#   transmute(
#     chrom  = paste0("chr", gsub('^"|"$', "", as.character(chrom))),
#     start  = as.integer(start),
#     end    = as.integer(end),
#     name   = as.character(GeneID),
#     score  = as.integer(pmin(1000, -10*log10(pmax(pvalue, 1e-300)))),
#     strand = as.character(Strand)
#   )
# # After liftover to hg19:
# meth_gse13402_bed_hg19 <- read.table(
#   "/Users/hossein.allahdadi/Downloads/ChIP/GSE213402/meth_gse13402.hg19.bed")

# ── A8: CpG-in-peak extraction (BG4 / Quadron) ───────────────────────────────
# Finds CG dinucleotides within ChIP-seq peak intervals using BSgenome.
# Saved to BED files; re-read in Sections 6-7 as cpg_df_e3PE_pure etc.
#
# bg4_gr <- makeGRangesFromDataFrame(BG4e5PE_hg19_pure, starts.in.df.are.0based = TRUE)
# bg4_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, bg4_gr)
# hits <- lapply(seq_along(bg4_seq), function(i) {
#   m <- matchPattern("CG", bg4_seq[[i]])
#   if (length(m) == 0) return(NULL)
#   data.frame(
#     chrom = as.character(seqnames(bg4_gr)[i]),
#     start = start(bg4_gr)[i] + start(m) - 2,   # convert to 0-based
#     end   = start(bg4_gr)[i] + start(m) - 1
#   )
# })
# cpg_BG4_peaks <- do.call(rbind, hits)
# write.table(cpg_BG4_peaks, "BG4_cpg_hg19.bed", sep="\t", quote=FALSE,
#             row.names=FALSE, col.names=FALSE)
