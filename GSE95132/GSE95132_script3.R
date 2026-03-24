# ══════════════════════════════════════════════════════════════════════════════
# G4 × METHYLATION × RNA-seq INTEGRATION PIPELINE
# Focus:   Visualise G4–DMC overlaps, test enrichment, integrate with
#          RNA-seq differential expression, pathway analysis
# Datasets: GSE95656 (RRBS), GSE95132 (RNA-seq), BG4/Quadron peaks (hg19)
# Author:  Hossein Allahdadi
#
# ── RELATIONSHIP TO PREVIOUS SCRIPTS ─────────────────────────────────────────
# This script is DOWNSTREAM of both previous pipelines:
#
#  methylation_pipeline_clean.R  →  supplies:
#    hyper_bed_new / hypo_bed_new (DMC BED files, 0-based, end+1 corrected)
#    meth_hyper_10 / meth_hypo_10 (data frames used as hyper/hypo here)
#    dds_gse95132 (DESeq2 object, Section 12A)
#    BG4_hypo_genes_inPromoter / BG4_hyper_genes_inPromoter (Section end)
#
#  chipseq_array_pipeline_clean.R  →  supplies:
#    BG4_e5_pe_peaks_hg19_pure (mock-subtracted peaks, named BG4e5PE_hg19_pure below)
#    quadron_BG4_e5_hg19 (Quadron predictions, named quadron_e5PE_peak_hg19 below)
#    anno_bg4_df (ChIPseeker annotation with peak_id column, Section 9)
#
# Run order:  methylation_pipeline_clean.R
#          →  chipseq_array_pipeline_clean.R
#          →  this script
#
# ── SHARED FUNCTIONS (defined ONCE here) ──────────────────────────────────────
#   to_gr()          also defined in methylation_pipeline_clean.R Section 6
#                    and repeated in the session restore blocks below.
#                    If running all three scripts in one session, skip
#                    re-definition.
#   shuffle_peaks()  the genome-wide (no CpG constraint) version. Archived in
#                    methylation_pipeline_clean.R A4. The better CpG-constrained
#                    version with parLapply lives in Section 8 of that script.
#                    Used here for a quick sanity check only.
#   count_hits()     unique to this script.
#   add_category()   unique to this script.
#
# Pipeline overview:
#   0.  Libraries & shared functions
#   1.  Paths & data setup (alias objects from upstream scripts)
#   2.  Per-peak DM CpG counting & categorisation
#   3.  Euler diagram  (set overlap)
#   4.  UpSet plot     (intersection breakdown)
#   5.  DM CpGs per peak — violin + stacked bar
#   6.  ChIPseeker annotation × methylation category
#   7.  CpG island overlap by methylation category
#   8.  Fisher's exact test   (hardcoded summary counts)
#   9.  Permutation test      (genome-wide shuffle, GRanges)
#  10.  Fisher + permutation visualisation (forest plot, obs vs expected)
#  11.  RNA-seq: DESeq2 paired analysis & robust DEG set
#  12.  RNA-seq: QC plots (PCA, outlier sensitivity)
#  13.  RNA-seq: Volcano, MA, heatmap
#  14.  RNA-seq: GO / KEGG enrichment
#  15.  G4 × RNA-seq integration
#  16.  BG4 promoter gene heatmap
#  ──  ARCHIVE: duplicate/exploratory blocks
# ══════════════════════════════════════════════════════════════════════════════


# ── SECTION 0: LIBRARIES & SHARED FUNCTIONS ───────────────────────────────────

library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(enrichplot)

library(DESeq2)
library(apeglm)
library(pheatmap)

library(UpSetR)
library(eulerr)

library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(tidyr)
library(scales)

# ── Shared function: convert data frame → GRanges ────────────────────────────
# zero_based = TRUE  → input uses 0-based BED start (NarrowPeak, BED files)
# zero_based = FALSE → input is already 1-based (methylKit output)
# NOTE: identical function defined in methylation_pipeline_clean.R Section 6.
#       Skip re-definition if running both scripts in the same R session.
to_gr <- function(df, zero_based = TRUE) {
  GRanges(
    seqnames = df[[1]],
    ranges   = IRanges(start = df[[2]] + as.integer(zero_based),
                       end   = df[[3]])
  )
}

# ── Count hyper / hypo DMC CpGs overlapping each peak ────────────────────────
count_hits <- function(peaks_gr, label, hyper_gr, hypo_gr) {
  df <- data.frame(
    peak_id = seq_along(peaks_gr),
    set     = label,
    hyper   = countOverlaps(peaks_gr, hyper_gr),
    hypo    = countOverlaps(peaks_gr, hypo_gr)
  )
  subset(df, hyper > 0 | hypo > 0)
}

# ── Label each peak by its methylation pattern ───────────────────────────────
add_category <- function(df) {
  df$category <- ifelse(
    df$hyper > 0 & df$hypo == 0, "Hyper only",
    ifelse(df$hypo > 0 & df$hyper == 0, "Hypo only", "Mixed")
  )
  df
}

# ── Chromosome sizes (hg19) ───────────────────────────────────────────────────
chr_sizes <- c(
  chr1  = 249250621, chr2  = 243199373, chr3  = 198022430, chr4  = 191154276,
  chr5  = 180915260, chr6  = 171115067, chr7  = 159138663, chr8  = 146364022,
  chr9  = 141213431, chr10 = 135534747, chr11 = 135006516, chr12 = 133851895,
  chr13 = 115169878, chr14 = 107349540, chr15 = 102531392, chr16 =  90354753,
  chr17 =  81195210, chr18 =  78077248, chr19 =  59128983, chr20 =  63025520,
  chr21 =  48129895, chr22 =  51304566, chrX  = 155270560, chrY  =  59373566
)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


# ── SECTION 1: PATHS & DATA SETUP ────────────────────────────────────────────
# Alias objects from upstream scripts to the names used in this pipeline.
# BG4e5PE_hg19_pure      = BG4_e5_pe_peaks_hg19_pure  (chipseq_array_pipeline_clean.R)
# quadron_e5PE_peak_hg19 = quadron_BG4_e5_hg19         (chipseq_array_pipeline_clean.R)
# hyper_bed_new / hypo_bed_new = end-corrected DMC BED  (methylation_pipeline_clean.R)
#   These have end = end + 1 applied so all widths = 1

# Alias objects from upstream scripts to the names used in this pipeline
BG4e5PE_hg19_pure      <- BG4_e5_pe_peaks_hg19_pure   # from Script 2
quadron_e5PE_peak_hg19 <- quadron_BG4_e5_hg19          # from Script 2

output_dir  <- "/Users/hossein.allahdadi/Downloads/ChIP/GSE95656"
rds_dds     <- "/Users/hossein.allahdadi/Downloads/ChIP/GSE95656/dds_results_gse95132(for 20 samples).rds"

# Create GRanges for all four sets
bg4     <- to_gr(BG4e5PE_hg19_pure,      zero_based = TRUE)
quadron <- to_gr(quadron_e5PE_peak_hg19, zero_based = TRUE)
hyper   <- to_gr(hyper_bed_new,          zero_based = FALSE)
hypo    <- to_gr(hypo_bed_new,           zero_based = FALSE)

cat("Peak set sizes:\n")
cat("  BG4:    ", length(bg4),     "peaks\n")
cat("  Quadron:", length(quadron), "peaks\n")
cat("  Hyper:  ", length(hyper),   "DMC CpGs\n")
cat("  Hypo:   ", length(hypo),    "DMC CpGs\n")


# ── SECTION 2: PER-PEAK DM CpG COUNTING & CATEGORISATION ─────────────────────
# For each BG4 and Quadron peak, count how many hyper/hypo CpGs it contains,
# then classify each peak as Hyper only / Hypo only / Mixed.

bg4_dm     <- count_hits(bg4,     "BG4",     hyper, hypo) |> add_category()
quadron_dm <- count_hits(quadron, "Quadron", hyper, hypo) |> add_category()

cat("BG4 peaks with DM CpGs:     ", nrow(bg4_dm),     "\n")
cat("Quadron peaks with DM CpGs: ", nrow(quadron_dm), "\n")

# Category breakdown
summary_df <- rbind(
  as.data.frame(table(bg4_dm$category))     |> transform(set = "BG4"),
  as.data.frame(table(quadron_dm$category)) |> transform(set = "Quadron")
)
names(summary_df)[1:2] <- c("category", "n")
print(summary_df)


# ── SECTION 3: EULER DIAGRAM ──────────────────────────────────────────────────
# Shows set sizes and overlaps as proportional circles.
# Counts below are from the intersection analysis in Section 2.
# BG4_GSE = BG4 peaks covered by RRBS data (71,204 CpGs in peaks).

euler_counts <- c(
  "BG4"                  = 3867,
  "BG4&Quadron"          =  815,
  "BG4&HypoDMC"          =  174,
  "BG4&HyperDMC"         =   65,
  "BG4&Quadron&HypoDMC"  =    7,
  "BG4&Quadron&HyperDMC" =    2
)

fit_euler <- euler(euler_counts, shape = "circle")
print(error_plot(fit_euler))   # diagnostic: lower residuals = better area fit

plot(
  fit_euler,
  fills     = list(fill  = c("#378ADD","#1D9E75","#7F77DD","#D85A30"),
                   alpha = 0.25),
  edges     = list(col   = c("#185FA5","#0F6E56","#534AB7","#993C1D"),
                   lwd   = 1.5),
  labels    = list(labels   = c("BG4 peaks",
                                "Quadron peaks",
                                "Hypo DMC\n174 peaks",
                                "Hyper DMC\n76 peaks"),
                   fontsize = 10,
                   col      = c("#0C447C","#085041","#3C3489","#712B13"),
                   fontface = "bold"),
  quantities = list(fontsize = 9, col = "#444441"),
  main      = list(label    = "G4 peaks \u00d7 differential methylation (CRC)",
                   fontsize = 13, fontface = "bold")
)


# ── SECTION 4: UPSET PLOT ─────────────────────────────────────────────────────
# UpSet plots handle higher-order intersections better than Venn/Euler.
# "Mixed" = peaks with BOTH hyper and hypo CpGs (11 peaks).

upset_data <- fromExpression(c(
  "BG4"                   = 3617,   # BG4 only (non-DM, non-Quadron)
  "BG4&Quadron"           =  806,   # Quadron but no DM
  "BG4&HyperDMC"          =   65,   # hyper only (no hypo overlap)
  "BG4&HypoDMC"           =  174,   # hypo only (no hyper overlap)
  "BG4&HyperDMC&HypoDMC"  =   11,   # mixed: peak has both hyper and hypo CpGs
  "BG4&Quadron&HyperDMC"  =    2,
  "BG4&Quadron&HypoDMC"   =    7
))

set_colors <- c(BG4 = "#378ADD", Quadron = "#1D9E75",
                HyperDMC = "#D85A30", HypoDMC = "#7F77DD")

upset(
  upset_data,
  sets            = c("BG4","Quadron","HyperDMC","HypoDMC"),
  keep.order      = TRUE,
  sets.bar.color  = unname(set_colors[c("BG4","Quadron","HyperDMC","HypoDMC")]),
  order.by        = "freq",
  decreasing      = TRUE,
  mainbar.y.label = "Intersection size (peaks)",
  sets.x.label    = "Set size (peaks)",
  text.scale      = c(1.4, 1.2, 1.2, 1.0, 1.3, 1.1),
  mb.ratio        = c(0.60, 0.40),
  point.size      = 4, line.size = 1,
  main.bar.color  = "#444441", matrix.color = "#444441",
  matrix.dot.alpha = 0.25,
  shade.color     = "#F1EFE8", shade.alpha = 0.4,
  show.numbers    = "yes", number.angles = 0
)
grid::grid.text(
  "G4 peaks \u00d7 differential methylation (CRC)",
  x = 0.65, y = 0.97, just = "top",
  gp = grid::gpar(fontsize = 13, fontface = "bold", col = "#2C2C2A")
)
grid::grid.text(
  "Quadron \u2282 BG4  \u00b7  mixed = peaks with both hyper and hypo CpGs",
  x = 0.65, y = 0.93, just = "top",
  gp = grid::gpar(fontsize = 9, col = "#73726c")
)


# ── SECTION 5: DM CpGs PER PEAK — VIOLIN & STACKED BAR ───────────────────────

# ── 5A: Violin — distribution of DM CpG counts per peak (log scale) ──────────
plot_dm_violin <- function(df, title) {
  long <- pivot_longer(df, c(hyper, hypo), names_to = "type", values_to = "count")
  long <- subset(long, count > 0)
  long$type <- factor(long$type, levels = c("hyper","hypo"),
                      labels = c("Hypermethylated","Hypomethylated"))

  ggplot(long, aes(x = type, y = count, fill = type)) +
    geom_violin(alpha = 0.4, trim = TRUE, linewidth = 0.4) +
    geom_boxplot(width = 0.12, outlier.size = 0.8, linewidth = 0.4) +
    scale_fill_manual(values = c("Hypermethylated" = "#D85A30",
                                 "Hypomethylated"  = "#7F77DD")) +
    scale_y_log10() +
    labs(title   = title,
         x       = NULL,
         y       = "CpGs per peak (log\u2081\u2080 scale)",
         caption = paste0("n = ", nrow(df), " peaks")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          plot.title      = element_text(face = "bold", size = 12))
}

plot_dm_violin(bg4_dm,     "BG4 peaks") +
  plot_dm_violin(quadron_dm, "Quadron peaks (subset of BG4)")
ggsave(file.path(output_dir, "dm_cpgs_per_peak.pdf"), width = 8, height = 5)

# ── 5B: Stacked bar — Hyper only / Hypo only / Mixed per peak set ────────────
ggplot(summary_df, aes(x = set, y = n, fill = category)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Hyper only" = "#D85A30",
                               "Mixed"      = "#B4506A",
                               "Hypo only"  = "#7F77DD")) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  labs(x = NULL, y = "Number of peaks", fill = NULL,
       title = "Peak methylation categories") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))


# ── SECTION 6: ChIPseeker ANNOTATION × METHYLATION CATEGORY ──────────────────
# Requires anno_bg4_df produced by chipseq_array_pipeline_clean.R Section 9.
# Peak IDs are sequential row indices matching seq_along(bg4).

anno_bg4_df$peak_id <- seq_len(nrow(anno_bg4_df))
anno_bg4_ann <- dplyr::left_join(anno_bg4_df,
                                  bg4_dm[, c("peak_id","category")],
                                  by = "peak_id")
anno_bg4_ann$category[is.na(anno_bg4_ann$category)] <- "No DM CpG"

# Simplify ChIPseeker annotation labels
anno_bg4_ann$annot_simple <- dplyr::case_when(
  grepl("Promoter",   anno_bg4_ann$annotation) ~ "Promoter",
  grepl("Exon",       anno_bg4_ann$annotation) ~ "Exon",
  grepl("Intron",     anno_bg4_ann$annotation) ~ "Intron",
  grepl("Downstream", anno_bg4_ann$annotation) ~ "Downstream",
  grepl("Distal",     anno_bg4_ann$annotation) ~ "Distal intergenic",
  TRUE ~ "Other"
)

# Bar: genomic region proportion per methylation category (DM peaks only)
ggplot(dplyr::filter(anno_bg4_ann, category != "No DM CpG"),
       aes(x = category, fill = annot_simple)) +
  geom_bar(position = "fill", width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Genomic annotation by methylation category (BG4)",
       x = NULL, y = "Proportion", fill = "Region") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Distance to TSS
anno_bg4_gr <- annotatePeak(bg4, TxDb = txdb, tssRegion = c(-2000,2000),
                             verbose = FALSE)
plotDistToTSS(anno_bg4_gr, title = "BG4 — distance to TSS")


# ── SECTION 7: CpG ISLAND OVERLAP BY METHYLATION CATEGORY ────────────────────
# Requires cpg_islands_hg19.bed — UCSC CpG Islands track (0-based).
# Download: UCSC Table Browser → CpG Islands → hg19 → BED format

cgi <- GRanges(
  read.table("/Users/hossein.allahdadi/Downloads/ChIP/files/cpgIslandExt.txt.gz",
             header = FALSE, sep = "\t")[, 2:4] |>
    setNames(c("chr","start","end")) |>
    transform(start = start + 1L)   # 0-based → 1-based
)

bg4_dm$cgi_overlap <- countOverlaps(bg4[bg4_dm$peak_id], cgi) > 0

bg4_dm |>
  dplyr::group_by(category, cgi_overlap) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(cgi_label = ifelse(cgi_overlap, "CGI", "Non-CGI")) |>
  ggplot(aes(x = category, y = n, fill = cgi_label)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("CGI" = "#1D9E75", "Non-CGI" = "#B4B2A9")) +
  labs(title = "CpG island overlap by methylation category",
       x = NULL, y = "Proportion", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))


# ── SECTION 8: FISHER'S EXACT TEST (SUMMARY COUNTS) ──────────────────────────
# These counts come from the per-CpG intersection in methylation_pipeline_clean.R.
# Figures are listed explicitly here for transparency and reproducibility.
#
# NOTE: This is a DIFFERENT Fisher test from:
#   - run_fisher_covered()    in methylation_pipeline_clean.R  (per-CpG enrichment)
#   - run_overlap_fisher()    in chipseq_array_pipeline_clean.R (hyper vs hypo rate)
# This one tests whether hyper/hypo CpGs are more/less frequent IN peaks
# than expected by random chance given the genome-wide proportion.

total_RRBS_CpGs <- 6780708   # total CpGs tested in RRBS
bg4_CpGs_total  <-   71204   # CpGs inside BG4 peaks (RRBS-covered)
bg4_hyper_CpGs  <-     359   # hyper CpGs inside BG4 peaks (10% threshold)
bg4_hypo_CpGs   <-     532   # hypo  CpGs inside BG4 peaks (10% threshold)
total_hyper_all <-  153453   # all hyper CpGs genome-wide
total_hypo_all  <-  422740   # all hypo  CpGs genome-wide

run_fisher_summary <- function(dm_in_peaks, total_in_peaks,
                               dm_total,    total_CpGs, label) {
  a <- dm_in_peaks
  b <- dm_total       - dm_in_peaks
  c <- total_in_peaks - dm_in_peaks
  d <- (total_CpGs - total_in_peaks) - (dm_total - dm_in_peaks)

  mat <- matrix(c(a, c, b, d), nrow = 2,
                dimnames = list(c("DM", "non-DM"),
                                c("Inside peaks", "Outside peaks")))
  cat("\n───", label, "───\n")
  print(mat)

  res <- fisher.test(mat, alternative = "two.sided")
  cat("OR =", round(res$estimate, 3),
      " 95% CI: [", round(res$conf.int[1], 3), ",",
      round(res$conf.int[2], 3), "]\n")
  cat("p-value =", signif(res$p.value, 4), "\n")
  cat("Direction:", ifelse(res$estimate < 1, "DEPLETED", "ENRICHED"), "\n")
  invisible(res)
}

fisher_hyper <- run_fisher_summary(bg4_hyper_CpGs, bg4_CpGs_total,
                                   total_hyper_all, total_RRBS_CpGs,
                                   "Hypermethylated CpGs in BG4 peaks")

fisher_hypo  <- run_fisher_summary(bg4_hypo_CpGs,  bg4_CpGs_total,
                                   total_hypo_all,  total_RRBS_CpGs,
                                   "Hypomethylated CpGs in BG4 peaks")


# ── SECTION 9: PERMUTATION TEST (GENOME-WIDE GRanges SHUFFLE) ────────────────
# IMPORTANT: This is the simple genome-wide shuffle — peaks can land anywhere.
# The CpG-constrained, parallel version using valr/parLapply is in
# methylation_pipeline_clean.R Section 8 and gives more biologically valid results.
# This version is kept as a fast sanity check.

set.seed(42)
n_perm <- 1000

shuffle_peaks <- function(peaks, chr_sizes) {
  chr_pick <- sample(names(chr_sizes), length(peaks), replace = TRUE,
                     prob = chr_sizes / sum(chr_sizes))
  peak_w   <- width(peaks)
  starts   <- mapply(function(chr, w) {
    max_start <- chr_sizes[chr] - w
    if (max_start < 1L) return(1L)
    sample.int(max_start, 1L)
  }, chr_pick, peak_w)
  GRanges(seqnames = chr_pick,
          ranges   = IRanges(start = as.integer(starts), width = peak_w))
}

observed_hyper <- sum(countOverlaps(bg4, hyper))
observed_hypo  <- sum(countOverlaps(bg4, hypo))
cat("Observed hyper overlaps:", observed_hyper, "\n")
cat("Observed hypo overlaps: ", observed_hypo,  "\n")
cat("Running", n_perm, "permutations...\n")

perm_hyper <- replicate(n_perm,
  sum(countOverlaps(shuffle_peaks(bg4, chr_sizes), hyper)))
perm_hypo  <- replicate(n_perm,
  sum(countOverlaps(shuffle_peaks(bg4, chr_sizes), hypo)))

p_hyper <- mean(perm_hyper >= observed_hyper)
p_hypo  <- mean(perm_hypo  >= observed_hypo)

cat("Permutation p (hyper):",
    ifelse(p_hyper == 0, paste0("< 1/", n_perm), p_hyper), "\n")
cat("Permutation p (hypo): ",
    ifelse(p_hypo  == 0, paste0("< 1/", n_perm), p_hypo),  "\n")


# ── SECTION 10: ENRICHMENT VISUALISATION ─────────────────────────────────────

# ── 10A: Permutation distribution ────────────────────────────────────────────
perm_df <- data.frame(
  count = c(perm_hyper, perm_hypo),
  type  = rep(c("Hypermethylated","Hypomethylated"), each = n_perm)
)
obs_df <- data.frame(
  count = c(observed_hyper, observed_hypo),
  type  = c("Hypermethylated","Hypomethylated"),
  label = c(paste0("Observed = ", observed_hyper, "\np = ",
                   ifelse(p_hyper == 0, paste0("< 1/", n_perm), p_hyper)),
            paste0("Observed = ", observed_hypo,  "\np = ",
                   ifelse(p_hypo  == 0, paste0("< 1/", n_perm), p_hypo)))
)

ggplot(perm_df, aes(x = count, fill = type)) +
  geom_histogram(bins = 40, alpha = 0.6, color = NA) +
  geom_vline(data = obs_df, aes(xintercept = count, color = type),
             linewidth = 1, linetype = "dashed") +
  geom_text(data = obs_df,
            aes(x = count, y = Inf, label = label, color = type),
            hjust = 1.1, vjust = 1.5, size = 3.2, inherit.aes = FALSE) +
  scale_fill_manual(values  = c("Hypermethylated" = "#D85A30",
                                "Hypomethylated"  = "#7F77DD")) +
  scale_color_manual(values = c("Hypermethylated" = "#993C1D",
                                "Hypomethylated"  = "#3C3489")) +
  facet_wrap(~type, scales = "free") +
  labs(title    = "Permutation test — DM CpG overlap with BG4 peaks",
       subtitle = paste0(n_perm, " random shuffles of BG4 peak positions"),
       x = "Overlap count (random peaks)", y = "Frequency") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "#73726c"))

ggsave(file.path(output_dir, "permutation_test_BG4.pdf"), width = 10, height = 5)

# ── 10B: Forest plot (Odds ratio) ────────────────────────────────────────────
or_df <- data.frame(
  type = c("Hypermethylated","Hypomethylated"),
  or   = c(0.217, 0.112),
  lo   = c(0.195, 0.103),
  hi   = c(0.241, 0.122)
)

p_forest <- ggplot(or_df, aes(x = or, y = type, color = type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#888780") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.8) +
  geom_point(size = 4) +
  scale_x_log10(limits = c(0.08, 1.5)) +
  scale_color_manual(values = c("Hypermethylated" = "#D85A30",
                                "Hypomethylated"  = "#7F77DD")) +
  labs(title = "Odds ratio (log scale) with 95% CI",
       x = "Odds ratio", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 11))

# ── 10C: Observed vs expected counts ─────────────────────────────────────────
peak_fraction <- bg4_CpGs_total / total_RRBS_CpGs

obs_exp_df <- data.frame(
  type     = rep(c("Hypermethylated","Hypomethylated"), each = 2),
  category = rep(c("Observed","Expected"), 2),
  count    = c(bg4_hyper_CpGs, round(total_hyper_all * peak_fraction),
               bg4_hypo_CpGs,  round(total_hypo_all  * peak_fraction))
)

p_obsexp <- ggplot(obs_exp_df,
                   aes(x = type, y = count,
                       fill = interaction(type, category))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = NA) +
  scale_fill_manual(
    values = c("Hypermethylated.Observed" = "#D85A30",
               "Hypermethylated.Expected" = "#F5C4B3",
               "Hypomethylated.Observed"  = "#7F77DD",
               "Hypomethylated.Expected"  = "#CECBF6"),
    labels = c("Hyper expected","Hyper observed",
               "Hypo expected", "Hypo observed")
  ) +
  geom_text(aes(label = count), position = position_dodge(0.6),
            vjust = -0.4, size = 3.2) +
  labs(title = "Observed vs expected CpGs in BG4 peaks",
       x = NULL, y = "CpG count", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11))

p_forest + p_obsexp
ggsave(file.path(output_dir, "fisher_test_plot.pdf"), width = 10, height = 4)


# ── SECTION 11: RNA-seq — DESeq2 PAIRED ANALYSIS & ROBUST DEG SET ────────────
# [UPSTREAM] The DESeq2 object was created by methylation_pipeline_clean.R
# Section 12A (dds_gse95132). This section extends that analysis:
#   1. Refit with the paired design (~ patient + condition)
#   2. Extract shrunken LFC results
#   3. Run sensitivity analysis removing 3 outlier patients
#   4. Define a ROBUST DEG set = genes significant in both full and clean analyses

dds <- readRDS(rds_dds)

# Refit with paired design to control for inter-patient variability
design(dds) <- ~ patient + condition
dds         <- DESeq(dds)

# Shrink LFC (apeglm — requires coef name, not contrast)
res <- lfcShrink(dds, coef = "condition_Tumor_vs_Normal", type = "apeglm")
summary(res)

# Build results data frame
res_df <- as.data.frame(res) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::filter(!is.na(padj))

# Gene symbols — rownames are already symbols in this dataset
res_df$gene_name <- res_df$gene_id

res_df$direction <- dplyr::case_when(
  res_df$padj < 0.05 & res_df$log2FoldChange >  1 ~ "Up",
  res_df$padj < 0.05 & res_df$log2FoldChange < -1 ~ "Down",
  TRUE ~ "NS"
)

cat("Full dataset DEGs — Up:", sum(res_df$direction == "Up"),
    "| Down:", sum(res_df$direction == "Down"), "\n")

# ── Sensitivity analysis: remove outlier patients 41, 42, 51 ─────────────────
# These patients cluster with the wrong condition in PCA (see Section 12).
dds_clean <- dds[, !colData(dds)$patient %in% c(41, 42, 51)]
dds_clean$patient   <- droplevels(dds_clean$patient)
dds_clean$condition <- droplevels(dds_clean$condition)
cat("Samples in clean dataset:", ncol(dds_clean), "\n")
print(table(dds_clean$condition))

dds_clean <- DESeq(dds_clean)
res_clean <- results(dds_clean, name = "condition_Tumor_vs_Normal", alpha = 0.05)

res_clean_df <- as.data.frame(res_clean) |>
  dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

cat("Clean dataset DEGs:", nrow(res_clean_df), "\n")

# ── Robust DEG set = significant in BOTH analyses ────────────────────────────
common <- intersect(rownames(res_clean_df),
                    res_df$gene_id[res_df$direction != "NS"])

res_df$robust  <- res_df$gene_id %in% common
robust_degs    <- dplyr::filter(res_df, robust, direction != "NS")

cat("Shared/robust DEGs: ", length(common),
    "| Up:", sum(robust_degs$direction == "Up"),
    "| Down:", sum(robust_degs$direction == "Down"), "\n")

# Pre-filtered subsets used in plots
df_ns     <- res_df[res_df$direction == "NS", ]
df_nonrob <- res_df[res_df$direction != "NS" & !res_df$robust, ]
df_robust <- res_df[res_df$direction != "NS" &  res_df$robust, ]

top_labels <- robust_degs |> dplyr::arrange(padj) |> dplyr::slice_head(n = 25)

# VST-normalised matrix for plotting
vsd <- vst(dds, blind = TRUE)

ann_col <- data.frame(
  condition = colData(dds)$condition,
  patient   = as.factor(colData(dds)$patient),
  row.names = colnames(dds)
)
ann_colors <- list(condition = c(Tumor = "#D85A30", Normal = "#378ADD"))


# ── SECTION 12: RNA-seq QC PLOTS ─────────────────────────────────────────────

# PCA
pca_data <- DESeq2::plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition,
                               label = as.character(colData(dds)$patient))) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c(Tumor = "#D85A30", Normal = "#378ADD")) +
  labs(title = "PCA — VST normalised counts",
       x = paste0("PC1: ", pct_var[1], "% variance"),
       y = paste0("PC2: ", pct_var[2], "% variance"),
       color = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
print(p_pca)


# ── SECTION 13: VOLCANO, MA, HEATMAP ─────────────────────────────────────────

# ── Volcano ───────────────────────────────────────────────────────────────────
p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = df_ns,     color = "#D3D1C7", size = 0.8, alpha = 0.5) +
  geom_point(data = df_nonrob, aes(color = direction), size = 1,   alpha = 0.35) +
  geom_point(data = df_robust, aes(color = direction), size = 1.8, alpha = 0.85) +
  geom_text_repel(data = top_labels,
                  aes(label = gene_id, color = direction),
                  size = 2.8, max.overlaps = 20, segment.size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "#888780", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "#888780", linewidth = 0.4) +
  scale_color_manual(values = c(Up = "#D85A30", Down = "#378ADD")) +
  annotate("text", x =  4.5, y = max(-log10(res_df$padj), na.rm = TRUE) * 0.9,
           size = 3.2, color = "#D85A30",
           label = paste0("Up: ", sum(robust_degs$direction == "Up"))) +
  annotate("text", x = -4.5, y = max(-log10(res_df$padj), na.rm = TRUE) * 0.9,
           size = 3.2, color = "#378ADD",
           label = paste0("Down: ", sum(robust_degs$direction == "Down"))) +
  labs(title    = "Volcano — Tumor vs Normal (paired design)",
       subtitle = "Bright = robust DEGs (significant in full + clean dataset)",
       x = expression(log[2]~fold~change),
       y = expression(-log[10]~adjusted~p-value),
       color = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "#73726c"),
        legend.position = "none")
print(p_volcano)

# ── MA plot ───────────────────────────────────────────────────────────────────
p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(data = dplyr::filter(res_df, direction == "NS"),
             color = "#D3D1C7", size = 0.8, alpha = 0.4) +
  geom_point(data = dplyr::filter(res_df, robust, direction != "NS"),
             aes(color = direction), size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = c(-1, 0, 1),
             linetype   = c("dashed","solid","dashed"),
             color      = "#888780", linewidth = 0.4) +
  scale_color_manual(values = c(Up = "#D85A30", Down = "#378ADD")) +
  labs(title = "MA plot — Tumor vs Normal",
       x = expression(log[10]~mean~expression),
       y = expression(log[2]~fold~change), color = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")
print(p_ma)

# ── Heatmap — top 50 robust DEGs ─────────────────────────────────────────────
top50_ids <- robust_degs |> dplyr::arrange(padj) |>
  dplyr::slice_head(n = 50) |> dplyr::pull(gene_id)

mat50       <- assay(vsd)[top50_ids, ]
mat50_scaled <- t(scale(t(mat50)))

pheatmap(mat50_scaled,
         annotation_col    = ann_col,
         annotation_colors = ann_colors,
         cluster_rows      = TRUE, cluster_cols = TRUE,
         show_rownames     = TRUE, show_colnames = FALSE,
         treeheight_row    = FALSE, treeheight_col = FALSE,
         color        = colorRampPalette(c("#378ADD","white","#D85A30"))(100),
         fontsize_row = 7, border_color = NA,
         main         = "Top 50 robust DEGs (z-score)")


# ── SECTION 14: GO / KEGG ENRICHMENT ─────────────────────────────────────────
# Uses robust DEG set only to reduce false positives.

get_entrez <- function(gene_symbols) {
  mapIds(org.Hs.eg.db,
         keys      = gene_symbols,
         column    = "ENTREZID",
         keytype   = "SYMBOL",
         multiVals = "first") |>
    na.omit() |> unique()
}

universe   <- get_entrez(res_df$gene_id)
up_entrez  <- get_entrez(robust_degs$gene_id[robust_degs$direction == "Up"])
dn_entrez  <- get_entrez(robust_degs$gene_id[robust_degs$direction == "Down"])

go_up  <- enrichGO(gene = up_entrez, universe = universe, OrgDb = org.Hs.eg.db,
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                   readable = TRUE)
go_dn  <- enrichGO(gene = dn_entrez, universe = universe, OrgDb = org.Hs.eg.db,
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                   readable = TRUE)
kegg_up <- enrichKEGG(gene = up_entrez, organism = "hsa", pvalueCutoff = 0.05)
kegg_dn <- enrichKEGG(gene = dn_entrez, organism = "hsa", pvalueCutoff = 0.05)

p_go_up   <- dotplot(go_up,   showCategory = 15) + ggtitle("GO:BP — upregulated (robust)")
p_go_dn   <- dotplot(go_dn,   showCategory = 15) + ggtitle("GO:BP — downregulated (robust)")
p_kegg_up <- dotplot(kegg_up, showCategory = 15) + ggtitle("KEGG — upregulated (robust)")
p_kegg_dn <- dotplot(kegg_dn, showCategory = 15) + ggtitle("KEGG — downregulated (robust)")

print(p_go_up); print(p_go_dn); print(p_kegg_up); print(p_kegg_dn)


# ── SECTION 15: G4 × RNA-seq INTEGRATION ─────────────────────────────────────
# Join G4 peak annotation (from ChIPseeker) with robust DEGs.
# anno_bg4_ann was built in Section 6 and includes peak_id + category.

g4_genes <- anno_bg4_ann |>
  dplyr::left_join(bg4_dm[, c("peak_id","category","hyper","hypo")],
                   by = "peak_id") |>
  dplyr::filter(!is.na(category)) |>
  dplyr::select(gene_name = geneId, category, hyper, hypo)

integrated <- robust_degs |>
  dplyr::inner_join(g4_genes, by = c("gene_id" = "gene_name"))

cat("Robust DEGs with nearby DM G4 peak:", nrow(integrated), "\n")

# Violin: expression log2FC by G4 category
p_violin_g4 <- ggplot(integrated,
                      aes(x = category, y = log2FoldChange, fill = category)) +
  geom_violin(alpha = 0.4, trim = TRUE, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.size = 0.8, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#888780") +
  scale_fill_manual(values = c("Hyper only" = "#D85A30",
                               "Mixed"      = "#B4506A",
                               "Hypo only"  = "#7F77DD")) +
  labs(title    = "Gene expression by G4 methylation category",
       subtitle = "Robust DEGs only",
       x = NULL, y = "log\u2082 fold change (Tumor vs Normal)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "#73726c"))
print(p_violin_g4)

# Volcano coloured by G4 status
res_df$g4_category <- g4_genes$category[
  match(res_df$gene_id, g4_genes$gene_name)]
res_df$g4_category[is.na(res_df$g4_category)] <- "No G4 peak"

p_volcano_g4 <- ggplot(res_df,
                        aes(x = log2FoldChange, y = -log10(padj),
                            color = g4_category, size = g4_category)) +
  geom_point(data = dplyr::filter(res_df, g4_category == "No G4 peak"),
             color = "#D3D1C7", size = 0.8, alpha = 0.4) +
  geom_point(data = dplyr::filter(res_df, g4_category != "No G4 peak"),
             alpha = 0.85) +
  geom_text_repel(
    data = dplyr::filter(res_df, padj < 0.05, abs(log2FoldChange) > 1,
                         g4_category != "No G4 peak", robust),
    aes(label = gene_name), size = 3, max.overlaps = 20, segment.size = 0.3
  ) +
  scale_color_manual(values = c("Hyper only" = "#D85A30", "Mixed" = "#B4506A",
                                "Hypo only"  = "#7F77DD", "No G4 peak" = "#D3D1C7")) +
  scale_size_manual( values = c("Hyper only" = 2.5, "Mixed" = 2.5,
                                "Hypo only"  = 2.5, "No G4 peak" = 0.8)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "#888780", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "#888780", linewidth = 0.4) +
  labs(title    = "Volcano — genes coloured by G4 methylation status",
       subtitle = "Labelled = robust DEGs with DM G4 peak",
       x = "log\u2082 fold change",
       y = "-log\u2081\u2080 adjusted p-value", color = NULL) +
  guides(size = "none") +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "#73726c"))
print(p_volcano_g4)

# Save all RNA-seq + G4 integration plots
pdf(file.path(output_dir, "RNAseq_G4_integration.pdf"),
    width = 10, height = 7, onefile = TRUE)
print(p_pca); print(p_volcano); print(p_ma)
print(p_go_up); print(p_go_dn); print(p_kegg_up); print(p_kegg_dn)
print(p_violin_g4); print(p_volcano_g4)
dev.off()
message("Saved: RNAseq_G4_integration.pdf")


# ── SECTION 16: BG4 PROMOTER GENE HEATMAP ────────────────────────────────────
# [UPSTREAM] BG4_hypo_genes_inPromoter and BG4_hyper_genes_inPromoter are
# produced by the findOverlaps → mapIds block in methylation_pipeline_clean.R.

all_promo_genes <- unique(na.omit(c(
  BG4_hypo_genes_inPromoter,
  BG4_hyper_genes_inPromoter
)))

genes_in_rna <- all_promo_genes[all_promo_genes %in% rownames(dds)]
missing      <- all_promo_genes[!all_promo_genes %in% rownames(dds)]

cat("Hypo promoter genes:       ", length(na.omit(BG4_hypo_genes_inPromoter)),  "\n")
cat("Hyper promoter genes:      ", length(na.omit(BG4_hyper_genes_inPromoter)), "\n")
cat("Total unique:              ", length(all_promo_genes), "\n")
cat("Found in RNA-seq:          ", length(genes_in_rna),   "\n")
cat("Missing:                   ", length(missing),         "\n")
if (length(missing) > 0) print(missing)

# Build scaled expression matrix
mat        <- assay(vsd)[genes_in_rna, ]
mat_scaled <- t(scale(t(mat)))

# Row annotation: G4 group + DEG status
ann_row <- data.frame(
  G4_group = dplyr::case_when(
    genes_in_rna %in% na.omit(BG4_hypo_genes_inPromoter) &
      genes_in_rna %in% na.omit(BG4_hyper_genes_inPromoter) ~ "Both",
    genes_in_rna %in% na.omit(BG4_hypo_genes_inPromoter)  ~ "Hypo only",
    genes_in_rna %in% na.omit(BG4_hyper_genes_inPromoter) ~ "Hyper only"
  ),
  DEG = dplyr::case_when(
    genes_in_rna %in% robust_degs$gene_id[robust_degs$direction == "Up"]   ~ "Up",
    genes_in_rna %in% robust_degs$gene_id[robust_degs$direction == "Down"] ~ "Down",
    TRUE ~ "NS"
  ),
  row.names = genes_in_rna
)

ann_colors_promo <- list(
  condition = c(Tumor = "#D85A30", Normal = "#378ADD"),
  G4_group  = c("Hypo only" = "#7F77DD", "Hyper only" = "#D85A30",
                "Both"      = "#B4506A"),
  DEG       = c(Up = "#D85A30", Down = "#378ADD", NS = "#D3D1C7")
)

# Manual row order: cluster within each G4 group, then concatenate
cluster_within <- function(names, mat) {
  if (length(names) < 2) return(names)
  names[hclust(dist(mat[names, ]))$order]
}
hypo_rows  <- genes_in_rna[ann_row$G4_group == "Hypo only"]
hyper_rows <- genes_in_rna[ann_row$G4_group == "Hyper only"]
both_rows  <- genes_in_rna[ann_row$G4_group == "Both"]
row_order  <- c(cluster_within(hypo_rows,  mat_scaled),
                cluster_within(both_rows,  mat_scaled),
                cluster_within(hyper_rows, mat_scaled))

# Manual column order: cluster within Normal, then within Tumor
normal_idx <- which(colData(dds)$condition == "Normal")
tumor_idx  <- which(colData(dds)$condition == "Tumor")
cluster_col_idx <- function(idx) {
  if (length(idx) < 2) return(idx)
  idx[hclust(dist(t(mat_scaled[, idx])))$order]
}
col_order <- c(cluster_col_idx(normal_idx), cluster_col_idx(tumor_idx))

pheatmap(
  mat_scaled[row_order, col_order],
  annotation_col    = ann_col,
  annotation_row    = ann_row,
  annotation_colors = ann_colors_promo,
  cluster_rows      = FALSE, cluster_cols  = FALSE,
  show_rownames     = TRUE,  show_colnames = FALSE,
  gaps_row  = c(length(hypo_rows), length(hypo_rows) + length(both_rows)),
  gaps_col  = length(normal_idx),
  color        = colorRampPalette(c("#378ADD","white","#D85A30"))(100),
  fontsize_row = 8, border_color = NA,
  main = paste0("BG4 promoter genes — hypo + hyper (n=", length(genes_in_rna), ")")
)

# DEG-only heatmap (genes with Up or Down status)
DEG_genes_up   <- rownames(ann_row[ann_row$DEG == "Up",   , drop = FALSE])
DEG_genes_down <- rownames(ann_row[ann_row$DEG == "Down", , drop = FALSE])
cat("DEG genes Up:", length(DEG_genes_up), "| Down:", length(DEG_genes_down), "\n")

up_order   <- if (length(DEG_genes_up)   >= 2) DEG_genes_up[hclust(dist(mat_scaled[DEG_genes_up,   ]))$order] else DEG_genes_up
down_order <- if (length(DEG_genes_down) >= 2) DEG_genes_down[hclust(dist(mat_scaled[DEG_genes_down,]))$order] else DEG_genes_down
degs_use   <- c(up_order, down_order)

pheatmap(
  mat_scaled[degs_use, col_order],
  annotation_col = ann_col,
  annotation_row = ann_row[degs_use, , drop = FALSE],
  cluster_rows   = FALSE, cluster_cols = FALSE,
  show_rownames  = TRUE,  show_colnames = FALSE,
  gaps_row       = length(up_order),
  gaps_col       = length(normal_idx),
  color        = colorRampPalette(c("#378ADD","white","#D85A30"))(100),
  fontsize_row = 8, border_color = NA,
  main = paste0("BG4 promoter DEGs (n=", length(degs_use), ")")
)


# ══════════════════════════════════════════════════════════════════════════════
# ARCHIVE — Duplicate / superseded blocks
# ══════════════════════════════════════════════════════════════════════════════

# ── A1: Second UpSet block with slightly different counts ─────────────────────
# The script contained two near-identical upset() calls. The first used
# BG4&HypoDMC = 178 / BG4&HyperDMC = 74 (no Mixed category);
# the second (kept in Section 4) added the Mixed = 11 term.
# Only the second version is biologically complete.
#
# upset_data_v1 <- fromExpression(c(
#   "BG4"                  = 3867,
#   "BG4&Quadron"          =  815,
#   "BG4&HypoDMC"          =  178,
#   "BG4&HyperDMC"         =   74,
#   "BG4&Quadron&HypoDMC"  =    7,
#   "BG4&Quadron&HyperDMC" =    2
# ))

# ── A2: Duplicate library loading blocks ─────────────────────────────────────
# GenomicRanges / ggplot2 / patchwork / dplyr were imported 4–5 times in
# the original script. Consolidated into Section 0.

# ── A3: Duplicate RNA-seq library block ──────────────────────────────────────
# DESeq2 / ggplot2 / pheatmap / ggrepel / clusterProfiler were listed twice
# in adjacent code blocks. First occurrence kept in Section 0.

# ── A4: Session restore block ─────────────────────────────────────────────────
# The original script contained a standalone "SESSION RESTORE" section that
# re-defined to_gr(), count_hits(), add_category(), and re-ran the full DESeq2
# pipeline from scratch. This is now redundant because the pipeline above is
# already self-contained. If you need a minimal restore after a crash, run
# Section 0 (libraries + functions), then Section 1 (alias objects), then
# Section 11 (re-run DESeq2).

# ── A5: Duplicate volcano annotation blocks ───────────────────────────────────
# The original volcano had duplicate annotate() calls placing "Up: N" text at
# both y=1 and y=max*0.9. Only the y=max*0.9 position is kept in Section 13.
#
# annotate("text", x = 4, y = 1, ...)        # ← duplicate, removed
# annotate("text", x = 4.5, y = max*0.9, ...) # ← kept

# ── A6: geom_bar with border = NA ─────────────────────────────────────────────
# border is not a valid ggplot2 geom_bar aesthetic; the correct argument is
# color = NA. Fixed in Section 10C.
#
# geom_bar(..., border = NA)   # ← invalid
# geom_bar(..., color  = NA)   # ← correct (Section 10C)
