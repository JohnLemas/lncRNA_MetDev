#!/usr/bin/env Rscript
# ==============================================================================
# lncRNA_circos.R
#
# Circos plot of mRNA, TE and lncRNA features (circlize).
#
# Track layout (outside -> inside)
#   0. chromosome axis + labels
#   1. mRNA density                     (heatmap)
#   2. TE density                       (heatmap)
#   3. lncRNA genes per window          (stacked bars: counts by Category)
#   4. lncRNA expression                (line: mean of per-gene TPM in each window)
#   centre: stacked bars of lncRNA gene counts by Category for each chromosome
#
# Inputs
#   1. karyotype file     (Chr, Start, End with a header; nuclear chromosomes only)
#   2. lncRNA annotation  (A.tha_final_annotation.tsv; needs Chr, Start, Stop,
#                          lncRNA_gene, Category; mean_TPM_* columns feed the TPM track)
#   3. mRNA annotation    (GFF3/GTF-style, 9 tab-separated columns)
#   4. TE annotation      (GFF3/GTF-style, 9 tab-separated columns)
#
# Outputs (in out_dir)
#   density_mRNA.tsv, density_TE.tsv, lncRNA_category_counts.tsv, tpm_lncRNA.tsv
#   circos_heatmap.svg / .pdf / .png
#
# Usage
#   Rscript lncRNA_circos.R --window 50000 --stack-window 500000 --tpm-window 250000
#   Rscript lncRNA_circos.R --help
#
# install.packages(c("tidyverse", "circlize", "optparse"))
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(circlize)
  library(optparse)
})

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-k", "--karyotype"), type = "character",
              default = "A.tha/01_inputs/A.tha.karyotype.txt",
              help = "karyotype file (Chr, Start, End) [default: %default]"),
  make_option(c("-l", "--lncrna"), type = "character",
              default = "A.tha/03_outputs/11_int_class/A.tha_final_annotation.tsv",
              help = "lncRNA annotation table [default: %default]"),
  make_option(c("-m", "--mrna"), type = "character",
              default = "A.tha/01_inputs/A.tha_mRNA.gtf",
              help = "mRNA annotation, GFF3/GTF [default: %default]"),
  make_option(c("-t", "--te"), type = "character",
              default = "A.tha/01_inputs/A.tha_TE.gff3",
              help = "TE annotation, GFF3/GTF [default: %default]"),
  make_option(c("-o", "--outdir"), type = "character",
              default = "A.tha/03_outputs/13_circos",
              help = "output directory [default: %default]"),
  make_option(c("-w", "--window"), type = "integer", default = 50000L,
              help = "window size (bp) for the mRNA and TE heatmaps [default: %default]"),
  make_option("--stack-window", dest = "stack_window", type = "integer", default = 500000L,
              help = "window size (bp) for the lncRNA stacked-count track [default: %default]"),
  make_option("--heat-quantile", dest = "heat_quantile", type = "double", default = 0.99,
              help = paste("heatmap colour scale tops out at this quantile of the non-zero",
                           "window values, so a few extreme windows don't wash out the rest [default: %default]")),
  make_option("--tpm-window", dest = "tpm_window", type = "integer", default = 250000L,
              help = paste("window size (bp) for the lncRNA TPM line; lncRNAs are sparse, so this",
                           "is usually wider than --window [default: %default]")),
  make_option("--tpm-stat", dest = "tpm_stat", type = "character", default = "mean",
              help = paste("per-transcript expression: 'mean' (average of the mean_TPM_* columns)",
                           "or 'max' (highest mean_TPM_* column). Per gene, the highest isoform is used [default: %default]")),
  make_option("--linear-tpm", dest = "linear_tpm", action = "store_true", default = FALSE,
              help = "plot raw TPM instead of log10(TPM + 1) [default: off]"),
  make_option("--mrna-feature", dest = "mrna_feature", type = "character", default = "^(mRNA|transcript)$",
              help = paste("regex for GFF/GTF column 3 in the mRNA file; AGAT's GTF writes",
                           "mRNAs as 'transcript' [default: %default]")),
  make_option("--te-feature", dest = "te_feature", type = "character", default = NULL,
              help = "regex for GFF/GTF column 3 in the TE file [default: use every row]"),
  make_option("--expressed-only", dest = "expressed_only", action = "store_true", default = FALSE,
              help = "restrict lncRNAs to the expressed set (TPM filter below) [default: off]"),
  make_option("--lower-tpm", dest = "lower_tpm", type = "double", default = 0.1,
              help = "with --expressed-only: keep max mean TPM > this [default: %default]"),
  make_option("--upper-tpm", dest = "upper_tpm", type = "double", default = 50,
              help = "with --expressed-only: keep max mean TPM <= this [default: %default]"),
  make_option("--workdir", type = "character", default = NULL,
              help = "setwd() here first; relative paths above are resolved from it")
)
opt <- parse_args(
  OptionParser(option_list = option_list,
               description = "Circos plot of mRNA, TE and lncRNA tracks (circlize)."),
  args = if (interactive()) character(0) else commandArgs(trailingOnly = TRUE)
)

if (is.null(opt$workdir) && interactive()) {
  opt$workdir <- "~/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/chapter4/lncRNA_met_dev"
}
if (!is.null(opt$workdir)) setwd(opt$workdir)

karyotype_file  <- opt$karyotype
lncRNA_file     <- opt$lncrna
mrna_file       <- opt$mrna
te_file         <- opt$te
out_dir         <- opt$outdir
window_bp       <- opt$window
stack_window_bp <- opt$stack_window
heat_quantile   <- opt$heat_quantile
tpm_window_bp   <- opt$tpm_window
tpm_stat        <- match.arg(opt$tpm_stat, c("mean", "max"))
log_tpm         <- !opt$linear_tpm
mrna_feature    <- opt$mrna_feature
te_feature      <- opt$te_feature
expressed_only  <- opt$expressed_only
LOWER_TPM       <- opt$lower_tpm
UPPER_TPM       <- opt$upper_tpm

for (f in c(karyotype_file, lncRNA_file, mrna_file, te_file)) {
  if (!file.exists(f)) stop("Input file not found: ", f, " (working directory: ", getwd(), ")")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------

# "1", "chr1" and "Chr1" all become "Chr1"
norm_chr <- function(x) paste0("Chr", sub("^[Cc][Hh][Rr]", "", trimws(x)))

# Minimal GFF/GTF reader: seqid, feature, start, end
read_gff <- function(path) {
  x <- read.delim(path, header = FALSE, comment.char = "#", quote = "",
                  stringsAsFactors = FALSE, fill = TRUE)
  if (ncol(x) < 5) stop(path, " does not look like a 9-column GFF/GTF file")
  out <- tibble(seqid   = as.character(x[[1]]),
                feature = as.character(x[[3]]),
                start   = suppressWarnings(as.numeric(x[[4]])),
                end     = suppressWarnings(as.numeric(x[[5]])))
  filter(out, !is.na(start), !is.na(end))
}

prep_features <- function(df, karyotype, feature_regex = NULL, label = "features") {
  if (!is.null(feature_regex)) {
    available <- sort(unique(df$feature))
    df <- filter(df, grepl(feature_regex, feature))
    if (nrow(df) == 0) {
      stop("No ", label, " matched feature regex '", feature_regex,
           "'. Feature types present: ", paste(available, collapse = ", "))
    }
  }
  df <- df %>%
    mutate(Chr = norm_chr(seqid), Start = pmin(start, end), End = pmax(start, end)) %>%
    filter(Chr %in% karyotype$Chr)
  if (nrow(df) == 0) stop("No ", label, " fall on the chromosomes in the karyotype; check sequence names.")
  df
}

# Counts of features per window (feature assigned to the window containing its midpoint)
window_density <- function(feat, karyotype, window_bp) {
  map_dfr(seq_len(nrow(karyotype)), function(i) {
    chr    <- karyotype$Chr[i]
    len    <- karyotype$End[i]
    starts <- seq(1, len, by = window_bp)
    ends   <- pmin(starts + window_bp - 1, len)
    f      <- feat[feat$Chr == chr, ]
    idx    <- findInterval((f$Start + f$End) / 2, starts)
    idx    <- idx[idx >= 1 & idx <= length(starts)]
    tibble(Chr = chr, Start = starts, End = ends,
           Value = tabulate(idx, nbins = length(starts)))
  })
}

# Per-window counts split by Category: one column per category, in the order given
window_category <- function(feat, karyotype, window_bp, cats) {
  map_dfr(seq_len(nrow(karyotype)), function(i) {
    chr    <- karyotype$Chr[i]
    len    <- karyotype$End[i]
    starts <- seq(1, len, by = window_bp)
    ends   <- pmin(starts + window_bp - 1, len)
    f      <- feat[feat$Chr == chr, ]
    idx    <- findInterval((f$Start + f$End) / 2, starts)
    ok     <- idx >= 1 & idx <= length(starts)
    out    <- tibble(Chr = chr, Start = starts, End = ends)
    for (cat in cats) out[[cat]] <- tabulate(idx[ok & f$Category == cat], nbins = length(starts))
    out
  })
}

# Mean of a per-feature value per window (windows with no features get 0)
window_mean <- function(feat, value_col, karyotype, window_bp) {
  map_dfr(seq_len(nrow(karyotype)), function(i) {
    chr    <- karyotype$Chr[i]
    len    <- karyotype$End[i]
    starts <- seq(1, len, by = window_bp)
    ends   <- pmin(starts + window_bp - 1, len)
    f      <- feat[feat$Chr == chr, ]
    idx    <- findInterval((f$Start + f$End) / 2, starts)
    ok     <- idx >= 1 & idx <= length(starts) & !is.na(f[[value_col]])
    idx    <- idx[ok]; v <- f[[value_col]][ok]
    n      <- tabulate(idx, nbins = length(starts))
    sums   <- numeric(length(starts))
    if (length(idx)) {
      agg <- tapply(v, idx, sum)
      sums[as.integer(names(agg))] <- agg
    }
    tibble(Chr = chr, Start = starts, End = ends,
           Value = ifelse(n > 0, sums / n, 0), n_genes = n)
  })
}

# ------------------------------------------------------------------------------
# 1. KARYOTYPE
# ------------------------------------------------------------------------------
karyotype <- read.delim(karyotype_file, sep = "\t", stringsAsFactors = FALSE)
if (!all(c("Chr", "Start", "End") %in% names(karyotype))) {
  stop(karyotype_file, " must have a header with columns Chr, Start, End; found: ",
       paste(names(karyotype), collapse = ", "))
}
karyotype <- karyotype %>%
  transmute(Chr = norm_chr(Chr), Start = as.numeric(Start), End = as.numeric(End)) %>%
  filter(!is.na(End))
karyotype <- karyotype[match(str_sort(karyotype$Chr, numeric = TRUE), karyotype$Chr), ]

if (nrow(karyotype) == 0) stop("No chromosomes read from ", karyotype_file)
if (anyDuplicated(karyotype$Chr)) stop("Duplicate chromosome names after renaming: ",
                                       paste(karyotype$Chr[duplicated(karyotype$Chr)], collapse = ", "))
message(sprintf("Karyotype: %d chromosomes (%s)",
                nrow(karyotype), paste(karyotype$Chr, collapse = ", ")))

# ------------------------------------------------------------------------------
# 2. LOAD FEATURES
# ------------------------------------------------------------------------------
mrna <- prep_features(read_gff(mrna_file), karyotype, mrna_feature, "mRNA features")
tes  <- prep_features(read_gff(te_file),   karyotype, te_feature,   "TE features")

lnc <- read.delim(lncRNA_file, sep = "\t", stringsAsFactors = FALSE)
tpm_cols <- grep("^mean_TPM_", names(lnc), value = TRUE)

if (expressed_only) {
  if (!length(tpm_cols)) stop("--expressed-only needs mean_TPM_* columns in ", lncRNA_file)
  mx <- apply(lnc[, tpm_cols, drop = FALSE], 1, max)
  lnc <- lnc[mx > LOWER_TPM & mx <= UPPER_TPM, ]
  message(sprintf("lncRNA filter: %g < max TPM <= %g -> %d transcripts retained", LOWER_TPM, UPPER_TPM, nrow(lnc)))
}

# Per-transcript expression (summary across the mean_TPM_* columns)
tpm_track <- length(tpm_cols) > 0
if (tpm_track) {
  tpm_mat <- as.matrix(lnc[, tpm_cols, drop = FALSE])
  lnc$tx_tpm <- if (tpm_stat == "mean") rowMeans(tpm_mat, na.rm = TRUE) else apply(tpm_mat, 1, max, na.rm = TRUE)
  message(sprintf("TPM track: %s of %d mean_TPM_* columns per transcript, highest isoform per gene%s",
                  tpm_stat, length(tpm_cols), if (log_tpm) ", log10(TPM + 1)" else ""))
} else {
  lnc$tx_tpm <- NA_real_
  message("No mean_TPM_* columns found; skipping the TPM track.")
}

chr_len <- setNames(karyotype$End, karyotype$Chr)

# One record per lncRNA GENE (isoforms merged: min start, max stop, most common Category)
lnc_genes <- lnc %>%
  mutate(Chr = norm_chr(Chr), Category = replace_na(Category, "Unclassified")) %>%
  filter(Chr %in% karyotype$Chr) %>%
  group_by(lncRNA_gene) %>%
  summarise(Chr   = first(Chr),
            Start = min(Start),
            End   = max(Stop),
            Category = names(sort(table(Category), decreasing = TRUE))[1],
            TPM   = if (all(is.na(tx_tpm))) NA_real_ else max(tx_tpm, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Start = pmax(Start, 1),
         End   = pmin(End, unname(chr_len[Chr])),
         TPM_plot = if (log_tpm) log10(TPM + 1) else TPM)

message(sprintf("Features on chromosomes: %d mRNA | %d TE records | %d lncRNA genes",
                nrow(mrna), nrow(tes), nrow(lnc_genes)))
print(table(Chr = lnc_genes$Chr, Category = lnc_genes$Category))

# ------------------------------------------------------------------------------
# 3. TRACK DATA
# ------------------------------------------------------------------------------
mrna_density <- window_density(mrna, karyotype, window_bp)
te_density   <- window_density(tes,  karyotype, window_bp)

write.table(mrna_density, file.path(out_dir, "density_mRNA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(te_density,   file.path(out_dir, "density_TE.tsv"),   sep = "\t", quote = FALSE, row.names = FALSE)

if (tpm_track) {
  lnc_tpm <- window_mean(lnc_genes, "TPM_plot", karyotype, tpm_window_bp)
  write.table(lnc_tpm, file.path(out_dir, "tpm_lncRNA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# lncRNA categories: colours fixed by alphabetical order (Okabe-Ito, colour-blind-safe),
# stacking order = most abundant category at the base of each bar
cat_levels <- sort(unique(lnc_genes$Category))
cat_cols   <- setNames(rep_len(c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                                 "#0072B2", "#D55E00", "#CC79A7", "#999999"), length(cat_levels)),
                       cat_levels)
if ("Unclassified" %in% cat_levels) cat_cols["Unclassified"] <- "#999999"
cat_order <- names(sort(table(lnc_genes$Category), decreasing = TRUE))

# stacked-count track (per window) and per-chromosome totals (centre of the plot)
lnc_stack <- window_category(lnc_genes, karyotype, stack_window_bp, cat_order)
stack_max <- max(rowSums(lnc_stack[, cat_order, drop = FALSE]), 1)
write.table(lnc_stack, file.path(out_dir, "lncRNA_category_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
chr_cat <- table(factor(lnc_genes$Chr, levels = karyotype$Chr),
                 factor(lnc_genes$Category, levels = cat_order))

# heatmap tracks: colour scale tops out at a high quantile of the non-zero windows
cap_of <- function(v) {
  v <- v[v > 0]
  if (!length(v)) return(1)
  max(quantile(v, heat_quantile, names = FALSE), 1)
}
tracks <- list(
  list(name = "mRNA", data = mrna_density, ramp = c("#f7f7f7", "#9ecae1", "#08519c")),
  list(name = "TE",   data = te_density,   ramp = c("#f7f7f7", "#fdae6b", "#d94801"))
)
for (i in seq_along(tracks)) tracks[[i]]$cap <- cap_of(tracks[[i]]$data$Value)

# ------------------------------------------------------------------------------
# 4. FIGURE
# ------------------------------------------------------------------------------
draw_circos <- function() {
  circos.clear()
  n <- nrow(karyotype)
  circos.par(start.degree = 90, gap.after = c(rep(3, n - 1), 10),
             cell.padding = c(0, 0, 0, 0), track.margin = c(0.004, 0.004),
             canvas.xlim = c(-1, 1.7), canvas.ylim = c(-1, 1))
  
  # outer track: chromosome axis + labels
  circos.genomicInitialize(as.data.frame(karyotype), plotType = c("axis", "labels"),
                           labels.cex = 0.9, axis.labels.cex = 0.5)
  first_sector <- karyotype$Chr[1]
  
  # 1-2: mRNA and TE density heatmaps
  for (tr in tracks) {
    col_fun <- colorRamp2(seq(0, tr$cap, length.out = 3), tr$ramp)
    circos.genomicTrack(as.data.frame(tr$data), ylim = c(0, 1), track.height = 0.06,
                        bg.border = NA,
                        panel.fun = function(region, value, ...) {
                          cols <- col_fun(value[[1]])
                          circos.genomicRect(region, value, ybottom = 0, ytop = 1,
                                             col = cols, border = cols, lwd = 0.1)
                        })
  }
  
  # 3: lncRNA genes per window, stacked counts by Category
  circos.genomicTrack(as.data.frame(lnc_stack), ylim = c(0, stack_max), track.height = 0.10,
                      bg.col = "#f5f5f5", bg.border = NA,
                      panel.fun = function(region, value, ...) {
                        cum <- rep(0, nrow(region))
                        for (j in seq_along(cat_order)) {
                          nxt <- cum + value[[j]]
                          circos.genomicRect(region, value, ybottom = cum, ytop = nxt,
                                             col = cat_cols[[cat_order[j]]], border = NA)
                          cum <- nxt
                        }
                        if (CELL_META$sector.index == first_sector) {
                          circos.yaxis(side = "left", at = c(0, stack_max), labels.cex = 0.35,
                                       tick.length = 0.5)
                        }
                      })
  
  # 4: lncRNA expression line
  tpm_max <- NA
  if (tpm_track) {
    tpm_max <- max(lnc_tpm$Value, 0.1)
    circos.genomicTrack(as.data.frame(lnc_tpm[, c("Chr", "Start", "End", "Value")]),
                        ylim = c(0, tpm_max), track.height = 0.10,
                        bg.col = "#f5f5f5", bg.border = NA,
                        panel.fun = function(region, value, ...) {
                          circos.genomicLines(region, value, type = "l", col = "#006d2c", lwd = 0.9)
                          if (CELL_META$sector.index == first_sector) {
                            circos.yaxis(side = "left", at = c(0, signif(tpm_max, 2)),
                                         labels.cex = 0.35, tick.length = 0.5)
                          }
                        })
  }
  
  # ---- centre: lncRNA gene counts per chromosome, stacked by Category ----
  r_in <- get.cell.meta.data("cell.bottom.radius", sector.index = first_sector,
                             track.index = get.current.track.index())
  n_chr <- nrow(karyotype)
  step  <- r_in / n_chr
  bh    <- step * 0.6
  xs    <- -0.4 * r_in
  wmax  <- 1.0 * r_in
  tot   <- rowSums(chr_cat)
  maxt  <- max(tot, 1)
  ycs   <- 0.5 * r_in - (seq_len(n_chr) - 0.5) * step
  text(0, 0.5 * r_in + 0.12 * r_in, "lncRNA genes per chromosome", cex = 0.55, font = 2, xpd = NA)
  for (i in seq_len(n_chr)) {
    x <- xs
    for (j in seq_along(cat_order)) {
      w <- chr_cat[i, j] / maxt * wmax
      if (w > 0) rect(x, ycs[i] - bh / 2, x + w, ycs[i] + bh / 2,
                      col = cat_cols[[cat_order[j]]], border = NA, xpd = NA)
      x <- x + w
    }
    text(xs - 0.02 * r_in, ycs[i], karyotype$Chr[i], adj = c(1, 0.5), cex = 0.6, xpd = NA)
    text(x + 0.02 * r_in,  ycs[i], tot[i],           adj = c(0, 0.5), cex = 0.55, xpd = NA)
  }
  segments(xs, ycs[n_chr] - bh, xs, ycs[1] + bh, col = "grey50", lwd = 0.6, xpd = NA)
  
  # ---- legend column (right of the circle) ----
  x0 <- 1.08; w <- 0.5
  y  <- 0.95
  for (i in seq_along(tracks)) {
    tr <- tracks[[i]]
    text(x0, y, sprintf("%d. %s per %s kb", i, tr$name, format(window_bp / 1e3, big.mark = ",")),
         adj = c(0, 0.5), cex = 0.7, xpd = NA)
    pal <- colorRampPalette(tr$ramp)(60)
    xs2 <- seq(x0, x0 + w, length.out = 61)
    rect(xs2[-61], y - 0.07, xs2[-1], y - 0.03, col = pal, border = NA, xpd = NA)
    text(x0,     y - 0.095, "0",                          adj = c(0, 0.5), cex = 0.55, xpd = NA)
    text(x0 + w, y - 0.095, paste0(">= ", signif(tr$cap, 3)), adj = c(1, 0.5), cex = 0.55, xpd = NA)
    y <- y - 0.20
  }
  k <- length(tracks) + 1
  leg <- legend(x = x0, y = y, legend = rev(cat_order), fill = unname(rev(cat_cols[cat_order])),
                title = sprintf("%d. lncRNA genes per %s kb (counts)", k,
                                format(stack_window_bp / 1e3, big.mark = ",")),
                bty = "n", cex = 0.65, xpd = NA, title.adj = 0)
  y <- leg$rect$top - leg$rect$h - 0.08
  k <- k + 1
  if (tpm_track) {
    text(x0, y, sprintf("%d. lncRNA %s (line)", k, if (log_tpm) "log10(TPM + 1)" else "TPM"),
         adj = c(0, 0.5), cex = 0.7, xpd = NA)
    text(x0, y - 0.05, sprintf("mean per gene in %s kb windows; axis 0-%s",
                               format(tpm_window_bp / 1e3, big.mark = ","), signif(tpm_max, 2)),
         adj = c(0, 0.5), cex = 0.5, xpd = NA)
  }
  
  circos.clear()
}

base <- file.path(out_dir, "circos_heatmap")
pdf(paste0(base, ".pdf"), width = 11, height = 8);                       draw_circos(); invisible(dev.off())
png(paste0(base, ".png"), width = 11, height = 8, units = "in", res = 300); draw_circos(); invisible(dev.off())
tryCatch({
  svg(paste0(base, ".svg"), width = 11, height = 8); draw_circos(); invisible(dev.off())
}, error = function(e) message("SVG output failed (needs cairo): ", conditionMessage(e)))

message("Done. Output written to: ", normalizePath(out_dir))