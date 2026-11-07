#!/usr/bin/env Rscript
# ==============================================================================
# lncRNA_circos.R
#
# Circos plot of mRNA, TE and lncRNA features using circlize.
# Each feature type is a HEATMAP track; the centre of the circle holds a
# stacked-bar chart of lncRNA classifications per chromosome.
#
# Track layout (outside -> inside)
#   0. chromosome axis + labels
#   1. mRNA heatmap     (genes per window)
#   2. TE heatmap       (% of window covered by TEs, or TE count per window)
#   3. lncRNA heatmap   (lncRNA genes per window)
#   4. (optional, --category-track) lncRNA genes as ticks coloured by Category
#   centre: stacked bars of lncRNA Category per chromosome
#
# Inputs
#   1. karyotype file     (Chr, Start, End with a header; nuclear chromosomes only)
#   2. lncRNA annotation  (needs Chr, Start, Stop, lncRNA_gene, Category and
#                          optionally mean_TPM_*)
#   3. mRNA annotation    (GFF3/GTF-style, 9 tab-separated columns)
#   4. TE annotation      (GFF3/GTF-style, 9 tab-separated columns; EDTA
#                          *.TEanno.gff3 or *.intact.gff3 both work)
#
# Outputs (in out_dir)
#   density_mRNA.tsv, density_TE.tsv, density_lncRNA.tsv   per-window values
#       (Value = what is plotted; density_TE.tsv also has Count and Coverage)
#   circos_heatmap.svg / .pdf / .png
#
# Usage
#   Rscript lncRNA_circos.R \
#     --karyotype A.tha/01_inputs/A.tha.karyotype.txt \
#     --lncrna    A.tha/03_outputs/11_int_class/A.tha_final_annotation.tsv \
#     --mrna      A.tha/01_inputs/A.tha_mRNA.gtf \
#     --te        A.tha/01_inputs/A.tha_TE.gff3 \
#     --outdir    A.tha/03_outputs/13_circos
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
  make_option(c("-w", "--window"), type = "integer", default = 500000L,
              help = "window size (bp) for the heatmap tracks [default: %default]"),
  make_option("--te-metric", dest = "te_metric", type = "character", default = "coverage",
              help = paste("what the TE heatmap shows: 'coverage' (fraction of each window covered",
                           "by TEs; overlaps merged) or 'count' (TE rows per window) [default: %default]")),
  make_option("--clip-quantile", dest = "clip_quantile", type = "double", default = 0.99,
              help = paste("heatmap colour scale tops out at this quantile of the window values, so a",
                           "single extreme window cannot wash out the rest; 1 = no clipping [default: %default]")),
  make_option("--center-chart", dest = "center_chart", type = "character", default = "counts",
              help = paste("centre chart of lncRNA Category per chromosome: 'counts' (stacked gene",
                           "counts), 'proportion' (100%% stacked) or 'none' [default: %default]")),
  make_option("--category-track", dest = "category_track", action = "store_true", default = FALSE,
              help = "also draw the innermost per-gene lncRNA Category tick track [default: off]"),
  make_option("--mark-width", dest = "mark_width", type = "integer", default = 30000L,
              help = "with --category-track: minimum drawn width (bp) of lncRNA ticks [default: %default]"),
  make_option("--mrna-feature", dest = "mrna_feature", type = "character", default = "^(mRNA|transcript)$",
              help = paste("regex for GFF/GTF column 3 in the mRNA file; AGAT's GTF writes",
                           "mRNAs as 'transcript' [default: %default]")),
  make_option("--te-feature", dest = "te_feature", type = "character", default = NULL,
              help = paste("regex for GFF/GTF column 3 in the TE file, e.g. 'Gypsy|Copia' to plot",
                           "only those superfamilies [default: use every row not excluded below]")),
  make_option("--te-exclude", dest = "te_exclude", type = "character",
              default = "^(repeat_region|long_terminal_repeat|target_site_duplication|terminal_inverted_repeat)$",
              help = paste("regex for TE sub-feature rows to drop (EDTA writes one parent row plus LTR/TSD/TIR",
                           "child rows per element, which inflates counts); 'none' keeps every row [default: %default]")),
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
               description = "Circos plot of mRNA, TE and lncRNA heatmaps (circlize)."),
  args = if (interactive()) character(0) else commandArgs(trailingOnly = TRUE)
)

if (is.null(opt$workdir) && interactive()) {
  opt$workdir <- "~/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/chapter4/lncRNA_met_dev"
}
if (!is.null(opt$workdir)) setwd(opt$workdir)

karyotype_file <- opt$karyotype
lncRNA_file    <- opt$lncrna
mrna_file      <- opt$mrna
te_file        <- opt$te
out_dir        <- opt$outdir
window_bp      <- opt$window
te_metric      <- match.arg(opt$te_metric, c("coverage", "count"))
clip_q         <- opt$clip_quantile
center_chart   <- match.arg(opt$center_chart, c("counts", "proportion", "none"))
category_track <- opt$category_track
mark_width     <- opt$mark_width
mrna_feature   <- opt$mrna_feature
te_feature     <- opt$te_feature
te_exclude     <- if (is.null(opt$te_exclude) || tolower(opt$te_exclude) %in% c("none", "")) NULL else opt$te_exclude
expressed_only <- opt$expressed_only
LOWER_TPM      <- opt$lower_tpm
UPPER_TPM      <- opt$upper_tpm

if (!is.finite(clip_q) || clip_q <= 0 || clip_q > 1) stop("--clip-quantile must be in (0, 1]")

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

prep_features <- function(df, karyotype, include = NULL, exclude = NULL, label = "features") {
  available <- sort(unique(df$feature))
  if (!is.null(include)) {
    df <- filter(df, grepl(include, feature))
    if (nrow(df) == 0) {
      stop("No ", label, " matched feature regex '", include,
           "'. Feature types present: ", paste(available, collapse = ", "))
    }
  }
  if (!is.null(exclude)) {
    df <- filter(df, !grepl(exclude, feature))
    if (nrow(df) == 0) {
      stop("Every ", label, " row was removed by the exclude regex '", exclude,
           "'. Feature types present: ", paste(available, collapse = ", "))
    }
  }
  df <- df %>%
    mutate(Chr = norm_chr(seqid), Start = pmin(start, end), End = pmax(start, end)) %>%
    filter(Chr %in% karyotype$Chr)
  if (nrow(df) == 0) stop("No ", label, " fall on the chromosomes in the karyotype; check sequence names.")
  df
}

window_table <- function(len, window_bp) {
  starts <- seq(1, len, by = window_bp)
  tibble(Start = starts, End = pmin(starts + window_bp - 1, len))
}

# Counts of features per window (feature assigned to the window containing its midpoint)
window_counts <- function(feat, karyotype, window_bp) {
  map_dfr(seq_len(nrow(karyotype)), function(i) {
    chr <- karyotype$Chr[i]
    w   <- window_table(karyotype$End[i], window_bp)
    f   <- feat[feat$Chr == chr, ]
    idx <- findInterval((f$Start + f$End) / 2, w$Start)
    idx <- idx[idx >= 1 & idx <= nrow(w)]
    tibble(Chr = chr, Start = w$Start, End = w$End,
           Count = tabulate(idx, nbins = nrow(w)))
  })
}

# Merge overlapping / adjacent intervals so bases are not counted twice
merge_intervals <- function(s, e) {
  if (length(s) == 0) return(list(s = numeric(0), e = numeric(0)))
  o  <- order(s); s <- s[o]; e <- e[o]
  ce <- cummax(e)
  new_grp <- c(TRUE, s[-1] > ce[-length(ce)] + 1)
  last <- c(which(new_grp)[-1] - 1L, length(s))
  list(s = s[new_grp], e = ce[last])
}

# Fraction of each window covered by at least one feature
window_coverage <- function(feat, karyotype, window_bp) {
  map_dfr(seq_len(nrow(karyotype)), function(i) {
    chr <- karyotype$Chr[i]
    len <- karyotype$End[i]
    w   <- window_table(len, window_bp)
    f   <- feat[feat$Chr == chr, ]
    m   <- merge_intervals(pmax(f$Start, 1), pmin(f$End, len))
    lo  <- findInterval(w$Start - 1, m$e) + 1L   # first merged interval ending inside/after window
    hi  <- findInterval(w$End, m$s)              # last merged interval starting inside/before window
    bp  <- vapply(seq_len(nrow(w)), function(j) {
      if (lo[j] > hi[j]) return(0)
      k <- lo[j]:hi[j]
      sum(pmin(m$e[k], w$End[j]) - pmax(m$s[k], w$Start[j]) + 1)
    }, numeric(1))
    tibble(Chr = chr, Start = w$Start, End = w$End,
           Coverage = bp / (w$End - w$Start + 1))
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
mrna <- prep_features(read_gff(mrna_file), karyotype, include = mrna_feature, label = "mRNA features")

te_raw <- read_gff(te_file)
message("TE feature types in ", basename(te_file), " (column 3):")
print(sort(table(te_raw$feature), decreasing = TRUE))
tes <- prep_features(te_raw, karyotype, include = te_feature, exclude = te_exclude, label = "TE features")

lnc <- read.delim(lncRNA_file, sep = "\t", stringsAsFactors = FALSE)
if (expressed_only) {
  tpm_cols <- grep("^mean_TPM_", names(lnc), value = TRUE)
  mx <- apply(lnc[, tpm_cols], 1, max)
  lnc <- lnc[mx > LOWER_TPM & mx <= UPPER_TPM, ]
  message(sprintf("lncRNA filter: %g < max TPM <= %g -> %d transcripts retained", LOWER_TPM, UPPER_TPM, nrow(lnc)))
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
            .groups = "drop") %>%
  mutate(Start = pmax(Start, 1),
         End   = pmin(End, unname(chr_len[Chr])))

message(sprintf("Features on chromosomes: %d mRNA | %d TE records | %d lncRNA genes",
                nrow(mrna), nrow(tes), nrow(lnc_genes)))
print(table(Chr = lnc_genes$Chr, Category = lnc_genes$Category))

# ------------------------------------------------------------------------------
# 3. WINDOW VALUES (what the heatmaps show)
# ------------------------------------------------------------------------------
mrna_density <- window_counts(mrna, karyotype, window_bp) %>% mutate(Value = Count)
lnc_density  <- window_counts(lnc_genes, karyotype, window_bp) %>% mutate(Value = Count)
te_density   <- window_counts(tes, karyotype, window_bp) %>%
  left_join(window_coverage(tes, karyotype, window_bp), by = c("Chr", "Start", "End")) %>%
  mutate(Value = if (te_metric == "coverage") Coverage else Count)

write.table(mrna_density, file.path(out_dir, "density_mRNA.tsv"),   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(te_density,   file.path(out_dir, "density_TE.tsv"),     sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lnc_density,  file.path(out_dir, "density_lncRNA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# lncRNA category colours (Okabe-Ito, colour-blind-safe), WITH '#'
cat_levels <- sort(unique(lnc_genes$Category))
cat_cols   <- setNames(rep_len(c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                                 "#0072B2", "#D55E00", "#CC79A7", "#999999"), length(cat_levels)),
                       cat_levels)
if ("Unclassified" %in% cat_levels) cat_cols["Unclassified"] <- "#999999"

# ticks (optional track 4) need a minimum width to be visible at genome scale
lnc_ticks <- lnc_genes %>%
  mutate(mid = (Start + End) / 2,
         S = pmax(1, pmin(Start, mid - mark_width / 2)),
         E = pmin(unname(chr_len[Chr]), pmax(End, mid + mark_width / 2))) %>%
  transmute(Chr, Start = S, End = E, col = unname(cat_cols[Category]))

# lncRNA class x chromosome table for the centre chart
class_tab <- table(factor(lnc_genes$Chr, levels = karyotype$Chr),
                   factor(lnc_genes$Category, levels = cat_levels))

# heatmap tracks: data, colour ramp, and a colour scale clipped at clip_q
wkb <- format(window_bp / 1e3, big.mark = ",")
tracks <- list(
  list(name = "mRNA", title = sprintf("1. mRNA genes per %s kb", wkb), data = mrna_density,
       ramp = c("#f7f7f7", "#9ecae1", "#08519c"), pct = FALSE),
  list(name = "TE",
       title = if (te_metric == "coverage") sprintf("2. TE coverage (%% of %s kb window)", wkb)
               else sprintf("2. TE records per %s kb", wkb),
       data = te_density, ramp = c("#f7f7f7", "#fdae6b", "#d94801"), pct = (te_metric == "coverage")),
  list(name = "lncRNA", title = sprintf("3. lncRNA genes per %s kb", wkb), data = lnc_density,
       ramp = c("#f7f7f7", "#a1d99b", "#006d2c"), pct = FALSE)
)
tracks <- lapply(tracks, function(tr) {
  v <- tr$data$Value
  vmax <- as.numeric(quantile(v, clip_q, na.rm = TRUE))
  if (!is.finite(vmax) || vmax <= 0) vmax <- if (max(v, na.rm = TRUE) > 0) max(v, na.rm = TRUE) else 1
  tr$vmax    <- vmax
  tr$clipped <- any(v > vmax, na.rm = TRUE)
  tr$col_fun <- colorRamp2(c(0, vmax / 2, vmax), tr$ramp)
  tr$fmt     <- if (tr$pct) function(x) paste0(round(x * 100), "%")
                else function(x) format(round(x, 1), big.mark = ",", trim = TRUE)
  tr
})

# ------------------------------------------------------------------------------
# 4. FIGURE
# ------------------------------------------------------------------------------

# horizontal colour bar drawn in plot coordinates
draw_colorbar <- function(x, y, w, h, tr) {
  n  <- 60
  xs <- seq(x, x + w, length.out = n + 1)
  rect(xs[-(n + 1)], y, xs[-1], y + h, col = tr$col_fun(seq(0, tr$vmax, length.out = n)),
       border = NA, xpd = NA)
  rect(x, y, x + w, y + h, border = "grey40", lwd = 0.4, xpd = NA)
  text(x, y + h + 0.02, tr$title, adj = c(0, 0), cex = 0.62, xpd = NA)
  text(c(x, x + w), y - 0.015,
       c(tr$fmt(0), paste0(tr$fmt(tr$vmax), if (tr$clipped) "+" else "")),
       adj = c(0.5, 1), cex = 0.5, xpd = NA)
}

# stacked bars of lncRNA Category per chromosome, fitted inside radius r_in
draw_center_chart <- function(r_in) {
  prop <- (center_chart == "proportion")
  tot  <- rowSums(class_tab)
  m    <- if (prop) sweep(class_tab, 1, pmax(tot, 1), "/") else class_tab
  xmax <- if (prop) 1 else max(tot, 1)
  
  r <- max(r_in - 0.02, 0.2)
  w <- 0.66 * r; h <- 0.68 * r          # half-width / half-height of the chart box
  W <- 2 * w;    H <- 2 * h
  x0 <- -w + 0.20 * W                   # room for chromosome labels on the left
  x1 <-  w - 0.14 * W                   # room for totals on the right
  y_top <- h - 0.17 * H                 # room for the title
  y_bot <- -h + 0.17 * H                # room for the axis
  nr    <- nrow(m)
  rowh  <- (y_top - y_bot) / nr
  bh    <- 0.72 * rowh
  cex_l <- if (nr > 8) 0.5 else 0.7
  x_of  <- function(v) x0 + v / xmax * (x1 - x0)
  
  text(0, h, "lncRNA classes per chromosome", font = 2, cex = 0.7, adj = c(0.5, 1), xpd = NA)
  
  for (i in seq_len(nr)) {
    yc   <- y_top - (i - 0.5) * rowh
    vals <- as.numeric(m[i, ])
    keep <- vals > 0
    cum  <- c(0, cumsum(vals))
    text(x0 - 0.012, yc, rownames(m)[i], adj = c(1, 0.5), cex = cex_l, xpd = NA)
    if (any(keep)) {
      rect(x_of(cum[-length(cum)][keep]), yc - bh / 2, x_of(cum[-1][keep]), yc + bh / 2,
           col = cat_cols[cat_levels][keep], border = "white", lwd = 0.2, xpd = NA)
    }
    text(x_of(sum(vals)) + 0.01, yc, format(tot[i], big.mark = ","),
         adj = c(0, 0.5), cex = cex_l * 0.9, col = "grey25", xpd = NA)
  }
  
  # axis
  tk <- if (prop) c(0, 0.5, 1) else pretty(c(0, xmax), n = 3)
  tk <- tk[tk <= xmax]
  segments(x0, y_bot, x1, y_bot, lwd = 0.5, col = "grey30", xpd = NA)
  segments(x_of(tk), y_bot, x_of(tk), y_bot - 0.012, lwd = 0.5, col = "grey30", xpd = NA)
  text(x_of(tk), y_bot - 0.016,
       if (prop) paste0(tk * 100, "%") else format(tk, big.mark = ",", trim = TRUE),
       adj = c(0.5, 1), cex = 0.45, xpd = NA)
  text((x0 + x1) / 2, -h, if (prop) "% of lncRNA genes (n at right)" else "lncRNA genes (total at right)",
       adj = c(0.5, 0), cex = 0.5, xpd = NA)
}

draw_circos <- function() {
  circos.clear()
  n <- nrow(karyotype)
  circos.par(start.degree = 90, gap.after = c(rep(3, n - 1), 10),
             cell.padding = c(0, 0, 0, 0), track.margin = c(0.004, 0.004),
             canvas.xlim = c(-1, 1.6), canvas.ylim = c(-1, 1))
  
  # outer track: chromosome axis + labels
  circos.genomicInitialize(as.data.frame(karyotype), plotType = c("axis", "labels"),
                           labels.cex = 0.9, axis.labels.cex = 0.5)
  first_sector <- karyotype$Chr[1]
  
  # heatmap tracks
  for (tr in tracks) {
    d       <- as.data.frame(tr$data[, c("Chr", "Start", "End", "Value")])
    col_fun <- tr$col_fun
    circos.genomicTrack(d, ylim = c(0, 1), track.height = 0.07, bg.border = NA,
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, ybottom = 0, ytop = 1,
                                             col = col_fun(value[[1]]), border = NA)
                        })
  }
  
  if (category_track) {
    circos.genomicTrack(as.data.frame(lnc_ticks), ylim = c(0, 1), track.height = 0.07,
                        bg.col = "#f5f5f5", bg.border = NA,
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, ybottom = 0, ytop = 1,
                                             col = value$col, border = value$col, lwd = 0.3)
                        })
  }
  
  # centre chart, sized to the free space inside the innermost track
  if (center_chart != "none") {
    r_in <- tryCatch(get.cell.meta.data("cell.bottom.radius", sector.index = first_sector,
                                        track.index = circos.par("track.index")),
                     error = function(e) NA_real_)
    if (length(r_in) != 1 || !is.finite(r_in)) r_in <- 0.6
    draw_center_chart(r_in)
  }
  
  # colour bars + class legend (right of the circle)
  bar_y <- c(0.80, 0.58, 0.36)
  for (i in seq_along(tracks)) draw_colorbar(1.08, bar_y[i], 0.42, 0.04, tracks[[i]])
  legend_title <- paste0("lncRNA class",
                         if (center_chart != "none" && category_track) " (centre chart, track 4)"
                         else if (center_chart != "none") " (centre chart)"
                         else if (category_track) " (track 4)" else "")
  legend(x = 1.08, y = 0.20, legend = names(cat_cols), fill = unname(cat_cols),
         title = legend_title, bty = "n", cex = 0.7, xpd = NA, title.adj = 0)
  
  circos.clear()
}

base <- file.path(out_dir, "circos_heatmap")
pdf(paste0(base, ".pdf"), width = 10, height = 7.7);           draw_circos(); invisible(dev.off())
png(paste0(base, ".png"), width = 10, height = 7.7, units = "in", res = 300); draw_circos(); invisible(dev.off())
tryCatch({
  svg(paste0(base, ".svg"), width = 10, height = 7.7); draw_circos(); invisible(dev.off())
}, error = function(e) message("SVG output failed (needs cairo): ", conditionMessage(e)))

message("Done. Output written to: ", normalizePath(out_dir))