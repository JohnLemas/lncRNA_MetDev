#!/usr/bin/env Rscript
# ==============================================================================
# lncRNA_circos.R
#
# Circos plot of mRNA, TE and lncRNA features (one track each) using circlize.
#
# Track layout (outside -> inside)
#   0. chromosome axis + labels
#   1. mRNA density
#   2. TE density
#   3. lncRNA gene density
#   4. lncRNA genes as ticks coloured by Category  (skip with --no-category-track)
#
# Inputs  (same as lncRNA_ideogram.R)
#   1. karyotype file     (Chr, Start, End with a header; nuclear chromosomes only)
#   2. lncRNA annotation  (A.tha_final_annotation.tsv; needs Chr, Start, Stop,
#                          lncRNA_gene, Category and optionally mean_TPM_*)
#   3. mRNA annotation    (GFF3/GTF-style, 9 tab-separated columns)
#   4. TE annotation      (GFF3/GTF-style, 9 tab-separated columns)
#
# Outputs (in out_dir)
#   density_*.tsv                  per-window feature counts used for the tracks
#   circos_<style>.svg / .pdf / .png
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
              help = "window size (bp) for the density tracks [default: %default]"),
  make_option("--style", type = "character", default = "area",
              help = "density track style: 'area' (filled line) or 'heatmap' [default: %default]"),
  make_option("--mark-width", dest = "mark_width", type = "integer", default = 30000L,
              help = "minimum drawn width (bp) of lncRNA ticks in the category track [default: %default]"),
  make_option("--no-category-track", dest = "no_category_track", action = "store_true", default = FALSE,
              help = "omit the innermost lncRNA-by-Category track"),
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

karyotype_file <- opt$karyotype
lncRNA_file    <- opt$lncrna
mrna_file      <- opt$mrna
te_file        <- opt$te
out_dir        <- opt$outdir
window_bp      <- opt$window
style          <- match.arg(opt$style, c("area", "heatmap"))
mark_width     <- opt$mark_width
category_track <- !opt$no_category_track
mrna_feature   <- opt$mrna_feature
te_feature     <- opt$te_feature
expressed_only <- opt$expressed_only
LOWER_TPM      <- opt$lower_tpm
UPPER_TPM      <- opt$upper_tpm

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
# 3. DENSITY TRACKS (counts per window)
# ------------------------------------------------------------------------------
mrna_density <- window_density(mrna,      karyotype, window_bp)
te_density   <- window_density(tes,       karyotype, window_bp)
lnc_density  <- window_density(lnc_genes, karyotype, window_bp)

write.table(mrna_density, file.path(out_dir, "density_mRNA.tsv"),   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(te_density,   file.path(out_dir, "density_TE.tsv"),     sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lnc_density,  file.path(out_dir, "density_lncRNA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# lncRNA category colours (Okabe-Ito, colour-blind-safe), WITH '#'
cat_levels <- sort(unique(lnc_genes$Category))
cat_cols   <- setNames(rep_len(c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                                 "#0072B2", "#D55E00", "#CC79A7", "#999999"), length(cat_levels)),
                       cat_levels)
if ("Unclassified" %in% cat_levels) cat_cols["Unclassified"] <- "#999999"

# ticks need a minimum width to be visible at genome scale
lnc_ticks <- lnc_genes %>%
  mutate(mid = (Start + End) / 2,
         S = pmax(1, pmin(Start, mid - mark_width / 2)),
         E = pmin(unname(chr_len[Chr]), pmax(End, mid + mark_width / 2))) %>%
  transmute(Chr, Start = S, End = E, col = unname(cat_cols[Category]))

# track definitions: data, fill colour, heatmap ramp
tracks <- list(
  list(name = "mRNA",  data = mrna_density, fill = "#3182bd", ramp = c("#f7f7f7", "#9ecae1", "#08519c")),
  list(name = "TE",    data = te_density,   fill = "#e6550d", ramp = c("#f7f7f7", "#fdae6b", "#d94801")),
  list(name = "lncRNA", data = lnc_density, fill = "#31a354", ramp = c("#f7f7f7", "#a1d99b", "#006d2c"))
)

# ------------------------------------------------------------------------------
# 4. FIGURE
# ------------------------------------------------------------------------------
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
  
  for (tr in tracks) {
    d    <- as.data.frame(tr$data)
    ymax <- max(d$Value, 1)
    
    if (style == "area") {
      circos.genomicTrack(d, ylim = c(0, ymax), track.height = 0.09,
                          bg.col = "#f5f5f5", bg.border = NA,
                          panel.fun = function(region, value, ...) {
                            circos.genomicLines(region, value, area = TRUE, baseline = 0,
                                                col = adjustcolor(tr$fill, 0.55),
                                                border = tr$fill, lwd = 0.6)
                            if (CELL_META$sector.index == first_sector) {
                              circos.yaxis(side = "left", at = c(0, ymax), labels.cex = 0.35,
                                           tick.length = 0.5)
                            }
                          })
    } else {
      col_fun <- colorRamp2(seq(0, ymax, length.out = 3), tr$ramp)
      circos.genomicTrack(d, ylim = c(0, 1), track.height = 0.06, bg.border = NA,
                          panel.fun = function(region, value, ...) {
                            circos.genomicRect(region, value, ybottom = 0, ytop = 1,
                                               col = col_fun(value[[1]]), border = NA)
                          })
    }
  }
  
  if (category_track) {
    circos.genomicTrack(as.data.frame(lnc_ticks), ylim = c(0, 1), track.height = 0.07,
                        bg.col = "#f5f5f5", bg.border = NA,
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, ybottom = 0, ytop = 1,
                                             col = value$col, border = value$col, lwd = 0.3)
                        })
  }
  
  # legend (right of the circle)
  legend_items <- c("Tracks, outer to inner:",
                    paste0(seq_along(tracks), ". ", vapply(tracks, `[[`, "", "name"),
                           " per ", format(window_bp / 1e3, big.mark = ","), " kb"),
                    if (category_track) paste0(length(tracks) + 1, ". lncRNA genes by Category"))
  legend(x = 1.05, y = 0.95, legend = legend_items, bty = "n", cex = 0.75, xpd = NA)
  legend(x = 1.05, y = 0.35, legend = names(cat_cols), fill = unname(cat_cols),
         title = "lncRNA Category", bty = "n", cex = 0.75, xpd = NA, title.adj = 0)
  
  circos.clear()
}

base <- file.path(out_dir, paste0("circos_", style))
pdf(paste0(base, ".pdf"), width = 10, height = 7.7);           draw_circos(); invisible(dev.off())
png(paste0(base, ".png"), width = 10, height = 7.7, units = "in", res = 300); draw_circos(); invisible(dev.off())
tryCatch({
  svg(paste0(base, ".svg"), width = 10, height = 7.7); draw_circos(); invisible(dev.off())
}, error = function(e) message("SVG output failed (needs cairo): ", conditionMessage(e)))

message("Done. Output written to: ", normalizePath(out_dir))
