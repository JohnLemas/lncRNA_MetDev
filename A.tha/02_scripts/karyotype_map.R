#!/usr/bin/env Rscript
# ==============================================================================
# lncRNA_ideogram.R
#
# Maps lncRNA, mRNA and TE locations onto a genome with RIdeogram.
#
# Inputs
#   1. karyotype file          (<SPEC>.karyotype.txt from datasets_cleanup.sb: tab-separated
#                               Chr, Start, End with a header; nuclear chromosomes only)
#   2. lncRNA annotation       (the A.tha_final_annotation.tsv table; needs columns
#                               Chr, Start, Stop, lncRNA_gene, Category and,
#                               optionally, mean_TPM_* for the expression filter)
#   3. mRNA annotation         (GFF3/GTF-style, 9 tab-separated columns)
#   4. TE annotation           (GFF3/GTF-style, 9 tab-separated columns)
#
# Outputs (in out_dir)
#   density_*.tsv                   per-window feature counts used for the heatmaps
#   ideogram_*.svg / .pdf / .png    the figures
#
# Usage
#   Command line (every flag is optional; defaults are in the CONFIG block below):
#     Rscript lncRNA_ideogram.R \
#       --karyotype A.tha/01_inputs/A.tha.karyotype.txt \
#       --lncrna    A.tha/03_outputs/11_int_class/A.tha_final_annotation.tsv \
#       --mrna      A.tha/01_inputs/A.tha_mRNA.gtf \
#       --te        A.tha/annotation/A.tha_TE.gff3 \
#       --outdir    A.tha/03_outputs/13_ideogram
#     Rscript lncRNA_ideogram.R --help        # all flags and their defaults
#   Interactive: edit the defaults in the CONFIG block and source the script.
#
# install.packages(c("tidyverse", "RIdeogram", "optparse", "rsvg"))   # rsvg: SVG -> PDF/PNG
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(RIdeogram)
  library(optparse)
})

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
# Defaults live here; any flag given on the command line overrides them.
option_list <- list(
  make_option(c("-k", "--karyotype"), type = "character",
              default = "A.tha/01_inputs/A.tha.karyotype.txt",
              help = "karyotype file from datasets_cleanup.sb (Chr, Start, End) [default: %default]"),
  make_option(c("-l", "--lncrna"), type = "character",
              default = "A.tha/03_outputs/11_int_class/A.tha_final_annotation.tsv",
              help = "lncRNA annotation table [default: %default]"),
  make_option(c("-m", "--mrna"), type = "character",
              default = "A.tha/01_inputs/A.tha_mRNA.gtf",
              help = "mRNA annotation, GFF3/GTF [default: %default]"),
  make_option(c("-t", "--te"), type = "character",
              default = "A.tha/annotation/A.tha_TE.gff3",
              help = "TE annotation, GFF3/GTF [default: %default]"),
  make_option(c("-o", "--outdir"), type = "character",
              default = "A.tha/03_outputs/13_ideogram",
              help = "output directory [default: %default]"),
  make_option(c("-w", "--window"), type = "integer", default = 500000L,
              help = "window size (bp) for the density heatmaps [default: %default]"),
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
               description = "Map lncRNA, mRNA and TE locations onto a genome with RIdeogram."),
  args = if (interactive()) character(0) else commandArgs(trailingOnly = TRUE)
)

# Interactive sessions (source / RStudio) start from the dissertation folder
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
mrna_feature   <- opt$mrna_feature
te_feature     <- opt$te_feature       # NULL = use every row
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

# Make sequence names comparable between the karyotype, the GFF/GTFs and the lncRNA table:
# "1", "chr1" and "Chr1" all become "Chr1".
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
  filter(out, !is.na(start), !is.na(end))   # drops stray ##FASTA sequence lines etc.
}

# Rename chromosomes, keep only those in the karyotype, optionally filter by feature type
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

# Draw one ideogram, then convert SVG -> PDF and PNG if rsvg is installed
make_ideogram <- function(name, overlaid, label = NULL, label_type = NULL,
                          colorset1, colorset2 = NULL) {
  svg_file <- file.path(out_dir, paste0(name, ".svg"))
  args <- list(karyotype = as.data.frame(karyotype), overlaid = as.data.frame(overlaid),
               colorset1 = colorset1, output = svg_file)
  if (!is.null(label)) {
    args$label      <- as.data.frame(label)
    args$label_type <- label_type
  }
  if (!is.null(colorset2)) args$colorset2 <- colorset2
  do.call(ideogram, args)
  
  tryCatch({
    convertSVG(svg_file, file = sub("\\.svg$", ".pdf", svg_file), device = "pdf")
    convertSVG(svg_file, file = sub("\\.svg$", ".png", svg_file), device = "png", dpi = 300)
  }, error = function(e) {
    message("SVG written but PDF/PNG conversion failed (is 'rsvg' installed?): ", conditionMessage(e))
  })
  invisible(svg_file)
}

# ------------------------------------------------------------------------------
# 1. KARYOTYPE (written by datasets_cleanup.sb; already nuclear chromosomes only)
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

write.table(mrna_density, file.path(out_dir, "density_mRNA.tsv"),  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(te_density,   file.path(out_dir, "density_TE.tsv"),    sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lnc_density,  file.path(out_dir, "density_lncRNA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# lncRNA markers coloured by Category (Okabe-Ito colour-blind-safe palette; hex WITHOUT '#')
cat_levels <- sort(unique(lnc_genes$Category))
cat_cols   <- setNames(rep_len(c("E69F00", "56B4E9", "009E73", "F0E442",
                                 "0072B2", "D55E00", "CC79A7", "999999"), length(cat_levels)),
                       cat_levels)
if ("Unclassified" %in% cat_levels) cat_cols["Unclassified"] <- "999999"

lnc_markers <- lnc_genes %>%
  transmute(Type = Category, Shape = "circle", Chr, Start, End,
            color = unname(cat_cols[Category]))

# ------------------------------------------------------------------------------
# 4. FIGURES
# ------------------------------------------------------------------------------
heat_mrna <- c("#f7f7f7", "#9ecae1", "#08519c")   # blue
heat_te   <- c("#f7f7f7", "#fdae6b", "#d94801")   # orange
heat_lnc  <- c("#f7f7f7", "#a1d99b", "#006d2c")   # green

# 1) mRNA density with lncRNA genes as markers (coloured by Category)
make_ideogram("ideogram_mRNA_density_lncRNA_markers",
              overlaid = mrna_density, label = lnc_markers, label_type = "marker",
              colorset1 = heat_mrna)

# 2) TE density with lncRNA genes as markers
make_ideogram("ideogram_TE_density_lncRNA_markers",
              overlaid = te_density, label = lnc_markers, label_type = "marker",
              colorset1 = heat_te)

# 3) mRNA density (on the chromosome) beside TE density (outer track)
make_ideogram("ideogram_mRNA_vs_TE_density",
              overlaid = mrna_density, label = te_density, label_type = "heatmap",
              colorset1 = heat_mrna, colorset2 = heat_te)

# 4) mRNA density beside lncRNA gene density
make_ideogram("ideogram_mRNA_vs_lncRNA_density",
              overlaid = mrna_density, label = lnc_density, label_type = "heatmap",
              colorset1 = heat_mrna, colorset2 = heat_lnc)

message("Done. Output written to: ", normalizePath(out_dir))