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
#   Interactive: edit the CONFIG block below and source the script.
#   Command line: Rscript lncRNA_ideogram.R karyotype.txt lncRNA.tsv mRNA.gff3 TE.gff3 [out_dir]
#
# install.packages(c("tidyverse", "RIdeogram", "rsvg"))   # rsvg: SVG -> PDF/PNG
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(RIdeogram)
})

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
setwd("~/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/chapter4/lncRNA_met_dev")

karyotype_file <- "A.tha/01_inputs/A.tha.karyotype.txt"                         # <- from datasets_cleanup.sb
lncRNA_file <- "A.tha/03_outputs/11_int_class/A.tha_final_annotation.tsv"
mrna_file   <- "A.tha/01_inputs/A.tha_mRNA.gtf"                                 # <- from datasets_cleanup.sb
te_file     <- "A.tha/annotation/A.tha_TE.gff3"                                 # <- edit
out_dir     <- "A.tha/03_outputs/13_ideogram"

# Command-line arguments override the paths above
cli <- commandArgs(trailingOnly = TRUE)
if (length(cli) >= 4) {
  karyotype_file <- cli[1]
  lncRNA_file <- cli[2]
  mrna_file   <- cli[3]
  te_file     <- cli[4]
  if (length(cli) >= 5) out_dir <- cli[5]
}

window_bp     <- 500000          # window size for mRNA / TE / lncRNA density heatmaps
mrna_feature  <- "^(mRNA|transcript)$"   # regex matched against GFF/GTF column 3 for the mRNA file
                                         # (AGAT's GTF writes mRNAs as 'transcript'; the gene lines are excluded)
te_feature    <- NULL            # regex for GFF column 3 in the TE file; NULL = use every row

# Optional: restrict the lncRNAs to the expressed set used in the TPM figures
expressed_only <- FALSE
LOWER_TPM <- 0.1
UPPER_TPM <- 50

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