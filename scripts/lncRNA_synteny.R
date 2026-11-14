#!/usr/bin/env Rscript
# ==============================================================================
# 21_lncRNA_synteny.R
#
# Positionally conserved lncRNAs between two species, shown on a RIdeogram
# synteny plot.
#
# Inputs (conventions follow the rest of the pipeline; every path can be overridden)
#   <base>/synteny/<sp1>_vs_<sp2>/<code1>_<code2>.collinearity   MCScanX blocks (20_synteny_mcscanx.sb)
#   <base>/synteny/<sp1>_vs_<sp2>/<code1>_<code2>.gff            MCScanX gene coordinates
#   <base>/<sp>/01_inputs/<sp>.karyotype.txt                     karyotype (Chr, Start, End)
#   <base>/<sp>/03_outputs/11_int_class/<sp>_final_annotation.tsv   lncRNA table
#                                       (Chr, Start, Stop, lncRNA_gene, Category)
#   optional --seq-hits: 2-column TSV (lncRNA gene in sp1, lncRNA gene in sp2) from a
#                        sequence search; adds a sequence-support flag only
#
# Method
#   An lncRNA is "positionally conserved" when it sits between two anchor (collinear
#   ortholog) genes in its own species and an lncRNA of the other species sits between
#   the partners of those same two anchors. Done in both directions; pairs found from
#   both sides are "reciprocal".
#
# Outputs (in <outdir>)
#   synteny_blocks.tsv        one row per block, with lncRNA counts
#   anchor_pairs.tsv          every collinear gene pair, with coordinates
#   conserved_lncRNA_pairs.tsv
#   lncRNA_status_<sp>.tsv    per lncRNA: conserved / in_block_no_partner / outside_blocks
#   synteny_lncRNA.svg / .pdf / .png
#
# Usage
#   Rscript 21_lncRNA_synteny.R --base . --sp1 A.tha --code1 at --sp2 O.sat --code2 os
#   Rscript 21_lncRNA_synteny.R --help
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(RIdeogram)
  library(optparse)
})

option_list <- list(
  make_option("--base", type = "character", default = ".",
              help = "project root holding <species>/ and synteny/ [default: %default]"),
  make_option("--sp1", type = "character", default = "A.tha",
              help = "species 1 (folder name under --base) [default: %default]"),
  make_option("--sp2", type = "character", default = NULL,
              help = "species 2 (folder name under --base)"),
  make_option("--code1", type = "character", default = "at",
              help = "2-letter MCScanX code of species 1 [default: %default]"),
  make_option("--code2", type = "character", default = NULL,
              help = "2-letter MCScanX code of species 2"),
  make_option("--lncrna1", type = "character", default = NULL,
              help = "lncRNA table, species 1 [default: <base>/<sp1>/03_outputs/11_int_class/<sp1>_final_annotation.tsv]"),
  make_option("--lncrna2", type = "character", default = NULL,
              help = "lncRNA table, species 2 [same convention]"),
  make_option("--seq-hits", dest = "seq_hits", type = "character", default = NULL,
              help = "optional TSV of sequence-similar lncRNA gene pairs (sp1 gene, sp2 gene)"),
  make_option("--max-interval", dest = "max_interval", type = "integer", default = 1000000L,
              help = "ignore a target interval wider than this many bp [default: %default]"),
  make_option("--min-width", dest = "min_width", type = "integer", default = 150000L,
              help = "minimum drawn width (bp) of an lncRNA ribbon so it is visible [default: %default]"),
  make_option("--chr-regex", dest = "chr_regex", type = "character", default = "^Chr[0-9]+$",
              help = paste("FIGURE ONLY: draw sequences whose name (as ChrN) matches this regex; unplaced",
                           "scaffolds still take part in blocks and lncRNA calls [default: %default]")),
  make_option("--n-chr", dest = "n_chr", type = "integer", default = NULL,
              help = "FIGURE ONLY, alternative to --chr-regex: draw the N longest sequences of each species"),
  make_option("--outdir", type = "character", default = NULL,
              help = "output directory [default: <base>/synteny/<sp1>_vs_<sp2>/lncRNA_synteny]")
)
opt <- parse_args(
  OptionParser(option_list = option_list,
               description = "Positionally conserved lncRNAs on MCScanX synteny blocks."),
  args = if (interactive()) character(0) else commandArgs(trailingOnly = TRUE)
)

if (is.null(opt$sp2) || is.null(opt$code2)) stop("--sp2 and --code2 are required (see --help)")
sp1 <- opt$sp1; sp2 <- opt$sp2; code1 <- opt$code1; code2 <- opt$code2
if (nchar(code1) != 2 || nchar(code2) != 2 || code1 == code2)
  stop("--code1 and --code2 must be two different 2-letter codes")

pair      <- paste0(sp1, "_vs_", sp2)
tag       <- paste0(code1, "_", code2)
syn_dir   <- file.path(opt$base, "synteny", pair)
coll_file <- file.path(syn_dir, paste0(tag, ".collinearity"))
gff_file  <- file.path(syn_dir, paste0(tag, ".gff"))
kar1_file <- file.path(opt$base, sp1, "01_inputs", paste0(sp1, ".karyotype.txt"))
kar2_file <- file.path(opt$base, sp2, "01_inputs", paste0(sp2, ".karyotype.txt"))
default_lnc <- function(sp) file.path(opt$base, sp, "03_outputs", "11_int_class", paste0(sp, "_final_annotation.tsv"))
lnc1_file <- if (is.null(opt$lncrna1)) default_lnc(sp1) else opt$lncrna1
lnc2_file <- if (is.null(opt$lncrna2)) default_lnc(sp2) else opt$lncrna2
out_dir   <- if (is.null(opt$outdir)) file.path(syn_dir, "lncRNA_synteny") else opt$outdir

for (f in c(coll_file, gff_file, kar1_file, kar2_file, lnc1_file, lnc2_file)) {
  if (!file.exists(f)) stop("Input file not found: ", f)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# colours (hex WITHOUT '#', as RIdeogram wants)
col_block    <- "d9d9d9"   # blocks with no conserved lncRNA
col_block_hl <- "9ecae1"   # blocks carrying >= 1 conserved lncRNA pair
col_lnc      <- "d55e00"   # the lncRNA pairs themselves

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------
norm_chr <- function(x) paste0("Chr", sub("^[Cc][Hh][Rr]", "", trimws(x)))
# MCScanX chromosome ("at1") -> "Chr1", comparable with the karyotype / lncRNA table
chr_norm <- function(x, code) norm_chr(sub(paste0("^", code), "", x))

read_karyotype <- function(path) {
  k <- read.delim(path, sep = "\t", stringsAsFactors = FALSE)
  if (!all(c("Chr", "Start", "End") %in% names(k))) stop(path, " needs columns Chr, Start, End")
  k <- k %>% transmute(Chr = norm_chr(Chr), Start = as.numeric(Start), End = as.numeric(End))
  k[match(str_sort(k$Chr, numeric = TRUE), k$Chr), ]
}

read_lnc <- function(path, karyotype) {
  lnc <- read.delim(path, sep = "\t", stringsAsFactors = FALSE)
  need <- c("Chr", "Start", "Stop", "lncRNA_gene", "Category")
  if (!all(need %in% names(lnc))) stop(path, " is missing: ", paste(setdiff(need, names(lnc)), collapse = ", "))
  lnc %>%
    mutate(Chr = norm_chr(Chr), Category = replace_na(Category, "Unclassified")) %>%
    filter(Chr %in% karyotype$Chr) %>%
    group_by(lncRNA_gene) %>%                       # collapse isoforms to one locus per gene
    summarise(Chr = first(Chr), Start = min(Start), End = max(Stop),
              Category = names(sort(table(Category), decreasing = TRUE))[1], .groups = "drop")
}

parse_collinearity <- function(path) {
  x   <- readLines(path)
  hdr <- grep("^## Alignment", x)
  if (!length(hdr)) stop("No '## Alignment' blocks in ", path)
  blocks <- tibble(line = hdr, txt = x[hdr]) %>%
    extract(txt, c("block", "score", "evalue", "n", "chr_a", "chr_b", "orient"),
            "^## Alignment ([0-9]+): score=([0-9.]+) e_value=(\\S+) N=([0-9]+) (\\S+)&(\\S+) (plus|minus)",
            convert = TRUE)
  pl    <- grep("^\\s*[0-9]+-\\s*[0-9]+:", x, perl = TRUE)
  pairs <- tibble(line = pl, txt = x[pl]) %>%
    extract(txt, c("geneX", "geneY", "pair_evalue"),
            "^\\s*[0-9]+-\\s*[0-9]+:\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s*$") %>%
    mutate(block = blocks$block[findInterval(line, blocks$line)]) %>%
    select(block, geneX, geneY, pair_evalue)
  list(blocks = select(blocks, -line), pairs = pairs)
}

# For every lncRNA (query species) inside a block's span, find the flanking anchor
# genes, take their partners in the target species, and list the target-species
# lncRNAs lying in that window.
#   anc: block, gq (query anchor gene), Chr_q, s_q, e_q, Chr_t, s_t, e_t
map_lnc <- function(lnc_q, lnc_t, anc, max_interval) {
  template <- tibble(block = integer(), lnc_q = character(), cat_q = character(),
                     left_q = character(), right_q = character(), win_chr_t = character(),
                     win_lo = numeric(), win_hi = numeric(), interval_ok = logical(),
                     lnc_t = character(), cat_t = character())
  lnc_t <- mutate(lnc_t, Mid = (Start + End) / 2)
  rows  <- list()
  for (b in unique(anc$block)) {
    ab <- anc[anc$block == b, ]
    ab <- ab[order(ab$s_q), ]
    n  <- nrow(ab)
    if (n < 2) next
    lq <- lnc_q[lnc_q$Chr == ab$Chr_q[1] & lnc_q$End >= min(ab$s_q) & lnc_q$Start <= max(ab$e_q), , drop = FALSE]
    if (nrow(lq) == 0) next
    mid    <- (lq$Start + lq$End) / 2
    i      <- pmin(pmax(findInterval(mid, ab$s_q), 1L), n - 1L)   # left flanking anchor
    win_lo <- pmin(ab$s_t[i], ab$s_t[i + 1L])
    win_hi <- pmax(ab$e_t[i], ab$e_t[i + 1L])
    ok     <- (win_hi - win_lo) <= max_interval
    lt_chr <- lnc_t[lnc_t$Chr == ab$Chr_t[1], , drop = FALSE]
    for (k in seq_len(nrow(lq))) {
      hit <- if (ok[k]) lt_chr[lt_chr$Mid >= win_lo[k] & lt_chr$Mid <= win_hi[k], , drop = FALSE] else lt_chr[0, ]
      rows[[length(rows) + 1L]] <- tibble(
        block = b, lnc_q = lq$lncRNA_gene[k], cat_q = lq$Category[k],
        left_q = ab$gq[i[k]], right_q = ab$gq[i[k] + 1L],
        win_chr_t = ab$Chr_t[1], win_lo = win_lo[k], win_hi = win_hi[k], interval_ok = ok[k],
        lnc_t = if (nrow(hit)) hit$lncRNA_gene else NA_character_,
        cat_t = if (nrow(hit)) hit$Category else NA_character_)
    }
  }
  if (length(rows)) bind_rows(rows) else template
}

# ------------------------------------------------------------------------------
# 1. LOAD
# ------------------------------------------------------------------------------
k1 <- read_karyotype(kar1_file); k2 <- read_karyotype(kar2_file)
lnc1 <- read_lnc(lnc1_file, k1); lnc2 <- read_lnc(lnc2_file, k2)
message(sprintf("lncRNA loci: %s = %d, %s = %d", sp1, nrow(lnc1), sp2, nrow(lnc2)))

col   <- parse_collinearity(coll_file)
genes <- read.delim(gff_file, header = FALSE, quote = "", stringsAsFactors = FALSE,
                    col.names = c("chr", "gene", "start", "end"),
                    colClasses = c("character", "character", "numeric", "numeric")) %>% as_tibble()

# ------------------------------------------------------------------------------
# 2. ANCHOR PAIRS AND BLOCKS (species 1 = left, species 2 = right)
# ------------------------------------------------------------------------------
p  <- col$pairs
sx <- substr(genes$chr[match(p$geneX, genes$gene)], 1, 2)
if (anyNA(sx)) warning(sum(is.na(sx)), " collinear genes are not in ", basename(gff_file))
p <- p %>% mutate(g1 = if_else(sx == code1, geneX, geneY),
                  g2 = if_else(sx == code1, geneY, geneX))
i1 <- match(p$g1, genes$gene); i2 <- match(p$g2, genes$gene)
anchors <- p %>%
  mutate(chr1 = genes$chr[i1], s1 = genes$start[i1], e1 = genes$end[i1],
         chr2 = genes$chr[i2], s2 = genes$start[i2], e2 = genes$end[i2]) %>%
  filter(substr(chr1, 1, 2) == code1, substr(chr2, 1, 2) == code2) %>%   # drops intra-species blocks
  mutate(Chr1 = chr_norm(chr1, code1), Chr2 = chr_norm(chr2, code2)) %>%
  select(block, g1, g2, Chr1, s1, e1, Chr2, s2, e2, pair_evalue)
if (nrow(anchors) == 0) stop("No inter-species anchor pairs; check --code1/--code2 against the gff")

blocks <- anchors %>%
  group_by(block) %>%
  summarise(Chr1 = first(Chr1), Start1 = min(s1), End1 = max(e1),
            Chr2 = first(Chr2), Start2 = min(s2), End2 = max(e2),
            n_genes = n(), .groups = "drop") %>%
  left_join(select(col$blocks, block, score, orient), by = "block")
message(sprintf("MCScanX: %d blocks, %d anchor pairs", nrow(blocks), nrow(anchors)))

# ------------------------------------------------------------------------------
# 3. POSITIONALLY CONSERVED lncRNAs (both directions)
# ------------------------------------------------------------------------------
anc12 <- anchors %>% transmute(block, gq = g1, Chr_q = Chr1, s_q = s1, e_q = e1, Chr_t = Chr2, s_t = s2, e_t = e2)
anc21 <- anchors %>% transmute(block, gq = g2, Chr_q = Chr2, s_q = s2, e_q = e2, Chr_t = Chr1, s_t = s1, e_t = e1)
d12 <- map_lnc(lnc1, lnc2, anc12, opt$max_interval)
d21 <- map_lnc(lnc2, lnc1, anc21, opt$max_interval)

pairs_long <- bind_rows(
  d12 %>% filter(!is.na(lnc_t)) %>% transmute(lnc1 = lnc_q, lnc2 = lnc_t, block, dir = "1to2", cat1 = cat_q, cat2 = cat_t),
  d21 %>% filter(!is.na(lnc_t)) %>% transmute(lnc1 = lnc_t, lnc2 = lnc_q, block, dir = "2to1", cat1 = cat_t, cat2 = cat_q))

pairs <- pairs_long %>%
  group_by(lnc1, lnc2) %>%
  summarise(blocks = paste(sort(unique(block)), collapse = ","), cat1 = first(cat1), cat2 = first(cat2),
            reciprocal = n_distinct(dir) == 2, .groups = "drop") %>%
  mutate(same_category = cat1 == cat2)

if (!is.null(opt$seq_hits)) {
  sh <- read_tsv(opt$seq_hits, col_names = c("lnc1", "lnc2"), col_types = "cc", show_col_types = FALSE)
  pairs$seq_support <- paste(pairs$lnc1, pairs$lnc2) %in% paste(sh$lnc1, sh$lnc2)
} else {
  pairs$seq_support <- NA
}
pairs <- pairs %>%
  mutate(support = paste0(if_else(reciprocal, "positional_reciprocal", "positional_one_way"),
                          if_else(!is.na(seq_support) & seq_support, "+sequence", ""))) %>%
  arrange(desc(reciprocal), desc(seq_support), lnc1)
message(sprintf("Conserved lncRNA pairs: %d (%d reciprocal)", nrow(pairs), sum(pairs$reciprocal)))

status_tbl <- function(lnc, d, pairs, me, other) {
  partners <- pairs %>% group_by(.data[[me]]) %>%
    summarise(partners = paste(.data[[other]], collapse = ","), .groups = "drop") %>%
    rename(lncRNA_gene = 1)
  lnc %>%
    left_join(partners, by = "lncRNA_gene") %>%
    mutate(in_block = lncRNA_gene %in% d$lnc_q,
           status = case_when(!is.na(partners) ~ "conserved",
                              in_block ~ "in_block_no_partner",
                              TRUE ~ "outside_blocks")) %>%
    select(-in_block)
}
st1 <- status_tbl(lnc1, d12, pairs, "lnc1", "lnc2")
st2 <- status_tbl(lnc2, d21, pairs, "lnc2", "lnc1")
for (s in list(list(sp1, st1), list(sp2, st2))) {
  message(sprintf("%s: %s", s[[1]], paste(names(table(s[[2]]$status)), table(s[[2]]$status), sep = " = ", collapse = ", ")))
}

blk <- blocks %>%
  left_join(d12 %>% distinct(block, lnc_q) %>% count(block, name = "n_lnc1"), by = "block") %>%
  left_join(d21 %>% distinct(block, lnc_q) %>% count(block, name = "n_lnc2"), by = "block") %>%
  left_join(pairs_long %>% distinct(block, lnc1, lnc2) %>% count(block, name = "n_conserved_pairs"), by = "block") %>%
  mutate(across(c(n_lnc1, n_lnc2, n_conserved_pairs), ~ replace_na(.x, 0L)))

write_tsv(blk,     file.path(out_dir, "synteny_blocks.tsv"))
write_tsv(anchors, file.path(out_dir, "anchor_pairs.tsv"))
write_tsv(pairs,   file.path(out_dir, "conserved_lncRNA_pairs.tsv"))
write_tsv(st1,     file.path(out_dir, paste0("lncRNA_status_", sp1, ".tsv")))
write_tsv(st2,     file.path(out_dir, paste0("lncRNA_status_", sp2, ".tsv")))

# ------------------------------------------------------------------------------
# 4. RIDEOGRAM SYNTENY FIGURE
# ------------------------------------------------------------------------------
# Everything above used ALL nuclear sequences (scaffolds included). Only the figure is limited
# to chromosomes; blocks / lncRNA pairs touching an unplaced scaffold stay in the tables.
pick_chr <- function(k) {
  kp <- if (!is.null(opt$n_chr)) slice_max(k, End, n = opt$n_chr, with_ties = FALSE) else filter(k, grepl(opt$chr_regex, Chr))
  if (nrow(kp) == 0) stop("No sequence matches --chr-regex '", opt$chr_regex, "'. Names present (first 10): ",
                          paste(head(k$Chr, 10), collapse = ", "), ". Adjust --chr-regex or use --n-chr.")
  kp[match(str_sort(kp$Chr, numeric = TRUE), kp$Chr), ]
}
k1p <- pick_chr(k1); k2p <- pick_chr(k2)
message(sprintf("Figure: %d/%d sequences of %s, %d/%d of %s", nrow(k1p), nrow(k1), sp1, nrow(k2p), nrow(k2), sp2))

karyo <- bind_rows(
  k1p %>% mutate(Chr = sub("^Chr", "", Chr), fill = "969696", species = sp1, size = 12, color = "ffffff"),
  k2p %>% mutate(Chr = sub("^Chr", "", Chr), fill = "969696", species = sp2, size = 12, color = "ffffff"))
karyo <- as.data.frame(karyo)

# Species_1 / Species_2 in the synteny table = chromosome index within that species'
# block of the karyotype (1, 2, 3 ...). CHECK against data(synteny_dual_comparison).
idx1 <- setNames(seq_len(nrow(k1p)), k1p$Chr)
idx2 <- setNames(seq_len(nrow(k2p)), k2p$Chr)

ribbon <- function(df, fill) {
  tibble(Species_1 = as.integer(unname(idx1[df$Chr1])), Start_1 = round(df$Start1), End_1 = round(df$End1),
         Species_2 = as.integer(unname(idx2[df$Chr2])), Start_2 = round(df$Start2), End_2 = round(df$End2),
         fill = rep(fill, nrow(df)))
}

# lncRNA pairs: widen to --min-width so a few-kb locus is visible on a ~100 mm chromosome
w    <- opt$min_width
len1 <- setNames(k1$End, k1$Chr); len2 <- setNames(k2$End, k2$Chr)
lp <- pairs %>%
  left_join(lnc1 %>% select(lnc1 = lncRNA_gene, Chr1 = Chr, S1 = Start, E1 = End), by = "lnc1") %>%
  left_join(lnc2 %>% select(lnc2 = lncRNA_gene, Chr2 = Chr, S2 = Start, E2 = End), by = "lnc2") %>%
  mutate(M1 = (S1 + E1) / 2, M2 = (S2 + E2) / 2,
         Start1 = pmax(M1 - w / 2, 1), End1 = pmin(M1 + w / 2, unname(len1[Chr1])),
         Start2 = pmax(M2 - w / 2, 1), End2 = pmin(M2 + w / 2, unname(len2[Chr2])))

syn <- bind_rows(
  ribbon(filter(blk, n_conserved_pairs == 0), col_block),      # drawn first (underneath)
  ribbon(filter(blk, n_conserved_pairs > 0),  col_block_hl),
  ribbon(lp, col_lnc)) %>%                                     # drawn last (on top)
  as.data.frame()
n_drawn <- sum(!is.na(syn$Species_1) & !is.na(syn$Species_2))
message(sprintf("Ribbons drawn: %d of %d (the rest touch sequences not drawn)", n_drawn, nrow(syn)))
syn <- syn[!is.na(syn$Species_1) & !is.na(syn$Species_2), ]

svg_file <- file.path(out_dir, "synteny_lncRNA.svg")
ideogram(karyotype = karyo, synteny = syn, output = svg_file)
tryCatch({
  rsvg::rsvg_pdf(svg_file, sub("\\.svg$", ".pdf", svg_file))
  rsvg::rsvg_png(svg_file, sub("\\.svg$", ".png", svg_file), width = 2480)
}, error = function(e) message("SVG written; PDF/PNG conversion skipped: ", conditionMessage(e)))

message("Done. Outputs in ", out_dir)
