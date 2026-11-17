#!/usr/bin/env Rscript
# 21_lncRNA_synteny.R
#
# Conserved lncRNA loci from MCScanX syntenic blocks (base R only, no packages needed).
#
# Idea: lncRNAs are not anchors in the MCScanX blocks (only mRNAs are). For every lncRNA we find the
# two mRNA genes that flank it, check whether it lies between two anchor genes of a syntenic block,
# follow those anchors to the other species, and ask which lncRNAs sit between the matching anchors
# there. Done in both directions; a pair found in BOTH directions is called "reciprocal".
#
# Inputs (all written by earlier steps):
#   $BASE/synteny/species/<SPEC>/<SPEC>.gff                       <code><chr> gene start end  (synteny prep)
#   $BASE/synteny/<SP1>_vs_<SP2>/<CODE1>_<CODE2>.collinearity     MCScanX output
#   $BASE/<SPEC>/03_outputs/12_finalize/<SPEC>_lncRNA_collapsed.gtf      (--isoforms collapsed, default)
#   $BASE/<SPEC>/03_outputs/12_finalize/<SPEC>_lncRNA_all_isoforms.gtf   (--isoforms all)
#   Isoforms are always collapsed to one locus per gene_id (min start, max end).
#
# Usage
#   One pair:
#     Rscript 21_lncRNA_synteny.R --base $BASE --sp1 A.tha --code1 at --sp2 G.max --code2 gm
#   All pairs of a species list + multi-species groups (default list = the 4-species run):
#     Rscript 21_lncRNA_synteny.R --base $BASE --species "A.tha:at,O.sat:os,G.max:gm,Z.mays:zm"
#
# Options
#   --isoforms collapsed|all   which lncRNA GTF to read (default collapsed)
#   --max-span N               ignore blocks where the two bracketing anchors are more than N genes apart
#                              in either species (default 20); keeps "conserved" meaningful
#   --include-nonreciprocal    all-species mode: also use pairs found in only one direction for grouping
#
# Outputs
#   pair mode:  $BASE/synteny/<SP1>_vs_<SP2>/lncRNA_synteny/
#                 <SP1>_vs_<SP2>.lncRNA_pairs.tsv    lncRNA pairs (one row per lnc1/lnc2 pair)
#                 <SP1>_vs_<SP2>.lncRNA_status.tsv   every lncRNA of both species and whether it sits in a block / has partners
#   all mode additionally: $BASE/synteny/lncRNA_conservation/
#                 lncRNA_conserved_groups.tsv        groups of lncRNAs linked across species
#                 lncRNA_group_members.tsv           one row per lncRNA that is in a group

options(stringsAsFactors = FALSE)

## ------------------------------------------------------------------ args
parse_args <- function(a) {
  out <- list(); i <- 1
  while (i <= length(a)) {
    if (grepl("^--", a[i])) {
      key <- sub("^--", "", a[i])
      if (i < length(a) && !grepl("^--", a[i + 1])) { out[[key]] <- a[i + 1]; i <- i + 2 }
      else { out[[key]] <- TRUE; i <- i + 1 }
    } else i <- i + 1
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
getopt <- function(k, d = NULL) if (!is.null(opt[[k]])) opt[[k]] else d

BASE     <- getopt("base");  if (is.null(BASE)) stop("--base is required")
ISO      <- getopt("isoforms", "collapsed")
if (!ISO %in% c("collapsed", "all")) stop("--isoforms must be 'collapsed' or 'all'")
MAX_SPAN <- as.integer(getopt("max-span", 20))
INCL_NR  <- isTRUE(opt[["include-nonreciprocal"]])

## ------------------------------------------------------------------ readers
read_lnc_gtf <- function(path, code) {
  x <- readLines(path, warn = FALSE)
  x <- x[nzchar(x) & !startsWith(x, "#")]
  d <- read.delim(text = x, header = FALSE, quote = "", fill = TRUE, stringsAsFactors = FALSE,
                  col.names = c("seqname", "source", "feature", "start", "end",
                                "score", "strand", "frame", "attr"))
  d <- d[d$feature %in% c("exon", "transcript"), , drop = FALSE]
  gid <- ifelse(grepl('gene_id "', d$attr, fixed = TRUE),
                sub('.*gene_id "([^"]+)".*', "\\1", d$attr), NA_character_)
  tid <- ifelse(grepl('transcript_id "', d$attr, fixed = TRUE),
                sub('.*transcript_id "([^"]+)".*', "\\1", d$attr), NA_character_)
  gid[is.na(gid)] <- tid[is.na(gid)]
  ok <- !is.na(gid)
  d <- d[ok, , drop = FALSE]; gid <- gid[ok]
  if (!nrow(d)) stop("no exon/transcript lines with gene_id in ", path)
  # same chromosome labelling as the MCScanX gff: <code> + name without leading chr/Chr
  chr <- paste0(code, sub("^[Cc][Hh][Rr]", "", d$seqname))
  s <- tapply(as.numeric(d$start), gid, min)
  e <- tapply(as.numeric(d$end),   gid, max)
  first <- !duplicated(gid)
  info <- data.frame(lnc = gid[first], chr = chr[first], strand = d$strand[first])
  info$start <- as.integer(s[info$lnc])
  info$end   <- as.integer(e[info$lnc])
  info
}

# flanking mRNA genes of every lncRNA (indices in the chromosome-ordered mRNA gene list)
add_flanks <- function(lnc, genes) {
  lnc$L <- NA_integer_; lnc$R <- NA_integer_; lnc$overlaps_mRNA <- NA
  gs <- split(genes, genes$chr)
  for (ch in intersect(unique(lnc$chr), names(gs))) {
    g <- gs[[ch]]; k <- which(lnc$chr == ch)
    L <- findInterval(lnc$start[k], g$start)        # last gene starting at/before the lncRNA start
    R <- findInterval(lnc$end[k],   g$start) + 1L   # first gene starting after the lncRNA end
    ov <- (R - L - 1L) > 0L | (L >= 1L & g$end[pmax(L, 1L)] >= lnc$start[k])
    L[L < 1L] <- NA; R[R > nrow(g)] <- NA
    lnc$L[k] <- L; lnc$R[k] <- R; lnc$overlaps_mRNA[k] <- ov
  }
  lnc
}

load_species <- function(base, spec, code, iso) {
  gff <- file.path(base, "synteny", "species", spec, paste0(spec, ".gff"))
  gtf <- file.path(base, spec, "03_outputs", "12_finalize",
                   paste0(spec, "_lncRNA_", if (iso == "all") "all_isoforms" else "collapsed", ".gtf"))
  for (f in c(gff, gtf)) if (!file.exists(f)) stop("missing file: ", f)
  genes <- read.table(gff, sep = "\t", header = FALSE, quote = "", comment.char = "",
                      col.names = c("chr", "gene", "start", "end"), stringsAsFactors = FALSE)
  genes <- genes[order(genes$chr, genes$start, genes$end), ]
  genes$idx <- ave(seq_len(nrow(genes)), genes$chr, FUN = seq_along)
  rownames(genes) <- NULL
  lnc <- read_lnc_gtf(gtf, code)
  shared <- intersect(unique(lnc$chr), unique(genes$chr))
  if (!length(shared))
    stop(spec, ": no lncRNA chromosome label matches the mRNA gff. lncRNA labels e.g. ",
         paste(head(unique(lnc$chr), 3), collapse = ","), " vs gff labels e.g. ",
         paste(head(unique(genes$chr), 3), collapse = ","))
  n_off <- sum(!lnc$chr %in% shared)
  message(sprintf("%s: %d mRNA genes, %d lncRNA loci (%s)%s", spec, nrow(genes), nrow(lnc), iso,
                  if (n_off) sprintf("; %d lncRNAs on sequences with no mRNA genes (ignored)", n_off) else ""))
  list(name = spec, code = code, genes = genes, lnc = add_flanks(lnc, genes))
}

read_collinearity <- function(path) {
  x <- readLines(path, warn = FALSE)
  h <- x[startsWith(x, "## Alignment")]
  hm <- regmatches(h, regexec(
    "^## Alignment ([0-9]+): score=([^[:space:]]+) e_value=([^[:space:]]+) N=([0-9]+) ([^&[:space:]]+)&([^[:space:]]+) (plus|minus)", h))
  if (!length(h) || !all(lengths(hm) == 8)) stop("could not parse block headers in ", path)
  blocks <- data.frame(block = as.integer(vapply(hm, `[`, "", 2)),
                       orient = vapply(hm, `[`, "", 8))
  a <- x[grepl("^[[:space:]]*[0-9]+-[[:space:]]*[0-9]+:", x)]
  am <- regmatches(a, regexec(
    "^[[:space:]]*([0-9]+)-[[:space:]]*([0-9]+):[[:space:]]+([^[:space:]]+)[[:space:]]+([^[:space:]]+)", a))
  anchors <- data.frame(block = as.integer(vapply(am, `[`, "", 2)),
                        g1 = vapply(am, `[`, "", 4), g2 = vapply(am, `[`, "", 5))
  list(blocks = blocks, anchors = anchors)
}

## ------------------------------------------------------------------ core: follow blocks from one species to the other
empty_inb <- data.frame(lnc = character(), block = integer(), anchor_left = character(),
                        anchor_right = character(), span_F = integer(), span_T = integer(),
                        too_wide = logical(), n_partners = integer())
empty_par <- data.frame(lncF = character(), lncT = character(), block = integer(),
                        span_F = numeric(), span_T = numeric(), strand_concordant = logical())

map_direction <- function(lncF, lncT, an, ori, max_span) {
  flip <- c("+" = "-", "-" = "+")
  lncF <- lncF[!is.na(lncF$L) & !is.na(lncF$R), , drop = FALSE]
  lncT <- lncT[!is.na(lncT$L) & !is.na(lncT$R), , drop = FALSE]
  fby <- split(seq_len(nrow(lncF)), lncF$chr)
  tby <- split(seq_len(nrow(lncT)), lncT$chr)
  inb <- list(); par <- list()
  ab <- split(an, an$block)
  for (b in names(ab)) {
    a <- ab[[b]]; a <- a[order(a$idxF), ]
    kF <- fby[[a$chrF[1]]]
    if (is.null(kF)) next
    sub <- lncF[kF, , drop = FALSE]
    # lncRNA must lie between two anchors of this block: lower anchor <= left flank, upper anchor >= right flank
    sel <- sub$L >= min(a$idxF) & sub$R <= max(a$idxF)
    if (!any(sel)) next
    sub <- sub[sel, , drop = FALSE]
    lo <- findInterval(sub$L, a$idxF)
    up <- findInterval(sub$R - 1L, a$idxF) + 1L
    tl <- a$idxT[lo]; tu <- a$idxT[up]
    lo2 <- pmin(tl, tu); hi2 <- pmax(tl, tu)
    spF <- a$idxF[up] - a$idxF[lo]; spT <- hi2 - lo2
    wide <- spF > max_span | spT > max_span
    kT <- tby[[a$chrT[1]]]
    cand <- if (is.null(kT)) lncT[0, , drop = FALSE] else lncT[kT, , drop = FALSE]
    npart <- integer(nrow(sub))
    for (i in seq_len(nrow(sub))) {
      if (wide[i] || nrow(cand) == 0) next
      m <- cand$L >= lo2[i] & cand$R <= hi2[i]     # partner lies between the matching anchors
      if (!any(m)) next
      npart[i] <- sum(m)
      pc <- cand[m, , drop = FALSE]
      expS <- if (ori[[b]] == "minus") unname(flip[sub$strand[i]]) else sub$strand[i]
      par[[length(par) + 1]] <- data.frame(
        lncF = sub$lnc[i], lncT = pc$lnc, block = as.integer(b), span_F = spF[i], span_T = spT[i],
        strand_concordant = if (!is.na(expS) && expS %in% c("+", "-")) pc$strand == expS else NA)
    }
    inb[[length(inb) + 1]] <- data.frame(
      lnc = sub$lnc, block = as.integer(b), anchor_left = a$gF[lo], anchor_right = a$gF[up],
      span_F = spF, span_T = spT, too_wide = wide, n_partners = npart)
  }
  list(inblock  = if (length(inb)) do.call(rbind, inb) else empty_inb,
       partners = if (length(par)) do.call(rbind, par) else empty_par)
}

as_pairs <- function(par, swap) {
  if (!nrow(par))
    return(data.frame(lnc1 = character(), lnc2 = character(), block = integer(),
                      span_1 = numeric(), span_2 = numeric(), strand_concordant = logical()))
  data.frame(lnc1 = if (swap) par$lncT else par$lncF,
             lnc2 = if (swap) par$lncF else par$lncT,
             block = par$block,
             span_1 = if (swap) par$span_T else par$span_F,
             span_2 = if (swap) par$span_F else par$span_T,
             strand_concordant = par$strand_concordant)
}
best_pair <- function(p) {   # one row per lnc pair: the tightest block
  if (!nrow(p)) return(p)
  p <- p[order(p$span_1 + p$span_2), ]
  p[!duplicated(paste(p$lnc1, p$lnc2)), ]
}

run_pair <- function(base, A, B, max_span) {
  pair <- paste0(A$name, "_vs_", B$name); tag <- paste0(A$code, "_", B$code)
  cf <- file.path(base, "synteny", pair, paste0(tag, ".collinearity"))
  if (!file.exists(cf)) {
    alt <- file.path(base, "synteny", paste0(B$name, "_vs_", A$name),
                     paste0(B$code, "_", A$code, ".collinearity"))
    if (file.exists(alt)) cf <- alt else stop("collinearity file not found: ", cf)
  }
  outdir <- file.path(base, "synteny", pair, "lncRNA_synteny")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("[", format(Sys.time()), "] ", pair)

  cl <- read_collinearity(cf)
  an <- cl$anchors
  if (!nrow(an)) stop("no anchor lines in ", cf)
  # the first gene column belongs to whichever species MCScanX listed first; orient to A/B by gene id
  if (mean(an$g1 %in% A$genes$gene) < mean(an$g1 %in% B$genes$gene)) {
    tmp <- an$g1; an$g1 <- an$g2; an$g2 <- tmp
  }
  ia <- match(an$g1, A$genes$gene); ib <- match(an$g2, B$genes$gene)
  ok <- !is.na(ia) & !is.na(ib)
  if (any(!ok)) message(sprintf("  %d of %d anchor pairs have a gene id missing from the gff (dropped)", sum(!ok), nrow(an)))
  if (!any(ok)) stop("no anchor gene ids match the species gffs; check species/code order")
  an <- data.frame(block = an$block[ok],
                   gA = an$g1[ok], chrA = A$genes$chr[ia[ok]], idxA = A$genes$idx[ia[ok]],
                   gB = an$g2[ok], chrB = B$genes$chr[ib[ok]], idxB = B$genes$idx[ib[ok]])
  ori <- setNames(cl$blocks$orient, cl$blocks$block)

  anAB <- data.frame(block = an$block, gF = an$gA, chrF = an$chrA, idxF = an$idxA,
                     gT = an$gB, chrT = an$chrB, idxT = an$idxB)
  anBA <- data.frame(block = an$block, gF = an$gB, chrF = an$chrB, idxF = an$idxB,
                     gT = an$gA, chrT = an$chrA, idxT = an$idxA)
  r12 <- map_direction(A$lnc, B$lnc, anAB, ori, max_span)
  r21 <- map_direction(B$lnc, A$lnc, anBA, ori, max_span)

  p12 <- best_pair(as_pairs(r12$partners, FALSE))
  p21 <- best_pair(as_pairs(r21$partners, TRUE))
  names(p12)[-(1:2)] <- paste0(names(p12)[-(1:2)], ".1to2")
  names(p21)[-(1:2)] <- paste0(names(p21)[-(1:2)], ".2to1")
  pairs <- merge(p12, p21, by = c("lnc1", "lnc2"), all = TRUE)
  pairs$in_1to2    <- !is.na(pairs$block.1to2)
  pairs$in_2to1    <- !is.na(pairs$block.2to1)
  pairs$reciprocal <- pairs$in_1to2 & pairs$in_2to1
  pairs <- pairs[order(!pairs$reciprocal, pairs$lnc1, pairs$lnc2), , drop = FALSE]
  pairs <- cbind(sp1 = rep(A$name, nrow(pairs)), sp2 = rep(B$name, nrow(pairs)), pairs)
  rownames(pairs) <- NULL

  status <- function(S, inb, col) {
    d <- S$lnc[, c("lnc", "chr", "start", "end", "strand", "L", "R", "overlaps_mRNA")]
    names(d)[names(d) == "L"] <- "flank_left_idx"; names(d)[names(d) == "R"] <- "flank_right_idx"
    d$in_block        <- d$lnc %in% inb$lnc
    d$in_narrow_block <- d$lnc %in% inb$lnc[!inb$too_wide]
    d$n_partners      <- as.integer(table(factor(pairs[[col]], levels = d$lnc)))
    d$n_reciprocal    <- as.integer(table(factor(pairs[[col]][pairs$reciprocal], levels = d$lnc)))
    cbind(species = rep(S$name, nrow(d)), d)
  }
  st <- rbind(status(A, r12$inblock, "lnc1"), status(B, r21$inblock, "lnc2"))

  write.table(pairs, file.path(outdir, paste0(pair, ".lncRNA_pairs.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(st, file.path(outdir, paste0(pair, ".lncRNA_status.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  sa <- st[st$species == A$name, ]; sb <- st[st$species == B$name, ]
  message(sprintf("  %s: %d lncRNAs, %d in blocks, %d with a partner, %d with a reciprocal partner",
                  A$name, nrow(sa), sum(sa$in_block), sum(sa$n_partners > 0), sum(sa$n_reciprocal > 0)))
  message(sprintf("  %s: %d lncRNAs, %d in blocks, %d with a partner, %d with a reciprocal partner",
                  B$name, nrow(sb), sum(sb$in_block), sum(sb$n_partners > 0), sum(sb$n_reciprocal > 0)))
  message(sprintf("  lncRNA pairs: %d (%d reciprocal)  -> %s", nrow(pairs), sum(pairs$reciprocal), outdir))
  pairs
}

## ------------------------------------------------------------------ multi-species groups
connected_components <- function(a, b) {
  nodes <- unique(c(a, b))
  parent <- seq_along(nodes)
  find <- function(x) {
    while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }
    x
  }
  ia <- match(a, nodes); ib <- match(b, nodes)
  for (k in seq_along(ia)) {
    ra <- find(ia[k]); rb <- find(ib[k])
    if (ra != rb) parent[rb] <- ra
  }
  data.frame(node = nodes, comp = vapply(seq_along(nodes), find, integer(1)))
}

combine_groups <- function(base, species_names, pair_list, include_nonrecip) {
  edges <- do.call(rbind, lapply(pair_list, function(p) {
    if (!include_nonrecip) p <- p[p$reciprocal, , drop = FALSE]
    if (!nrow(p)) return(NULL)
    data.frame(a = paste(p$sp1, p$lnc1, sep = "|"), b = paste(p$sp2, p$lnc2, sep = "|"))
  }))
  if (is.null(edges) || !nrow(edges)) { message("no lncRNA pairs to group"); return(invisible(NULL)) }
  nodes <- connected_components(edges$a, edges$b)
  nodes$species <- sub("\\|.*", "", nodes$node)
  nodes$lnc     <- sub("^[^|]*\\|", "", nodes$node)
  comps <- unique(nodes$comp)
  groups <- data.frame(comp = comps,
                       n_species = vapply(comps, function(cc) length(unique(nodes$species[nodes$comp == cc])), 1L),
                       n_lnc     = vapply(comps, function(cc) sum(nodes$comp == cc), 1L))
  groups <- groups[order(-groups$n_species, -groups$n_lnc), ]
  groups$group <- sprintf("LNCG%05d", seq_len(nrow(groups)))
  groups$species <- vapply(groups$comp, function(cc)
    paste(species_names[species_names %in% nodes$species[nodes$comp == cc]], collapse = ","), "")
  for (s in species_names)
    groups[[s]] <- vapply(groups$comp, function(cc)
      paste(sort(nodes$lnc[nodes$comp == cc & nodes$species == s]), collapse = ";"), "")
  nodes$group <- groups$group[match(nodes$comp, groups$comp)]
  nodes <- nodes[order(nodes$group, nodes$species), c("group", "species", "lnc")]

  outdir <- file.path(base, "synteny", "lncRNA_conservation")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  groups <- groups[, c("group", "n_species", "species", "n_lnc", species_names)]
  write.table(groups, file.path(outdir, "lncRNA_conserved_groups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(nodes,  file.path(outdir, "lncRNA_group_members.tsv"),   sep = "\t", quote = FALSE, row.names = FALSE)
  message("lncRNA groups by number of species:")
  print(table(groups$n_species))
  message("-> ", outdir)
  invisible(groups)
}

## ------------------------------------------------------------------ main
if (!is.null(opt[["sp1"]])) {
  # single pair
  for (k in c("sp2", "code1", "code2")) if (is.null(opt[[k]])) stop("--", k, " is required in pair mode")
  A <- load_species(BASE, opt[["sp1"]], opt[["code1"]], ISO)
  B <- load_species(BASE, opt[["sp2"]], opt[["code2"]], ISO)
  invisible(run_pair(BASE, A, B, MAX_SPAN))
} else {
  # all unordered pairs of the species list, then multi-species groups
  sp_arg <- getopt("species", "A.tha:at,O.sat:os,G.max:gm,Z.mays:zm")
  ent <- strsplit(strsplit(sp_arg, ",")[[1]], ":")
  SPN <- vapply(ent, `[`, "", 1); CODES <- vapply(ent, `[`, "", 2)
  if (anyNA(CODES) || anyDuplicated(CODES)) stop("--species needs unique name:code entries")
  SP <- setNames(Map(function(s, c) load_species(BASE, s, c, ISO), SPN, CODES), SPN)
  pair_list <- list()
  for (i in seq_along(SPN)) for (j in seq_along(SPN)) if (i < j)
    pair_list[[length(pair_list) + 1]] <- run_pair(BASE, SP[[i]], SP[[j]], MAX_SPAN)
  combine_groups(BASE, SPN, pair_list, INCL_NR)
}
