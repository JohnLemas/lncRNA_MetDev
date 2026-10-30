# Install required CRAN packages if not already installed:
# install.packages(c("tidyverse", "pheatmap", "ggridges", "extrafont", "ggVennDiagram", "UpSetR"))

library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggridges)
library(extrafont)
library(ggVennDiagram)
library(UpSetR)

# Load fonts and register Times New Roman for PDF device
# Note: Run font_import(pattern = "Times") once if Times New Roman isn't registered on your system.
loadfonts(device = "pdf", quiet = TRUE)
font_family <- "Times New Roman"

# ------------------------------------------------------------------------------
# 0. EXPRESSION THRESHOLDS (edit here; applied consistently to every figure)
# ------------------------------------------------------------------------------
# Lower: permissive, because these are ab initio predictions (not ground truth).
# Upper: lncRNAs are typically lowly expressed. Transcripts whose max tissue mean
#        TPM exceeds this are treated as likely mis-classified (e.g. mRNA-like
#        or unannotated protein-coding) and removed. In this dataset 50 TPM
#        sits above the ~97th percentile of expressed transcripts and removes
#        only a handful of extreme outliers (e.g. MSTRG.2072.1 at ~600-800 TPM).
LOWER_TPM <- 0.1
UPPER_TPM <- 50

# ------------------------------------------------------------------------------
# 1. LOAD & PREPARE DATA
# ------------------------------------------------------------------------------

setwd("~/Library/CloudStorage/OneDrive-MichiganStateUniversity/Dissertation/chapter4/lncRNA_met_dev")
# Output directory for all figures/tables (created if it doesn't exist)
out_dir <- "species/03_outputs/13_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

annotation_data <- read.delim("species/03_outputs/11_int_class/species_final_annotation.tsv", sep = "\t", stringsAsFactors = FALSE)

# These columns need to come from colnames() of the annotation table
# below are examples from Arabidopsis thaliana
tissue_cols <- c(
  seedling       = "mean_TPM_seedling",
  Inflorescences = "mean_TPM_Inflorescences",
  leaf           = "mean_TPM_leaf",
  silique        = "mean_TPM_silique",
  rosette        = "mean_TPM_rosette",
  root           = "mean_TPM_root"
)

# Max TPM across tissues for every transcript, then apply BOTH thresholds
annotation_data$max_tissue_TPM <- apply(annotation_data[, tissue_cols], 1, max)

keep <- annotation_data$max_tissue_TPM > LOWER_TPM &
  annotation_data$max_tissue_TPM <= UPPER_TPM

annotation_filt <- annotation_data[keep, ]
rownames(annotation_filt) <- NULL  # drop carried-over row names (column_to_rownames() requires none)

message(sprintf(
  "Transcripts: %d total | %d with max TPM > %g | %d removed by upper cutoff (> %g) | %d retained",
  nrow(annotation_data),
  sum(annotation_data$max_tissue_TPM > LOWER_TPM),
  LOWER_TPM,
  sum(annotation_data$max_tissue_TPM > UPPER_TPM),
  UPPER_TPM,
  nrow(annotation_filt)
))

# Transcript x tissue matrix (ID as row names)
tpm_expressed <- annotation_filt %>%
  select(ID, all_of(tissue_cols)) %>%
  column_to_rownames("ID") %>%
  as.matrix()

# Log2 transformation: log2(TPM + 1)
log_tpm <- log2(tpm_expressed + 1)

# ------------------------------------------------------------------------------
# 2. SAVE HEATMAPS TO PDF
# ------------------------------------------------------------------------------

# Heatmap 1: Absolute Log2(TPM + 1) Expression
pdf(file.path(out_dir, "heatmap_absolute_expression.pdf"), width = 7, height = 8, family = font_family)
pheatmap(
  log_tpm,
  color = colorRampPalette(c("white", "orange", "darkred"))(100),
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = sprintf("lncRNA Expression Across Tissues (log2 TPM+1)\n%g < max TPM ≤ %g", LOWER_TPM, UPPER_TPM),
  fontsize = 11,
  fontsize_col = 10,
  fontfamily = font_family
)
dev.off()

# Heatmap 2: Relative Expression / Row Z-score
pdf(file.path(out_dir, "heatmap_relative_zscore.pdf"), width = 7, height = 8, family = font_family)
pheatmap(
  log_tpm,
  scale = "row",
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Relative lncRNA Expression Across Tissues (Row Z-score)",
  fontsize = 11,
  fontsize_col = 10,
  fontfamily = font_family
)
dev.off()

# ------------------------------------------------------------------------------
# 3. SAVE HISTOGRAMS & DISTRIBUTION PLOTS TO PDF
# ------------------------------------------------------------------------------

# Define common Times New Roman theme for ggplot2
times_theme <- theme_minimal(base_family = font_family, base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, family = font_family),
    axis.title = element_text(size = 12, family = font_family),
    axis.text  = element_text(size = 10, family = font_family),
    strip.text = element_text(size = 11, family = font_family)
  )

# Long format from the FILTERED transcripts. A transcript is only counted in a
# tissue if its TPM in that tissue is above the lower cutoff (otherwise the
# zeros produce a huge spike at 0 that hides the real distribution).
tpm_long <- annotation_filt %>%
  select(ID, Category, all_of(tissue_cols)) %>%
  # select() with the named vector already renamed the columns to the tissue names
  pivot_longer(
    cols = all_of(names(tissue_cols)),
    names_to = "Tissue",
    values_to = "TPM"
  ) %>%
  mutate(log2_TPM = log2(TPM + 1)) %>%
  filter(TPM > LOWER_TPM)

subtitle_txt <- sprintf("Expressed lncRNAs: %g < TPM; transcripts with max TPM > %g excluded", LOWER_TPM, UPPER_TPM)

# Plot 1: Faceted Histograms per Tissue Type
p1 <- ggplot(tpm_long, aes(x = log2_TPM)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  facet_wrap(~ Tissue, scales = "free_y") +
  times_theme +
  labs(
    title = "Distribution of lncRNA Expression by Tissue Type",
    subtitle = subtitle_txt,
    x = expression(log[2](TPM + 1)),
    y = "Count"
  )

ggsave(file.path(out_dir, "histogram_expression_by_tissue.pdf"), plot = p1, width = 8, height = 6, device = cairo_pdf)

# Plot 2: Density Ridges across Tissues
p2 <- ggplot(tpm_long, aes(x = log2_TPM, y = Tissue, fill = Tissue)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) +
  times_theme +
  theme(legend.position = "none") +
  labs(
    title = "lncRNA Expression Density Across Tissues",
    subtitle = subtitle_txt,
    x = expression(log[2](TPM + 1)),
    y = "Tissue"
  )

ggsave(file.path(out_dir, "density_ridges_by_tissue.pdf"), plot = p2, width = 7, height = 6, device = cairo_pdf)

# Plot 3: Histograms grouped by Genomic Category across Tissues
p3 <- ggplot(tpm_long, aes(x = log2_TPM, fill = Category)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.5) +
  facet_wrap(~ Tissue) +
  times_theme +
  labs(
    title = "lncRNA Expression by Tissue and Genomic Category",
    subtitle = subtitle_txt,
    x = expression(log[2](TPM + 1)),
    y = "Count"
  )

ggsave(file.path(out_dir, "histogram_by_category_and_tissue.pdf"), plot = p3, width = 9, height = 7, device = cairo_pdf)

# ------------------------------------------------------------------------------
# 4. SHARED lncRNA GENES ACROSS TISSUES (VENN / UPSET)
# ------------------------------------------------------------------------------
# Gene-level view: a lncRNA gene (lncRNA_gene column, e.g. MSTRG.334) has several
# isoforms (MSTRG.334.1, .2, ...). Gene expression in a tissue = sum of its
# retained isoforms' mean TPM; the gene is "detected" in that tissue if the
# sum is >= LOWER_TPM. Genes removed by the upper cutoff are excluded entirely
# (only genes with at least one retained transcript are considered).

gene_tpm <- annotation_filt %>%
  group_by(lncRNA_gene) %>%
  summarise(across(all_of(unname(tissue_cols)), sum), .groups = "drop")
colnames(gene_tpm) <- c("lncRNA_gene", names(tissue_cols))

gene_sets <- lapply(names(tissue_cols), function(tis) {
  gene_tpm$lncRNA_gene[gene_tpm[[tis]] >= LOWER_TPM]
})
names(gene_sets) <- names(tissue_cols)

message("Genes detected per tissue:")
print(sapply(gene_sets, length))

# Make ggVennDiagram's text labels use Times New Roman too
update_geom_defaults("text",  list(family = font_family))
update_geom_defaults("label", list(family = font_family))

p_venn <- ggVennDiagram(gene_sets, label = "count", label_alpha = 0, set_size = 5, label_size = 3.5) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Genes") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  theme(text = element_text(family = font_family),
        plot.title = element_text(face = "bold", size = 14, family = font_family)) +
  labs(title = "lncRNA Genes Detected Across Tissues",
       subtitle = sprintf("Gene detected if summed isoform TPM ≥ %g; genes with transcripts > %g TPM excluded", LOWER_TPM, UPPER_TPM))

ggsave(file.path(out_dir, "venn_lncRNA_genes_by_tissue.pdf"), plot = p_venn, width = 9, height = 7, device = cairo_pdf)

# A 6-set Venn is crowded for slides; the UpSet plot shows the same intersections
# as a sorted bar chart and is usually easier to read.
upset_file <- file.path(out_dir, "upset_lncRNA_genes_by_tissue.pdf")
pdf(upset_file, width = 9, height = 6, family = font_family)
# newpage = FALSE stops UpSetR from calling grid.newpage(), which otherwise
# leaves a blank first page in the PDF
print(upset(
  fromList(gene_sets),
  nsets = length(gene_sets),
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  text.scale = c(1.4, 1.3, 1.2, 1.1, 1.4, 1.2),
  mainbar.y.label = "Shared lncRNA genes",
  sets.x.label = "Genes per tissue"
), newpage = FALSE)
dev.off()

# UpSetR leaves a blank first page in PDF output; keep only the last page.
# (Needs: install.packages("pdftools"))
n_pages <- pdftools::pdf_length(upset_file)
if (n_pages > 1) {
  tmp_pdf <- tempfile(fileext = ".pdf")
  pdftools::pdf_subset(upset_file, pages = n_pages, output = tmp_pdf)
  file.copy(tmp_pdf, upset_file, overwrite = TRUE)
}

# Table of which tissues each gene is detected in (so shared IDs can be pulled out)
gene_membership <- gene_tpm %>%
  mutate(across(all_of(names(tissue_cols)), ~ . >= LOWER_TPM, .names = "det_{.col}")) %>%
  rowwise() %>%
  mutate(
    n_tissues = sum(c_across(starts_with("det_"))),
    tissues   = paste(names(tissue_cols)[c_across(starts_with("det_"))], collapse = ",")
  ) %>%
  ungroup() %>%
  filter(n_tissues > 0) %>%
  arrange(desc(n_tissues), tissues)

write.table(gene_membership, file.path(out_dir, "lncRNA_gene_tissue_membership.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. BAR CHARTS: lncRNA CLASSIFICATION AND ASSOCIATED TEs
# ------------------------------------------------------------------------------
# Counts are per TRANSCRIPT (rows of the annotation table), shown for all predicted
# transcripts and for the expressed set that passes both TPM thresholds.
# "NA" in Category (e.g. convergent/intronic cases the classifier left unassigned)
# is shown as "Unclassified".

set_all  <- "All predicted transcripts"
set_expr <- sprintf("Expressed (%g < TPM ≤ %g)", LOWER_TPM, UPPER_TPM)
set_cols <- setNames(c("grey65", "steelblue"), c(set_all, set_expr))

class_df <- bind_rows(
  annotation_data %>% mutate(set = set_all),
  annotation_filt %>% mutate(set = set_expr)
) %>%
  mutate(
    set      = factor(set, levels = c(set_all, set_expr)),
    Category = replace_na(Category, "Unclassified"),
    TE_types = replace_na(TE_types, "None")
  )

# 5a. Classification features (one panel each)
class_counts <- class_df %>%
  select(set, type, direction, subtype, location, Category) %>%
  pivot_longer(-set, names_to = "feature", values_to = "class") %>%
  mutate(feature = factor(
    feature,
    levels = c("type", "direction", "subtype", "location", "Category"),
    labels = c("Genic / intergenic", "Strand vs. partner gene", "Subtype", "Location", "Category")
  )) %>%
  count(set, feature, class)

p_class <- ggplot(class_counts, aes(x = class, y = n, fill = set)) +
  geom_col(position = position_dodge(width = 0.9), color = "black", linewidth = 0.2) +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 2.8) +
  facet_wrap(~ feature, scales = "free_x", nrow = 2) +
  scale_fill_manual(values = set_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  times_theme +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "top") +
  labs(title = "lncRNA Classification",
       x = NULL, y = "Number of transcripts")

ggsave(file.path(out_dir, "barchart_lncRNA_classification.pdf"), plot = p_class,
       width = 12, height = 8, device = cairo_pdf)

# 5b. Overlapping TE classes (a transcript with several TE classes counts once per class)
te_counts <- class_df %>%
  select(set, ID, TE_types) %>%
  separate_rows(TE_types, sep = ";") %>%
  count(set, TE_types) %>%
  mutate(TE_types = fct_reorder(TE_types, n, .fun = sum))

p_te <- ggplot(te_counts, aes(x = TE_types, y = n, fill = set)) +
  geom_col(position = position_dodge(width = 0.9), color = "black", linewidth = 0.2) +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_manual(values = set_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  times_theme +
  theme(legend.position = "top") +
  labs(title = "Transposable Element Classes Overlapping lncRNAs",
       subtitle = "\"None\" = no TE overlap",
       x = "TE class", y = "Number of transcripts")

ggsave(file.path(out_dir, "barchart_lncRNA_TE_types.pdf"), plot = p_te,
       width = 8, height = 6, device = cairo_pdf)

# 5c. Fraction of each Category that overlaps a TE
te_by_cat <- class_df %>%
  mutate(TE_status = ifelse(TE_types == "None", "No TE overlap", "TE overlap")) %>%
  count(set, Category, TE_status)

p_te_cat <- ggplot(te_by_cat, aes(x = Category, y = n, fill = TE_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.2) +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 3) +
  facet_wrap(~ set) +
  scale_fill_manual(values = c("No TE overlap" = "grey80", "TE overlap" = "darkorange"), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  times_theme +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "top") +
  labs(title = "TE Overlap by lncRNA Category",
       x = NULL, y = "Proportion of transcripts")

ggsave(file.path(out_dir, "barchart_TE_overlap_by_category.pdf"), plot = p_te_cat,
       width = 10, height = 6, device = cairo_pdf)