# lncRNA prediction benchmarking workflow

Run: `snakemake -j 16` (needs snakemake, minimap2, gffread, TransDecoder, BUSCO, python with pandas/numpy/matplotlib; `upsetplot` optional).

## Sets
| Set | Contents |
|---|---|
| C | candidates entering the ML tools |
| P | final predictions = FEELnc ∩ LncFinder ∩ LncBoost ∩ CPAT ∩ DIAMOND-pass |
| T | PLncDB, all sources, cd-hit condensed |
| E | EVLncDB subset, kept aside and matched independently of T |
| A / N | reference-annotated lncRNAs / mRNAs (N = presumed negatives) |
| T_C | T entries matched by at least one candidate = the ceiling given your input data |

## Matching
minimap2 (`-x asm20 -c`, secondary alignments kept) of C against T, A, N, E, and of A against T. Pairs are
collapsed to identity and coverage (union of aligned intervals), with coverage taken on the shorter
sequence of the pair. Two thresholds are applied downstream (strict 0.95/0.90, relaxed 0.90/0.80), so changing
them never re-runs minimap2. Only sense-strand alignments count by default (`require_same_strand`).

E is matched directly rather than inferred from T, because cd-hit may have collapsed an EVLncDB entry into a
non-EVLncDB representative.

## Where each of your layers is (results/summary/)
Everything is in `performance_all_species.tsv`, filtered by `call_set` and `label` (T, A, TA, E), once per threshold.

| Layer | Rows |
|---|---|
| 1 candidates in T | `call_set=candidates`, `label=T`: `n_calls_matching_label` (candidates) and `n_ref_recovered` (= \|T_C\|) |
| 2 final predictions in T | `call_set=final`, `label=T` |
| 3 annotated lncRNAs in T | `annot_vs_truth_all_species.tsv` |
| 4 annotated lncRNAs in C | `call_set=candidates`, `label=A` |
| 5 predictions matching annotated | `call_set=final`, `label=A` |
| 6 tool stringency | each tool / `ge{k}_tools` / `final`, labels T and A; `cum*` rows give the stage-wise loss |
| 7 EVLncDB | same rows with `label=E` |
| 8 BUSCO | `busco_summary_all_species.tsv` |
| 9 T_C benchmark | `tc_benchmark_all_species.tsv` (`T_C_size`, `T_C_entries_recovered`, `recall_vs_ceiling`) |

Column meanings: `n_calls_matching_label` counts calls (candidate level). `n_ref_recovered` counts unique truth entries
(entry level). `recall_vs_ceiling` = recovered / \|T_C\|; `recall_vs_all_ref` = recovered / \|T\|.
`precision`, `specificity`, `F1` are computed inside the labelled subset only: positives = candidates matching the label set,
negatives = candidates matching an annotated mRNA but no lncRNA set. Calls in neither group are `n_calls_unlabelled`
(novel, not false positives). `candidate_flags.tsv` has one row per candidate with every call and match flag.

## Inputs
`config.yaml` per species: genome, lncRNA and mRNA GTFs, PLncDB and EVLncDB FASTAs, candidate FASTA, full transcriptome FASTA,
DIAMOND outfmt-6 table, BUSCO lineage, and one ID list per tool (IDs the tool called lncRNA/non-coding).
Candidate IDs must be identical across the FASTA, tool lists and DIAMOND table.

## Caveats
- BUSCO runs on the full transcriptome via TransDecoder (`--single_best_only`, no Pfam/BLAST retention). This underestimates
  completeness slightly, so the nucleotide `-m transcriptome` run is included as a cross-check.
- minimap2 can miss very short or divergent entries; if `A` vs `T` overlap looks surprisingly low, relax the thresholds.
- `identity` sums matches over all alignments of a pair; overlapping secondary alignments can inflate it slightly.
