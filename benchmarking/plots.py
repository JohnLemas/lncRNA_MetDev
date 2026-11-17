#!/usr/bin/env python3
"""Summary figures: recall vs T_C per call set, stage-wise recall, UpSet of tool agreement."""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True)
ap.add_argument("--species", nargs="+", required=True)
ap.add_argument("--tools", nargs="+", required=True)
a = ap.parse_args()

pdir = os.path.join(a.results, "summary", "plots")
os.makedirs(pdir, exist_ok=True)
perf = pd.read_csv(f"{a.results}/summary/performance_all_species.tsv", sep="\t")

for thr in perf.threshold.unique():
    d = perf[(perf.threshold == thr) & (perf.label == "T")]

    # 1. recall of the custom truth (T_C) by tool / final
    order = list(d[d.kind.isin(["tool", "filter", "final"])].call_set.drop_duplicates())
    pv = d[d.call_set.isin(order)].pivot(index="call_set", columns="species",
                                         values="recall_vs_ceiling").reindex(order)
    ax = pv.plot.bar(figsize=(8, 4))
    ax.set_ylabel("Fraction of T_C recovered")
    ax.set_xlabel("")
    ax.set_title(f"Recall against custom truth (T_C), threshold: {thr}")
    plt.tight_layout()
    plt.savefig(f"{pdir}/recall_vs_TC_{thr}.png", dpi=200)
    plt.close()

    # 2. stage-wise loss
    order = list(d[d.kind.isin(["baseline", "cumulative"])].call_set.drop_duplicates())
    pv = d[d.call_set.isin(order)].pivot(index="call_set", columns="species",
                                         values="recall_vs_ceiling").reindex(order)
    ax = pv.plot(marker="o", figsize=(8, 4))
    ax.set_ylabel("Fraction of T_C retained")
    ax.set_xlabel("")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=30, ha="right")
    ax.set_title(f"Stage-wise recall of T_C, threshold: {thr}")
    plt.tight_layout()
    plt.savefig(f"{pdir}/stagewise_recall_{thr}.png", dpi=200)
    plt.close()

# 3. UpSet of tool agreement (optional dependency)
try:
    from upsetplot import UpSet, from_indicators
    for sp in a.species:
        fl = pd.read_csv(f"{a.results}/{sp}/metrics/candidate_flags.tsv", sep="\t", index_col=0)
        data = from_indicators(a.tools, data=fl[a.tools].astype(bool))
        UpSet(data, subset_size="count", show_counts=True).plot()
        plt.savefig(f"{pdir}/upset_{sp}.png", dpi=200, bbox_inches="tight")
        plt.close()
except ImportError:
    print("upsetplot not installed; skipping UpSet plots (pip install upsetplot)")
