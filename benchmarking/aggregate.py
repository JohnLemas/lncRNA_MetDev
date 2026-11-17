#!/usr/bin/env python3
"""Concatenate per-species tables into results/summary/*_all_species.tsv."""
import argparse
import os

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--results", required=True)
ap.add_argument("--species", nargs="+", required=True)
a = ap.parse_args()

out = os.path.join(a.results, "summary")
os.makedirs(out, exist_ok=True)

for name, sub in [("performance", "metrics"), ("tc_benchmark", "metrics"),
                  ("annot_vs_truth", "metrics"), ("set_sizes", "metrics"),
                  ("busco_summary", "busco")]:
    dfs = [pd.read_csv(f"{a.results}/{sp}/{sub}/{name}.tsv", sep="\t") for sp in a.species]
    pd.concat(dfs, ignore_index=True).to_csv(f"{out}/{name}_all_species.tsv", sep="\t", index=False)
