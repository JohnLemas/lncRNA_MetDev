#!/usr/bin/env python3
"""Parse BUSCO short_summary files from one or more run directories into one table."""
import argparse
import glob
import os
import re

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--species", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("dirs", nargs="+")
a = ap.parse_args()

rows = []
for d in a.dirs:
    mode = os.path.basename(d.rstrip("/"))
    for f in glob.glob(os.path.join(d, "**", "short_summary*.txt"), recursive=True):
        with open(f) as fh:
            for line in fh:
                s = line.strip()
                if s.startswith("C:"):
                    v = dict(re.findall(r"([A-Za-z]+):([\d.]+)", s))
                    rows.append(dict(species=a.species, mode=mode, file=os.path.basename(f),
                                     complete=v.get("C"), single=v.get("S"), duplicated=v.get("D"),
                                     fragmented=v.get("F"), missing=v.get("M"), n_buscos=v.get("n")))
                    break
pd.DataFrame(rows, columns=["species", "mode", "file", "complete", "single", "duplicated",
                            "fragmented", "missing", "n_buscos"]).to_csv(a.out, sep="\t", index=False)
