#!/usr/bin/env python3
"""Collapse a minimap2 PAF (run with -c) into one row per (query, target) pair.

Identity  = sum(residue matches) / sum(alignment block length) over the pair's alignments.
qcov/tcov = fraction of query/target covered by the UNION of aligned intervals.
scov      = coverage of the SHORTER of the two sequences (handles fragmented StringTie
            transcripts contained in longer reference entries, and vice versa).

No thresholds are applied here; compute_metrics.py applies them so thresholds can be
changed without re-running minimap2.
"""
import argparse


def union_len(intervals):
    total, cs, ce = 0, None, None
    for s, e in sorted(intervals):
        if cs is None:
            cs, ce = s, e
        elif s <= ce:
            ce = max(ce, e)
        else:
            total += ce - cs
            cs, ce = s, e
    if cs is not None:
        total += ce - cs
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paf")
    ap.add_argument("out")
    ap.add_argument("--same-strand", action="store_true",
                    help="only keep '+' alignments (sense-strand matches)")
    a = ap.parse_args()

    pairs = {}
    with open(a.paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            if a.same_strand and f[4] != "+":
                continue
            key = (f[0], f[5])
            r = pairs.get(key)
            if r is None:
                r = pairs[key] = {"ql": int(f[1]), "tl": int(f[6]),
                                  "qi": [], "ti": [], "m": 0, "b": 0}
            r["qi"].append((int(f[2]), int(f[3])))
            r["ti"].append((int(f[7]), int(f[8])))
            r["m"] += int(f[9])
            r["b"] += int(f[10])

    with open(a.out, "w") as o:
        o.write("q\tt\tqlen\ttlen\tidentity\tqcov\ttcov\tscov\n")
        for (q, t), r in pairs.items():
            qc = min(1.0, union_len(r["qi"]) / r["ql"])
            tc = min(1.0, union_len(r["ti"]) / r["tl"])
            ident = r["m"] / r["b"] if r["b"] else 0.0
            sc = qc if r["ql"] <= r["tl"] else tc
            o.write(f"{q}\t{t}\t{r['ql']}\t{r['tl']}\t{ident:.4f}\t{qc:.4f}\t{tc:.4f}\t{sc:.4f}\n")


if __name__ == "__main__":
    main()
