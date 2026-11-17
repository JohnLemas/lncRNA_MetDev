#!/usr/bin/env python3
"""Per-species benchmarking metrics.

Sets
  C  candidates            T  PLncDB (cd-hit condensed, all sources)
  A  annotated lncRNAs     N  annotated mRNAs (presumed negatives)
  E  EVLncDB subset of PLncDB (kept aside, matched independently)
  T_C = T entries matched by >=1 candidate (the ceiling given your input data)

Call sets: each tool, DIAMOND-pass, final (all tools AND DIAMOND-pass),
cumulative stage-wise intersections, and ">= k tools agree".

Two views of every call set:
  * entry level     : unique reference entries recovered (recall vs all of T, and vs T_C)
  * candidate level : calls matching a reference (precision/recall/specificity within the
                      labelled subset: positives = candidates matching the label set,
                      negatives = candidates matching an annotated mRNA but no lncRNA set)
"""
import argparse
import os

import numpy as np
import pandas as pd


def fasta_ids(path):
    with open(path) as fh:
        return [l[1:].split()[0] for l in fh if l.startswith(">")]


def read_ids(path):
    with open(path) as fh:
        return {l.split()[0] for l in fh if l.strip()}


def diamond_hits(path, evalue):
    hits = set()
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return hits
    with open(path) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) >= 11 and float(f[10]) <= evalue:  # outfmt 6: col 11 = evalue
                hits.add(f[0])
    return hits


def div(a, b):
    return a / b if b else np.nan


def kv(items):
    return dict(i.split("=", 1) for i in items)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--tool", action="append", required=True, help="name=path (IDs called lncRNA)")
    ap.add_argument("--diamond", required=True)
    ap.add_argument("--diamond-evalue", type=float, default=1e-5)
    ap.add_argument("--fasta", action="append", required=True, help="T|A|N|E=path")
    ap.add_argument("--pairs", action="append", required=True, help="C_vs_T=path ... A_vs_T=path")
    ap.add_argument("--threshold", action="append", required=True, help="name:min_identity:min_scov")
    ap.add_argument("--stage-order", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    sp = args.species

    cands = fasta_ids(args.candidates)
    cand_set = set(cands)
    tools = {k: read_ids(v) & cand_set for k, v in kv(args.tool).items()}
    tool_names = list(tools)
    ref_total = {k: len(set(fasta_ids(v))) for k, v in kv(args.fasta).items()}
    pairs_raw = {k: pd.read_csv(v, sep="\t") for k, v in kv(args.pairs).items()}
    thresholds = []
    for t in args.threshold:
        n, mi, ms = t.split(":")
        thresholds.append((n, float(mi), float(ms)))

    # ---- call sets -------------------------------------------------------------
    dia_hit = diamond_hits(args.diamond, args.diamond_evalue) & cand_set
    dia_pass = cand_set - dia_hit
    final = set.intersection(*tools.values()) & dia_pass

    call_sets = {"candidates": ("baseline", cand_set)}
    for t in tool_names:
        call_sets[t] = ("tool", tools[t])
    call_sets["diamond_pass"] = ("filter", dia_pass)
    call_sets["final"] = ("final", final)

    cum = None
    for i, st in enumerate(args.stage_order, 1):
        s = dia_pass if st == "diamond" else tools[st]
        cum = set(s) if cum is None else cum & s
        call_sets[f"cum{i}_" + "+".join(args.stage_order[:i])] = ("cumulative", set(cum))

    ntools = {c: sum(c in tools[t] for t in tool_names) for c in cands}
    for k in range(1, len(tool_names) + 1):
        call_sets[f"ge{k}_tools"] = ("n_tools", {c for c, n in ntools.items() if n >= k})

    # ---- flag table ------------------------------------------------------------
    flags = pd.DataFrame(index=pd.Index(cands, name="candidate"))
    for t in tool_names:
        flags[t] = flags.index.isin(tools[t])
    flags["diamond_hit"] = flags.index.isin(dia_hit)
    flags["n_tools"] = pd.Series(ntools)
    flags["final"] = flags.index.isin(final)

    # ---- metrics ---------------------------------------------------------------
    rows, a_vs_t_rows = [], []
    for thr, mi, ms in thresholds:
        flt = {}
        for r in "TANE":
            p = pairs_raw[f"C_vs_{r}"]
            flt[r] = p[(p.identity >= mi) & (p.scov >= ms)]
        cm = {r: set(flt[r].q) for r in flt}
        ceiling = {r: flt[r].t.nunique() for r in "TAE"}   # T_C, A_C, E_C
        for r in "TANE":
            flags[f"in_{r}_{thr}"] = flags.index.isin(cm[r])
        flags[f"in_T_C_{thr}"] = flags[f"in_T_{thr}"]  # candidate belongs to custom truth
        pos = {"T": cm["T"], "A": cm["A"], "TA": cm["T"] | cm["A"], "E": cm["E"]}
        neg = cm["N"] - (cm["T"] | cm["A"] | cm["E"])

        for name, (kind, S) in call_sets.items():
            for label in ["T", "A", "TA", "E"]:
                P_ = pos[label]
                tp, fp = len(S & P_), len(S & neg)
                fn, tn = len(P_) - tp, len(neg) - fp
                row = dict(species=sp, threshold=thr, min_identity=mi, min_scov=ms,
                           call_set=name, kind=kind, label=label,
                           n_calls=len(S),
                           n_calls_matching_label=tp,
                           n_calls_matching_N_only=fp,
                           n_calls_unlabelled=len(S) - tp - fp,
                           n_labelled_pos_total=len(P_), n_labelled_neg_total=len(neg),
                           precision=div(tp, tp + fp),
                           cand_recall=div(tp, len(P_)),
                           specificity=div(tn, tn + fp),
                           F1=div(2 * tp, 2 * tp + fp + fn))
                if label != "TA":
                    h = flt[label][flt[label].q.isin(S)]
                    rec = h.t.nunique()
                    row.update(n_ref_total=ref_total[label],
                               n_ref_ceiling=ceiling[label],      # T_C / A_C / E_C
                               n_ref_recovered=rec,
                               recall_vs_all_ref=div(rec, ref_total[label]),
                               recall_vs_ceiling=div(rec, ceiling[label]))
                rows.append(row)

        # layer 3: annotated lncRNAs vs PLncDB truth
        p = pairs_raw["A_vs_T"]
        p = p[(p.identity >= mi) & (p.scov >= ms)]
        a_vs_t_rows.append(dict(species=sp, threshold=thr,
                                n_A_total=ref_total["A"], n_A_in_T=p.q.nunique(),
                                n_T_total=ref_total["T"], n_T_matched_by_A=p.t.nunique(),
                                frac_A_in_T=div(p.q.nunique(), ref_total["A"])))

    perf = pd.DataFrame(rows)
    perf.to_csv(f"{args.outdir}/performance.tsv", sep="\t", index=False)

    tc = perf[perf.label == "T"][[
        "species", "threshold", "call_set", "kind", "n_calls", "n_calls_matching_label",
        "n_ref_ceiling", "n_ref_recovered", "recall_vs_ceiling", "recall_vs_all_ref",
        "precision", "specificity", "F1"]].rename(columns={
            "n_calls_matching_label": "n_calls_in_T_C",
            "n_ref_ceiling": "T_C_size",
            "n_ref_recovered": "T_C_entries_recovered"})
    tc.to_csv(f"{args.outdir}/tc_benchmark.tsv", sep="\t", index=False)

    pd.DataFrame(a_vs_t_rows).to_csv(f"{args.outdir}/annot_vs_truth.tsv", sep="\t", index=False)
    flags.to_csv(f"{args.outdir}/candidate_flags.tsv", sep="\t")

    sizes = [dict(species=sp, set="candidates_C", n=len(cands))]
    sizes += [dict(species=sp, set=f"ref_{k}", n=v) for k, v in ref_total.items()]
    sizes += [dict(species=sp, set=f"call_{k}", n=len(v[1])) for k, v in call_sets.items()]
    pd.DataFrame(sizes).to_csv(f"{args.outdir}/set_sizes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
