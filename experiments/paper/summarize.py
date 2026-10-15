"""
Summarize the CSV output of the paper's simulation study as Markdown tables.

Usage (from the repository root):

    python3 experiments/paper/summarize.py results/paper

Prints mean (SD) across trajectories for every experiment whose CSV is
present. Metrics are first averaged over time within a trajectory (by the
experiment scripts), then summarized across trajectories.
"""

import csv
import glob
import math
import os
import statistics
import sys


def read(pattern):
    rows = []
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            rows.extend(csv.DictReader(f))
    return rows


def values(rows, key):
    out = []
    for r in rows:
        try:
            x = float(r[key])
        except (KeyError, ValueError):
            continue
        if math.isfinite(x):
            out.append(x)
    return out


def cell(rows, key, digits=2):
    v = values(rows, key)
    if not v:
        return "–"
    if len(v) == 1:
        return f"{v[0]:.{digits}f}"
    return f"{statistics.mean(v):.{digits}f} ({statistics.stdev(v):.{digits}f})"


def table(title, groups, columns):
    """groups: list of (label, rows); columns: list of (header, key, digits)."""
    groups = [(label, rows) for label, rows in groups if rows]
    if not groups:
        return
    print(f"\n#### {title}\n")
    print("| | n | " + " | ".join(h for h, _, _ in columns) + " |")
    print("|---|---|" + "---|" * len(columns))
    for label, rows in groups:
        seeds = len({r["seed"] for r in rows})
        print(f"| {label} | {seeds} | " +
              " | ".join(cell(rows, key, d) for _, key, d in columns) + " |")


PRECISION = [
    ("RV", "rv", 2), ("RFE", "rfe", 2), ("MSE", "mse", 3), ("q", "norm_ratio", 2),
    ("RFE adj.", "rfe_adjusted", 2), ("RV off-diag.", "offdiag_rv", 2),
    ("RV partial corr.", "partial_rv", 2),
]
BLOCK = [("block RV", "block_rv", 2), ("block RFE", "block_rfe", 2),
         ("block MSE", "block_mse", 3)]
SHOCKS = [
    ("precision", "precision", 2), ("recall", "recall", 2),
    ("specificity", "specificity", 2), ("F1", "f1", 2), ("AUROC", "auroc", 2),
    ("Brier", "brier", 3), ("γ true shocks", "gamma_true_shocks", 2),
    ("γ non-shocks", "gamma_non_shocks", 2),
]


def by_method(rows, label=""):
    return [(f"{label}{m}", [r for r in rows if r["method"] == m]) for m in ("em", "mcmc")]


def main(directory):
    p = lambda name: os.path.join(directory, name)

    for dim in (10, 50):
        rows = read(p(f"experiment1_p{dim}*.csv"))
        rows = [r for r in rows if "t" not in r or r.get("mean_gamma") is None]
        table(f"Experiment I, p = {dim}: precision recovery", by_method(rows),
              PRECISION)
        table(f"Experiment I, p = {dim}: mean γ at the changepoint", by_method(rows),
              [("γ at t* = 26", "gamma_at_change", 2), ("γ elsewhere", "gamma_elsewhere", 2)])

    for dim in (10, 50):
        rows = [r for r in read(p(f"experiment2_p{dim}.csv")) +
                read(p(f"experiment2_p{dim}_s*.csv")) if "mean_gamma" not in r]
        table(f"Experiment II, p = {dim}: precision recovery (full active sets)",
              by_method(rows), PRECISION)
        table(f"Experiment II, p = {dim}: common node-time block", by_method(rows),
              BLOCK + [("block length", "block_end", 0), ("block vars", "block_vars", 1)])
        table(f"Experiment II, p = {dim}: continuation shocks", by_method(rows), SHOCKS)

    full = read(p("experiment2_p10.csv"))
    ablation = read(p("experiment2_p10_ablation.csv"))
    groups = []
    for label, rows in (("full", full), ("a_b = a_d = 0", ablation)):
        groups += by_method(rows, label + ", ")
    for subset in ("event_", "nonevent_"):
        table(f"Ablation, p = 10: shocks at {'event' if subset == 'event_' else 'non-event'} times",
              groups,
              [(h, subset + k, d) for h, k, d in SHOCKS[:6]])

    rank_groups = []
    for k in (2, 5, 8, 12):
        rows = full if k == 5 else read(p(f"experiment2_p10_k{k}.csv"))
        rank_groups += by_method(rows, f"k_fit = {k}, ")
    table("Rank sensitivity, p = 10, k_true = 5", rank_groups,
          [("RV", "rv", 2), ("RV off-diag.", "offdiag_rv", 2), ("RFE", "rfe", 2),
           ("AUROC", "auroc", 2), ("recall", "recall", 2), ("Brier", "brier", 3),
           ("entropy H_γ", "entropy", 3), ("γ near 0/1", "fraction_extreme", 2)])

    mech_groups = []
    for mech in ("history", "stochastic", "contemporaneous"):
        mech_groups += by_method(read(p(f"experiment3_{mech}.csv")), f"{mech}, ")
    table("Experiments III and IV: precision recovery", mech_groups,
          [("exit rate", "exit_rate", 3), ("RV", "rv", 2), ("RV off-diag.", "offdiag_rv", 2),
           ("RV partial corr.", "partial_rv", 2), ("q", "norm_ratio", 2),
           ("RFE adj.", "rfe_adjusted", 2)])
    table("Experiments III and IV: common block, selection diagnostics", mech_groups,
          [("block RV", "block_rv", 2), ("block RV, marginal", "blockmarg_rv", 2),
           ("block RFE, marginal", "blockmarg_rfe", 2), ("Δ bias", "delta_bias", 3),
           ("Δ relative bias", "delta_relative_bias", 2),
           ("RV off-diag., far from boundary", "safe_offdiag_rv", 2),
           ("AUROC", "auroc", 2)])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "results/paper")
