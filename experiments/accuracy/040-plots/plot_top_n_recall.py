#!/usr/bin/env python3
"""
Plot how the quality of the heuristic's top-N motif pairs varies with some
budget parameter (max time per run, number of runs, heuristic output size k,
etc.). The y-axis is the percentage of those top N whose best match in the
reference method has similarity at least s*. The x-axis is the chosen budget
variable.

Input: cmm_eval best-match TSVs, one per (dataset, x-value). For each TSV,
solution A should be the heuristic output (or the aggregated heuristic output
in the multi-run case) and solution B should be the reference, so each TSV
row is a heuristic pair with its best match in the reference. The first
--top-n rows of each TSV (= heuristic's top N by score) are evaluated.

Set --x-label to describe what the x-axis represents:
  --x-label 'Max time (min)'               # RQ1: vary max time per run
  --x-label 'Number of runs'               # RQ2: vary number of aggregated runs
  --x-label '$k$ (heuristic output size)'  # vary heuristic's truncation k

Curve shapes worth knowing:
  - Rising curve: extra budget buys the heuristic better exploration, so
    its top N approximates the reference better.
  - Flat curve: budget is just a truncation parameter; the heuristic's top
    N is unchanged by spending more.
  - Saturating curve: rising at first, then plateauing -- the heuristic has
    a reachable space it converges into and can't go further.
"""

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path

## expected columns from output TSV of 030-eval
COLS = ["rank", "X", "Y", "score",
        "best_X", "best_Y", "best_score", "similarity"]

def parse_input_spec(spec):
    """Split 'dataset:x:path' on the first two colons. Path may contain ':'.
    The middle x token must be numeric (int or float)."""
    parts = spec.split(":", 2)
    if len(parts) != 3:
        raise ValueError(f"Expected 'dataset:x:path', got: {spec}")
    dataset, x_str, path = parts
    try:
        x = float(x_str)
    except ValueError:
        raise ValueError(f"x must be numeric, got: '{x_str}' in {spec}")
    if x == int(x):
        x = int(x)  # keep tidy integer formatting on the x-axis
    return dataset, x, path

def compute_recall(df, top_n, threshold):
    """Take the first top_n rows of df, return fraction whose 'similarity'
    column meets or exceeds threshold."""
    top = df.head(top_n)
    n = len(top)
    if n < top_n:
        print(f"Warning: requested top-{top_n} but only {n} rows available")
    if n == 0:
        return 0.0
    matched = int((top["similarity"] >= threshold).sum())
    return matched / n

def plot(specs, output_path, top_n, threshold, log_scale, x_label):
    # Group by dataset, preserving first-seen order in specs.
    by_dataset = {}
    for dataset, x, path in specs:
        df = pd.read_csv(path, sep="\t", header=None, names=COLS)
        r = compute_recall(df, top_n, threshold)
        by_dataset.setdefault(dataset, []).append((x, r))
        print(f"{dataset} x={x}: {100*r:.1f}% from {path}")

    palette = sns.color_palette("tab10", n_colors=max(len(by_dataset), 3))
    markers = ["o", "s", "^", "D", "v", "P", "*"]  # cycled per dataset

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)

    all_xs = sorted({x for points in by_dataset.values() for x, _ in points})

    for idx, ((dataset, points), color) in enumerate(zip(by_dataset.items(), palette)):
        points.sort()  # ascending x
        xs   = [x for x, _ in points]
        pcts = [100.0 * r for _, r in points]
        ax.plot(xs, pcts, marker=markers[idx % len(markers)],
                linewidth=2.0, markersize=8,
                label=dataset, color=color)

    if log_scale:
        ax.set_xscale("log")
        ax.set_xticks(all_xs)
        ax.set_xticklabels([str(x) for x in all_xs])
        ax.minorticks_off()

    ax.set_xlabel(x_label, fontsize=20)
    ax.set_ylabel(f"% top {top_n} with $s_h \\geq$ {threshold}", fontsize=20)
    ax.set_ylim(-2, 102)
    ax.tick_params(axis="both", labelsize=18)
    ax.grid(linestyle="--", linewidth=0.6, alpha=0.5)
    if len(by_dataset) > 1:
        ax.legend(fontsize=14, loc="best")
    sns.despine(ax=ax)

    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, nargs="+",
                        help="One or more 'dataset:x:path' triples. x is the "
                             "x-axis value for this data point (e.g., 10 for "
                             "max_time=10min, or 5 for n_runs=5). Path is the "
                             "cmm_eval TSV.")
    parser.add_argument("--output", required=True,
                        help="Output stem (e.g. rq1_recall)")
    parser.add_argument("--top-n", type=int, default=1000,
                        help="Top N rows to evaluate per TSV (default 1000, "
                             "matching the reference method's solution size)")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Similarity threshold s* (default 0.5)")
    parser.add_argument("--x-label", required=True,
                        help="X-axis label (LaTeX accepted). "
                             "E.g. 'Max time (min)' or 'Number of runs'.")
    parser.add_argument("--log-scale", action="store_true",
                        help="Use log scale on x-axis (default linear)")
    args = parser.parse_args()

    sns.set_theme(style="white", font_scale=1.0)

    specs = [parse_input_spec(s) for s in args.input]
    stem  = str(Path(args.output).with_suffix(""))
    plot(specs, f"{stem}.pdf",
         top_n=args.top_n, threshold=args.threshold,
         log_scale=args.log_scale, x_label=args.x_label)

if __name__ == "__main__":
    main()
