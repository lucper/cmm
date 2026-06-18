#!/usr/bin/env python3
"""
Grouped box plot of per-(dataset, method) normalized score distributions from
030-score files.

X-axis: datasets. Within each dataset, one box per method (color-coded).

Input is one or more 'dataset:method:path' triples, where each path points at a
030-score output file. Those files are ALREADY normalized: column 3 is the
per-rank ratio (slider_score / exact_score * 100). This script reads column 3
directly and plots the distribution -- no baseline pairing is done here, because
030-score already performed that division.

(The earlier version of this script read raw 025-dedup solution files and did
the baseline normalization itself. That responsibility now lives upstream in
030-score, so the normalization machinery has been removed.)

The dataset and method names cannot contain colons.
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns
from pathlib import Path


def parse_input_spec(spec):
    """Split 'dataset:method:path' on the first two colons.
    Path may contain colons (only first two ':' are separators)."""
    parts = spec.split(":", 2)
    if len(parts) != 3:
        raise ValueError(f"Expected 'dataset:method:path', got: {spec}")
    return parts[0], parts[1], parts[2]


def read_scores(path):
    """Read a 030-score file. Returns the list of column-3 values (the
    normalized ratio) in file order. Lines whose third token isn't parseable
    as a float are skipped -- handles a 'slider cmm ratio' (or 'X Y x2')
    header transparently."""
    scores = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                score = float(parts[2])
            except ValueError:
                continue
            scores.append(score)
    return scores


def load_long_dataframe(specs):
    """specs: list of (dataset, method, path). Returns a long-form DataFrame
    with columns (dataset, method, score), plus first-seen orderings for both
    categorical axes (so the plot respects the CLI input order rather than
    sorting alphabetically)."""
    datasets_order = []
    methods_order  = []
    rows = []
    for dataset, method, path in specs:
        if dataset not in datasets_order: datasets_order.append(dataset)
        if method  not in methods_order:  methods_order.append(method)
        scores = read_scores(path)
        if not scores:
            print(f"Warning: no rows read from {path}")
            continue
        for s in scores:
            rows.append({"dataset": dataset, "method": method, "score": s})
        print(f"Loaded {len(scores)} scores from {path} "
              f"(dataset={dataset}, method={method})")
    df = pd.DataFrame(rows)
    if df.empty:
        return df, [], []
    df["dataset"] = pd.Categorical(df["dataset"],
                                   categories=datasets_order, ordered=True)
    df["method"]  = pd.Categorical(df["method"],
                                   categories=methods_order, ordered=True)
    return df, datasets_order, methods_order


def plot_box(df, datasets, methods, output_path, score_label,
             reference_line=None, reference_label=None):
    n_methods = len(methods)

    palette = sns.color_palette("tab10", n_colors=max(n_methods, 1))

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)

    sns.boxplot(data=df, x="dataset", y="score", hue="method",
                ax=ax, palette=palette, width=0.8, gap=0.05,
                flierprops={"marker": "o", "markersize": 3, "alpha": 0.5},
                showmeans=True,
                meanprops={"marker": "D", "markerfacecolor": "white",
                           "markeredgecolor": "black", "markersize": 5})

    if reference_line is not None:
        ax.axhline(reference_line, color="red", linestyle="--",
                   linewidth=2.0, alpha=0.8, zorder=0)
        # Keep 0 and the reference line visible even when all data is below.
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(bottom=min(ymin, 0),
                    top=max(ymax, reference_line * 1.05))
        # Y ticks formatted as percentages (0%, 20%, ..., 100%).
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))

    ax.set_xlabel("")
    ax.set_ylabel(score_label, fontsize=20)
    ax.tick_params(axis="both", labelsize=18)
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)

    # Legend: when a reference line is shown, build a custom legend that names
    # both the box color(s) and the red reference line. Otherwise fall back to
    # the default per-method legend (or remove it entirely for a single method).
    if reference_line is not None and reference_label is not None:
        handles = [Patch(facecolor=palette[i], edgecolor="black", label=m)
                   for i, m in enumerate(methods)]
        handles.append(Line2D([0], [0], color="red", linestyle="--",
                              linewidth=2.0, label=reference_label))
        ax.legend(handles=handles, fontsize=14, ncol=len(handles),
                  loc="lower center", bbox_to_anchor=(0.5, 1.02),
                  frameon=False)
    elif n_methods > 1:
        ax.legend(fontsize=14, ncol=n_methods,
                  loc="lower center", bbox_to_anchor=(0.5, 1.02),
                  frameon=False)
    elif ax.get_legend() is not None:
        ax.get_legend().remove()
    sns.despine(ax=ax)

    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, nargs="+",
                        help="One or more 'dataset:method:path' triples, each "
                             "path a 030-score file (column 3 = normalized "
                             "ratio). Order determines x-axis (datasets) and "
                             "legend (methods) ordering.")
    parser.add_argument("--output", required=True,
                        help="Output stem (e.g. score_box)")
    parser.add_argument("--score-label", default="Score / exact (%)",
                        help="Y-axis label.")
    parser.add_argument("--reference-label", default="Exact algorithm",
                        help="Legend label for the red 100%% reference line.")
    parser.add_argument("--no-reference-line", action="store_true",
                        help="Suppress the red 100%% reference line and the "
                             "percent y-axis formatting (use if column 3 is "
                             "not a percentage).")
    args = parser.parse_args()

    sns.set_theme(style="white", font_scale=1.0)

    specs = [parse_input_spec(s) for s in args.input]
    df, datasets, methods = load_long_dataframe(specs)
    if df.empty:
        raise SystemExit("Error: no data read from any input file")

    reference_line = None if args.no_reference_line else 100.0

    stem = str(Path(args.output).with_suffix(""))
    plot_box(df, datasets, methods,
             output_path=f"{stem}.pdf",
             score_label=args.score_label,
             reference_line=reference_line,
             reference_label=args.reference_label)


if __name__ == "__main__":
    main()
