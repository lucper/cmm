#!/usr/bin/env python3
"""
Grouped box plot of per-(dataset, method) score distributions from main_algo
solution files (typically after deduplication via cmm_dedup).

X-axis: datasets. Within each dataset, one box per method (color-coded).
N (count of dedup'd pairs) is annotated above each box. Whiskers use Tukey's
1.5*IQR rule; points outside are drawn as outliers.

Input is one or more 'dataset:method:path' triples; the dataset and method
names cannot contain colons.
"""

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path

def parse_input_spec(spec):
    """Split 'dataset:method:path' on the first two colons.
    Path may contain colons (only first two ':' are separators)."""
    parts = spec.split(":", 2)
    if len(parts) != 3:
        raise ValueError(f"Expected 'dataset:method:path', got: {spec}")
    return parts[0], parts[1], parts[2]

def read_solution(path):
    """Read a raw main_algo solution file. Returns a list of scores in file
    order. Lines whose third token isn't parseable as a float are skipped --
    handles the 'X Y x2' header transparently."""
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
        scores = read_solution(path)
        if not scores:
            print(f"Warning: no rows read from {path}")
            continue
        for s in scores:
            rows.append({"dataset": dataset, "method": method, "score": s})
        print(f"Loaded {len(scores)} pairs from {path} "
              f"(dataset={dataset}, method={method})")
    df = pd.DataFrame(rows)
    df["dataset"] = pd.Categorical(df["dataset"],
                                   categories=datasets_order, ordered=True)
    df["method"]  = pd.Categorical(df["method"],
                                   categories=methods_order, ordered=True)
    return df, datasets_order, methods_order

def plot_box(df, datasets, methods, output_path, score_label, log_scale):
    n_datasets = len(datasets)
    n_methods  = len(methods)

    palette = sns.color_palette("tab10", n_colors=n_methods)

    # Width scales with number of datasets so boxes don't get crushed.
    fig_width = max(7, 2.5 * n_datasets)
    fig, ax = plt.subplots(figsize=(fig_width, 4.5), constrained_layout=True)

    sns.boxplot(data=df, x="dataset", y="score", hue="method",
                ax=ax, palette=palette, width=0.8, gap=0.05,
                flierprops={"marker": "o", "markersize": 3, "alpha": 0.5},
                showmeans=True,
                meanprops={"marker": "D", "markerfacecolor": "white",
                           "markeredgecolor": "black", "markersize": 5})

    # Width scales with number of datasets so boxes don't get crushed.
    box_width = 0.8
    dodge = box_width / n_methods

    for d_idx, dataset in enumerate(datasets):
        for m_idx, method in enumerate(methods):
            sub = df[(df["dataset"] == dataset) & (df["method"] == method)]
            if len(sub) == 0:
                continue

            n = len(sub)
            top = sub["score"].max()
            mean_val = sub["score"].mean()

            # Calculate the exact center position of this specific box
            x_center = d_idx - box_width / 2 + (m_idx + 0.5) * dodge

            # 1. Label 1: Sample size (N) placed higher up (e.g., 18 points above the max data point)
            ax.annotate(f"$M={n}$", xy=(x_center, top), xytext=(1, 35),
                        textcoords="offset points",
                        ha="left", va="center",
                        rotation=30,
                        fontsize=12, color="gray",
                        fontfamily="monospace")

            # 2. Label 2: Mean value placed directly under the N label (e.g., 4 points above the max data point)
            # Color-coded to match the box method so it's instantly recognizable
            ax.annotate(f"$\mu={mean_val:.2f}$", xy=(x_center, top), xytext=(0, 30),
                        textcoords="offset points",
                        ha="left", va="center",
                        rotation=30,
                        fontsize=12, color=palette[m_idx],
                        fontweight="bold",
                        fontfamily="monospace")

    ax.set_xlabel("")
    ax.set_ylabel(score_label, fontsize=20)
    ax.tick_params(axis="both", labelsize=18)
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)
    if n_methods > 1:
        ax.legend(fontsize=18, title=None,
                  loc="upper left", bbox_to_anchor=(1.02, 1),
                  borderaxespad=0)
    sns.despine(ax=ax)

    if log_scale:
        ax.set_yscale("log")

    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, nargs="+",
                        help="One or more 'dataset:method:path' triples. "
                             "Order determines x-axis (datasets) and legend "
                             "(methods) ordering.")
    parser.add_argument("--output", required=True,
                        help="Output stem (e.g. score_box)")
    parser.add_argument("--score-label", default="Score",
                        help="Y-axis label, e.g. '$\\chi^2$'")
    parser.add_argument("--log-scale", action="store_true",
                        help="Use log scale on y-axis")
    args = parser.parse_args()

    sns.set_theme(style="white", font_scale=1.0)

    specs = [parse_input_spec(s) for s in args.input]
    df, datasets, methods = load_long_dataframe(specs)
    stem  = str(Path(args.output).with_suffix(""))

    plot_box(df, datasets, methods,
             output_path=f"{stem}.pdf",
             score_label=args.score_label,
             log_scale=args.log_scale)

if __name__ == "__main__":
    main()
