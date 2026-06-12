#!/usr/bin/env python3
"""
Grouped box plot of per-(dataset, method) score distributions from main_algo
solution files (typically after deduplication via cmm_dedup).

X-axis: datasets. Within each dataset, one box per method (color-coded).

When --baseline METHOD is given, each non-baseline method's scores are paired
by rank (file order) with the baseline's, per dataset, truncated to the
shorter of the two, and converted to (method_score / baseline_score * 100).
The baseline method drops out of the plot; the resulting boxes show "% of
baseline score" at each rank position.

Input is one or more 'dataset:method:path' triples; the dataset and method
names cannot contain colons.
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
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
    sorting alphabetically). Row order within each (dataset, method) group is
    rank order (file order); the normalization step relies on this."""
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


def normalize_against_baseline(df, baseline_method):
    """For each dataset, pair each non-baseline method's scores with the
    baseline's by rank (file order), truncate both to min(len_a, len_b), and
    replace 'score' with (method_score / baseline_score * 100). The baseline
    method is dropped from the output. Pairs where the baseline score is 0
    (undefined normalization) are dropped with a note."""
    out_rows = []

    for dataset in df["dataset"].cat.categories:
        ds_df = df[df["dataset"] == dataset]
        baseline_df = ds_df[ds_df["method"] == baseline_method]

        if baseline_df.empty:
            print(f"Warning: dataset '{dataset}' has no entries for baseline "
                  f"method '{baseline_method}', skipping")
            continue

        baseline_scores = baseline_df["score"].to_numpy()

        for method in df["method"].cat.categories:
            if method == baseline_method:
                continue
            method_df = ds_df[ds_df["method"] == method]
            if method_df.empty:
                continue

            method_scores = method_df["score"].to_numpy()
            n = min(len(baseline_scores), len(method_scores))

            if len(baseline_scores) != len(method_scores):
                print(f"Warning: ({dataset}, {method}): {len(method_scores)} "
                      f"entries vs baseline's {len(baseline_scores)}, "
                      f"truncating both to {n}")

            a = baseline_scores[:n]
            b = method_scores[:n]

            with np.errstate(divide="ignore", invalid="ignore"):
                ratios = (b / a) * 100.0
            ratios = np.where(a > 0, ratios, np.nan)

            dropped = int(np.isnan(ratios).sum())
            if dropped:
                print(f"Note: dropped {dropped} pair(s) where baseline score "
                      f"== 0 for ({dataset}, {method})")

            for r in ratios:
                if np.isnan(r):
                    continue
                out_rows.append({"dataset": dataset, "method": method,
                                 "score": float(r)})

    new_df = pd.DataFrame(out_rows)
    if new_df.empty:
        return new_df, [], []

    new_datasets = [d for d in df["dataset"].cat.categories
                    if d in set(new_df["dataset"].unique())]
    new_methods  = [m for m in df["method"].cat.categories
                    if m != baseline_method
                    and m in set(new_df["method"].unique())]
    new_df["dataset"] = pd.Categorical(new_df["dataset"],
                                       categories=new_datasets, ordered=True)
    new_df["method"]  = pd.Categorical(new_df["method"],
                                       categories=new_methods, ordered=True)
    return new_df, new_datasets, new_methods


def plot_box(df, datasets, methods, output_path, score_label,
             reference_line=None):
    n_datasets = len(datasets)
    n_methods  = len(methods)

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
    if n_methods > 1:
        ax.legend(fontsize=14, loc="best")
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
                        help="One or more 'dataset:method:path' triples. "
                             "Order determines x-axis (datasets) and legend "
                             "(methods) ordering.")
    parser.add_argument("--output", required=True,
                        help="Output stem (e.g. score_box)")
    parser.add_argument("--score-label", default="Score",
                        help="Y-axis label, e.g. 'Score / exact (\\%%)' when "
                             "normalizing.")
    parser.add_argument("--baseline", default=None,
                        help="Method name (matching one in --input) to "
                             "normalize against. When set, each other "
                             "method's scores are paired by rank with the "
                             "baseline's and shown as a percentage. The "
                             "baseline drops out of the plot.")
    args = parser.parse_args()

    sns.set_theme(style="white", font_scale=1.0)

    specs = [parse_input_spec(s) for s in args.input]
    df, datasets, methods = load_long_dataframe(specs)

    reference_line = None
    if args.baseline is not None:
        if args.baseline not in methods:
            raise SystemExit(
                f"Error: --baseline '{args.baseline}' not found among input "
                f"methods: {methods}")
        df, datasets, methods = normalize_against_baseline(df, args.baseline)
        if df.empty:
            raise SystemExit(
                "Error: no non-baseline data remained after normalization")
        reference_line = 100.0

    stem = str(Path(args.output).with_suffix(""))
    plot_box(df, datasets, methods,
             output_path=f"{stem}.pdf",
             score_label=args.score_label,
             reference_line=reference_line)


if __name__ == "__main__":
    main()
