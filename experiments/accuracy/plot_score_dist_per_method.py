#!/usr/bin/env python3
"""
For one or more raw solution files (the format produced by main_algo: a header
line followed by 'X Y score' rows, sorted by score descending), plot the score
curve against rank. Each input becomes one line, with rank 1 (the
highest-scoring pair) on the far left and the last rank on the far right.
"""

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from pathlib import Path

def parse_input_spec(spec):
    """Split a 'path[:label]' spec. If no label is given, use the file stem."""
    if ":" in spec:
        path, label = spec.rsplit(":", 1)
        if not label:
            label = Path(path).stem
    else:
        path  = spec
        label = Path(spec).stem
    return path, label

def read_solution(path):
    """Read a raw main_algo solution file. Returns a list of scores in file order
    (which is already sorted by score descending). Lines whose third token can't
    be parsed as a float are skipped, which transparently handles the 'X Y x2'
    header row."""
    scores = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                score = float(parts[2])
            except ValueError:
                continue # likely the header row
            scores.append(score)
    return scores

def plot_score_rank(specs, output_path, score_label, log_scale):
    palette = sns.color_palette("tab10", n_colors=max(len(specs), 1))

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)

    max_rank = 0
    for color, (path, label) in zip(palette, specs):
        scores = read_solution(path)
        if not scores:
            print(f"Warning: no rows read from {path}")
            continue
        ranks = list(range(1, len(scores) + 1))
        ax.plot(ranks, scores, linewidth=2.0, color=color, label=label)
        max_rank = max(max_rank, len(scores))
        print(f"Loaded {len(scores)} pairs from {path} (label: {label})")

    ax.set_xlabel("Rank", fontsize=20)
    ax.set_ylabel(score_label, fontsize=20)
    ax.set_xlim(0.5, max_rank + 0.5)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.tick_params(axis="both", labelsize=18)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)
    sns.despine(ax=ax)

    if len(specs) > 1:
        ax.legend(fontsize=18)

    if log_scale:
        ax.set_yscale("log")

    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",       required=True, nargs="+",
                        help="One or more 'path[:label]' inputs. The label "
                             "becomes the legend entry; if omitted, the file "
                             "stem is used.")
    parser.add_argument("--output",      required=True,
                        help="Output stem (e.g. rank)")
    parser.add_argument("--score-label", default="Score",
                        help="Y-axis label, e.g. '$\\chi^2$'")
    parser.add_argument("--log-scale",   action="store_true",
                        help="Use log scale on y-axis")
    args = parser.parse_args()

    sns.set_theme(style="white", font_scale=1.0)

    specs = [parse_input_spec(s) for s in args.input]
    stem  = str(Path(args.output).with_suffix(""))

    plot_score_rank(
        specs       = specs,
        output_path = f"{stem}.pdf",
        score_label = args.score_label,
        log_scale   = args.log_scale,
    )

if __name__ == "__main__":
    main()
