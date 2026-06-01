#!/usr/bin/env python3
"""
For each row in the input TSV -- a motif pair from one solution and the closest
counterpart in the other -- plot (distance, score). Points in the upper-right
corner (high score, high distance) are high-scoring pairs whose best match in
the other solution is far away.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

EXPECTED_COLS = {
    "rank", "X", "Y", "score",
    "best_X", "best_Y", "best_score",
    "similarity", "distance",
}

def load(path):
    df = pd.read_csv(path, sep="\t")
    missing = EXPECTED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"TSV missing expected columns: {sorted(missing)}")
    return df

def plot_score_distance(df, output_path, score_label, log_scale, threshold):
    palette = sns.color_palette("tab10", n_colors=10)

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)

    if threshold is not None:
        below = df[df["distance"] <  threshold]
        above = df[df["distance"] >= threshold]
        n_above = len(above)
        n_total = len(df)
        pct = 100.0 * n_above / n_total if n_total else 0.0

        ax.scatter(below["distance"], below["score"],
                   s=30, alpha=0.5, color=palette[0], edgecolor="none",
                   label=f"$d$ < {threshold}")
        ax.scatter(above["distance"], above["score"],
                   s=30, alpha=0.7, color=palette[3], edgecolor="none",
                   label=f"$d$ $\\geq$ {threshold} "
                         f"({n_above}/{n_total}, {pct:.1f}%)")
        ax.axvline(x=threshold, color="gray",
                   linestyle="--", linewidth=1.0, alpha=0.7)
        ax.legend(fontsize=14, loc="best")
    else:
        ax.scatter(df["distance"], df["score"],
                   s=30, alpha=0.6, color=palette[0], edgecolor="none")

    ax.set_xlabel("Distance ($1-J_h$)", fontsize=20)
    ax.set_ylabel(score_label, fontsize=20)
    ax.set_xlim(-0.02, 1.02)
    ax.tick_params(axis="both", labelsize=18)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)
    sns.despine(ax=ax)

    if log_scale:
        ax.set_yscale("log")

    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",       required=True,
                        help="cmm_eval best-match TSV (A->B or B->A)")
    parser.add_argument("--output",      required=True,
                        help="Output stem (e.g. a2b_scatter)")
    parser.add_argument("--score-label", default="Score",
                        help="Y-axis label, e.g. '$\\chi^2$'")
    parser.add_argument("--threshold",   type=float, default=None,
                        help="Distance threshold d* in [0,1] to mark with a vertical line")
    parser.add_argument("--log-scale",   action="store_true",
                        help="Use log scale on y-axis")
    args = parser.parse_args()

    sns.set_theme(style="white", font_scale=1.0)

    df = load(args.input)
    stem = str(Path(args.output).with_suffix(""))

    plot_score_distance(
        df=df,
        output_path=f"{stem}.pdf",
        score_label=args.score_label,
        log_scale=args.log_scale,
        threshold=args.threshold,
    )

if __name__ == "__main__":
    main()
