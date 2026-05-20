#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from pathlib import Path

# TSV columns
COLS = [
    "trial", "instance", "edge_density", "V", "N",
    "ell", "d", "max_rank", "pruning_cnt",
    "num_threads_requested", "num_threads_spawned",
    "time_ms", "peak_ram_kb"
]

def load(path):
    df = pd.read_csv(path, sep="\t", header=None, names=COLS)
    df["time_s"]    = df["time_ms"] / 1000.0
    df["peak_ram_mb"] = df["peak_ram_kb"] / 1024.0
    return df

def aggregate(df, group_cols, value_col, error):
    agg = df.groupby(group_cols)[value_col].agg(
        median="median",
        vmin="min",
        vmax="max"
    ).reset_index()
    return agg

def v_to_n(df):
    """Return a dict mapping V -> median N (N is fixed per V by seed)."""
    return df.groupby("V")["N"].median().astype(int).to_dict()

def req_to_spawned(df):
    """Return a dict mapping num_threads_requested -> median num_threads_spawned."""
    return df.groupby("num_threads_requested")["num_threads_spawned"].median().astype(int).to_dict()

def make_secondary_xaxis(ax, primary_ticks, mapping, label, color="gray"):
    """Add a secondary x-axis on top with mapped tick labels."""
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(primary_ticks)
    ax2.set_xticklabels(
        [str(mapping.get(t, "")) for t in primary_ticks],
        color=color, fontsize=8
    )
    ax2.set_xlabel(label, color=color, fontsize=9)
    ax2.tick_params(axis="x", colors=color)
    return ax2

def plot_rq(df, x_col, x_label, y_col, y_label, title,
            error, secondary_map, secondary_label,
            output_path):
    """
    Generic plot: subplots per ell, lines per d.
    """
    ells = sorted(df["ell"].unique())
    n_subplots = len(ells)

    fig, axes = plt.subplots(1, n_subplots, figsize=(5 * n_subplots, 4),
                             sharey=False, constrained_layout=True)
    if n_subplots == 1:
        axes = [axes]

    palette = sns.color_palette("tab10")

    for ax, ell in zip(axes, ells):
        sub = df[df["ell"] == ell]
        ds  = sorted(sub["d"].unique())

        group_cols = [x_col, "ell", "d"]
        agg = aggregate(sub, group_cols, y_col, error)

        for i, d in enumerate(ds):
            row = agg[agg["d"] == d].sort_values(x_col)
            color = palette[i % len(palette)]
            ax.plot(row[x_col], row["median"],
                    marker="o", label=f"d={d}", color=color)
            if error:
                ax.fill_between(row[x_col], row["vmin"], row["vmax"],
                                alpha=0.2, color=color)

        ax.set_xlabel(x_label, fontsize=10)
        ax.set_ylabel(y_label, fontsize=10)
        ax.set_title(f"ell={ell}", fontsize=11)
        ax.legend(fontsize=8)
        sns.despine(ax=ax)

        # Secondary x-axis
        if secondary_map is not None:
            primary_ticks = sorted(sub[x_col].unique())
            make_secondary_xaxis(ax, primary_ticks, secondary_map,
                                 secondary_label, color="gray")

    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  required=True, help="Input TSV file")
    parser.add_argument("--rq",     required=True, choices=["rq1", "rq2", "rq3"],
                        help="Which RQ to plot")
    parser.add_argument("--output", required=True, help="Output PDF file")
    parser.add_argument("--error",  action="store_true",
                        help="Show min/max error bands")
    args = parser.parse_args()

    sns.set_theme(style="whitegrid", font_scale=1.0)

    df = load(args.input)

    if args.rq == "rq1":
        v2n = v_to_n(df)
        plot_rq(
            df,
            x_col="V",
            x_label="V (number of nodes)",
            y_col="time_s",
            y_label="Time (s)",
            title="RQ1: Impact of V and N on running time",
            error=args.error,
            secondary_map=v2n,
            secondary_label="N (total sequence length)",
            output_path=args.output
        )

    elif args.rq == "rq2":
        df["edge_density_pct"] = (df["edge_density"] * 100).round().astype(int)
        plot_rq(
            df,
            x_col="edge_density_pct",
            x_label="Edge density (%)",
            y_col="time_s",
            y_label="Time (s)",
            title="RQ2: Impact of edge density on running time",
            error=args.error,
            secondary_map=None,
            secondary_label=None,
            output_path=args.output
        )

    elif args.rq == "rq3":
        req2spawned = req_to_spawned(df)
        stem = Path(args.output).stem
        parent = Path(args.output).parent

        # Time plot
        plot_rq(
            df,
            x_col="num_threads_requested",
            x_label="Threads requested",
            y_col="time_s",
            y_label="Time (s)",
            title="RQ3: Scaling — running time",
            error=args.error,
            secondary_map=req2spawned,
            secondary_label="Threads spawned",
            output_path=str(parent / f"{stem}_time.pdf")
        )

        # Memory plot
        plot_rq(
            df,
            x_col="num_threads_requested",
            x_label="Threads requested",
            y_col="peak_ram_mb",
            y_label="Peak RAM (MB)",
            title="RQ3: Scaling — peak memory",
            error=args.error,
            secondary_map=req2spawned,
            secondary_label="Threads spawned",
            output_path=str(parent / f"{stem}_memory.pdf")
        )

if __name__ == "__main__":
    main()
