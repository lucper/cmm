#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# TSV columns
COLS = [
    "trial", "instance", "edge_density", "V", "N",
    "ell", "d", "max_rank", "pruning_cnt",
    "num_threads_requested", "num_threads_spawned",
    "time_ms", "peak_ram_kb"
]

ELLS = [5, 8]

def load(path):
    df = pd.read_csv(path, sep="\t", header=None, names=COLS)
    df = df[df["ell"].isin(ELLS)]
    df["time_s"]      = df["time_ms"] / 1000.0
    df["peak_ram_mb"] = df["peak_ram_kb"] / 1024.0
    return df

def v_to_n(df):
    return df.drop_duplicates("V").set_index("V")["N"].to_dict()

def req_to_spawned(df):
    return df.drop_duplicates("num_threads_requested") \
             .set_index("num_threads_requested")["num_threads_spawned"].to_dict()

def plot_one(df_ell, x_col, x_label, y_col, y_label, title,
             error, secondary_map, secondary_label, output_path,
             ell=None, show_max_rank=False):
    ds      = sorted(df_ell["d"].unique())
    palette = sns.color_palette("tab10", n_colors=len(ds))

    x_vals      = sorted(df_ell[x_col].unique())
    x_to_pos    = {v: i for i, v in enumerate(x_vals)}
    x_positions = list(range(len(x_vals)))

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)

    for color, d in zip(palette, ds):
        sub = df_ell[df_ell["d"] == d].groupby(x_col).agg(
            median=(y_col, "median"),
            vmin=(y_col, "min"),
            vmax=(y_col, "max"),
            max_rank=("max_rank", "first")
        ).reset_index().sort_values(x_col)

        xs = [x_to_pos[v] for v in sub[x_col]]

        label = f"($\\ell$={ell}, $d$={d})" if ell else f"$d$={d}"
        ax.plot(xs, sub["median"],
                marker="o", linewidth=2.0, markersize=7,
                label=label, color=color)
        if error:
            ax.fill_between(xs, sub["vmin"], sub["vmax"],
                            alpha=0.2, color=color)

        if show_max_rank:
            for x, y, r in zip(xs, sub["median"], sub["max_rank"]):
                ax.annotate(
                    f"{int(r):,}",
                    xy=(x, y),
                    xytext=(0, 6),
                    textcoords="offset points",
                    ha="center", va="bottom",
                    fontsize=7, color="gray",
                    fontfamily="monospace"
                )

    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(v) for v in x_vals], fontsize=10)
    ax.set_xlim(-0.5, len(x_vals) - 0.5)

    ax.set_xlabel(x_label, fontsize=11)
    ax.set_ylabel(y_label, fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.tick_params(axis="y", labelsize=10)
    ax.legend(fontsize=9)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)
    sns.despine(ax=ax)

    if secondary_map is not None:
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(x_positions)
        ax2.set_xticklabels(
            [str(secondary_map.get(v, "")) for v in x_vals],
            color="gray", fontsize=8
        )
        ax2.set_xlabel(secondary_label, color="gray", fontsize=9)
        ax2.tick_params(axis="x", colors="gray")

    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")

def run_rq1(df, stem, error, max_rank):
    v2n = v_to_n(df)
    for ell in ELLS:
        plot_one(
            df_ell=df[df["ell"] == ell],
            x_col="V",
            x_label="$V$",
            y_col="time_s",
            y_label="Time (s)",
            title="",
            error=error,
            secondary_map=v2n,
            secondary_label="$N$",
            output_path=f"{stem}_ell{ell}.pdf",
            ell=ell,
            show_max_rank=max_rank
        )

def run_rq2(df, stem, error, max_rank):
    df["edge_density_pct"] = (df["edge_density"] * 100).round().astype(int)
    for ell in ELLS:
        plot_one(
            df_ell=df[df["ell"] == ell],
            x_col="edge_density_pct",
            x_label="Edge density (%)",
            y_col="time_s",
            y_label="Time (s)",
            title="",
            error=error,
            secondary_map=None,
            secondary_label=None,
            output_path=f"{stem}_ell{ell}.pdf",
            ell=ell,
            show_max_rank=max_rank
        )

def run_rq3(df, stem, error, max_rank):
    req2spawned = req_to_spawned(df)

    # Time plot
    plot_one(
        df_ell=df,
        x_col="num_threads_requested",
        x_label="Threads requested",
        y_col="time_s",
        y_label="Time (s)",
        title="",
        error=error,
        secondary_map=req2spawned,
        secondary_label="Threads spawned",
        output_path=f"{stem}_time.pdf",
        show_max_rank=max_rank
    )

    # Memory plot
    plot_one(
        df_ell=df,
        x_col="num_threads_requested",
        x_label="Threads requested",
        y_col="peak_ram_mb",
        y_label="Peak RAM (MB)",
        title="",
        error=error,
        secondary_map=req2spawned,
        secondary_label="Threads spawned",
        output_path=f"{stem}_memory.pdf",
        show_max_rank=False  # max_rank not meaningful for memory plot
    )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",    required=True)
    parser.add_argument("--rq",       required=True, choices=["rq1", "rq2", "rq3"])
    parser.add_argument("--output",   required=True,
                        help="Output stem (e.g. rq1)")
    parser.add_argument("--error",    action="store_true",
                        help="Show min/max error bands")
    parser.add_argument("--max-rank", action="store_true",
                        help="Annotate each point with its max_rank value")
    args = parser.parse_args()

    sns.set_theme(style="white", font_scale=1.0)

    df   = load(args.input)
    stem = str(Path(args.output).with_suffix(""))  # strip extension if given

    if args.rq == "rq1":
        run_rq1(df, stem, args.error, args.max_rank)
    elif args.rq == "rq2":
        run_rq2(df, stem, args.error, args.max_rank)
    elif args.rq == "rq3":
        run_rq3(df, stem, args.error, args.max_rank)

if __name__ == "__main__":
    main()
