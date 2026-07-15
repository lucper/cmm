#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# TSV columns
COLS = [
    "trial", "instance", "edge_density", "V", "N",
    "ell", "d", "max_rank", "min_rank", "avg_rank", "pruning_cnt",
    "num_threads_requested", "num_threads_spawned",
    "time_ms", "peak_ram_kb"
]

ELLS = [5, 8]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

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

def d_groups(ell):
    """Split d values for a given ell into lower and upper halves."""
    half = ell // 2
    return [
        (f"d0_{half-1}",     list(range(0, half))),
        (f"d{half}_{ell-1}", list(range(half, ell))),
    ]

def ell_d_pair(s):
    """Argparse type: parse 'ell:d' into an (int, int) tuple."""
    try:
        ell_str, d_str = s.split(":")
        return (int(ell_str), int(d_str))
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid value: {s!r}. Expected format ell:d, e.g. 8:3"
        )

def save_legend(handles, labels, output_path, max_cols=4):
    """Save a standalone legend PDF, fixed canvas size for uniform scaling.
    Uses min(len(labels), max_cols) columns; the legend itself centers in
    the canvas so fewer entries cluster in the middle with whitespace.
    """
    fig = plt.figure(figsize=(11, 0.6))  # wide enough for 4 long entries
    fig.legend(handles, labels, loc="center", ncol=min(len(labels), max_cols),
               fontsize=18, frameon=False, handlelength=2.5,
               columnspacing=1.5)
    fig.savefig(output_path, format="pdf")  # no bbox_inches="tight"
    plt.close(fig)
    print(f"Saved: {output_path}")

def plot_one(df_ell, x_col, x_label, y_col, y_label, title,
             error, secondary_map, secondary_label, output_path,
             ell=None, annotate_pairs=None, log_scale=False,
             show_legend=True):
    """Plot one (ell, x, y) view with lines per d.
    annotate_pairs: set of (ell, d) tuples to annotate with max_rank, or None.
    """
    ds      = sorted(df_ell["d"].unique())
    palette = sns.color_palette("tab10", n_colors=len(ds))

    x_vals      = sorted(df_ell[x_col].unique())
    x_to_pos    = {v: i for i, v in enumerate(x_vals)}
    x_positions = list(range(len(x_vals)))

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)

    for i, (color, d) in enumerate(zip(palette, ds)):
        marker = MARKERS[i % len(MARKERS)]
        sub = df_ell[df_ell["d"] == d].groupby(x_col).agg(
            median=(y_col, "median"),
            vmin=(y_col, "min"),
            vmax=(y_col, "max"),
            max_rank=("max_rank", "first")
        ).reset_index().sort_values(x_col)

        xs = [x_to_pos[v] for v in sub[x_col]]

        label = f"($\\ell={ell}$, $d$={d})"
        ax.plot(xs, sub["median"],
                marker=marker, linewidth=2.0, markersize=7,
                label=label, color=color)
        if error:
            ax.fill_between(xs, sub["vmin"], sub["vmax"],
                            alpha=0.2, color=color)

        if annotate_pairs is not None and ell is not None and (ell, d) in annotate_pairs:
            for x, y, r in zip(xs, sub["median"], sub["max_rank"]):
                ax.annotate(
                    f"{int(r):,}",
                    xy=(x, y),
                    xytext=(0, 6),
                    textcoords="offset points",
                    ha="center", va="bottom",
                    fontsize=14, color="gray",
                    fontfamily="monospace"
                )

    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(v) for v in x_vals], fontsize=18)
    ax.set_xlim(-0.5, len(x_vals) - 0.5)

    if log_scale:
        ax.set_yscale("log")
        y_label = y_label.replace(")", ", log scale)")

    ax.set_xlabel(x_label, fontsize=20)
    ax.set_ylabel(y_label, fontsize=20)
    ax.set_title(title, fontsize=12)
    ax.tick_params(axis="y", labelsize=18)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)
    sns.despine(ax=ax)

    # Capture handles before twin axis so the legend only carries primary lines
    handles, labels = ax.get_legend_handles_labels()
    if show_legend:
        ax.legend(handles, labels, fontsize=18)

    if secondary_map is not None:
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(x_positions)
        ax2.set_xticklabels([str(secondary_map.get(v, "")) for v in x_vals],
                            fontsize=18, rotation=30, ha="left")
        ax2.set_xlabel(secondary_label, fontsize=20)
        ax2.tick_params(axis="x", labelsize=18)

    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")

    return handles, labels

def run_rq3(df, stem, error, annotate_pairs, log_scale):
    v2n = v_to_n(df)
    for ell in ELLS:
        for group_name, ds in d_groups(ell):
            sub  = df[(df["ell"] == ell) & (df["d"].isin(ds))]
            base = f"{stem}_ell{ell}_{group_name}"

            h, l = plot_one(
                df_ell=sub,
                x_col="V",
                x_label="$V$",
                y_col="time_s",
                y_label="Time (s)",
                title="",
                error=error,
                secondary_map=v2n,
                secondary_label="$N$",
                output_path=f"{base}_time.pdf",
                ell=ell,
                annotate_pairs=annotate_pairs,
                log_scale=log_scale,
                show_legend=False
            )
            plot_one(
                df_ell=sub,
                x_col="V",
                x_label="$V$",
                y_col="peak_ram_mb",
                y_label="Peak RAM (MB)",
                title="",
                error=error,
                secondary_map=v2n,
                secondary_label="$N$",
                output_path=f"{base}_memory.pdf",
                ell=ell,
                annotate_pairs=None,
                log_scale=False,
                show_legend=False
            )
            save_legend(h, l, f"{base}_legend.pdf")

def run_rq4(df, stem, error, annotate_pairs, log_scale):
    df["edge_density_pct"] = (df["edge_density"] * 100).round().astype(int)
    for ell in ELLS:
        for group_name, ds in d_groups(ell):
            sub  = df[(df["ell"] == ell) & (df["d"].isin(ds))]
            base = f"{stem}_ell{ell}_{group_name}"

            h, l = plot_one(
                df_ell=sub,
                x_col="edge_density_pct",
                x_label="Edge density (%)",
                y_col="time_s",
                y_label="Time (s)",
                title="",
                error=error,
                secondary_map=None,
                secondary_label=None,
                output_path=f"{base}_time.pdf",
                ell=ell,
                annotate_pairs=annotate_pairs,
                log_scale=log_scale,
                show_legend=False
            )
            plot_one(
                df_ell=sub,
                x_col="edge_density_pct",
                x_label="Edge density (%)",
                y_col="peak_ram_mb",
                y_label="Peak RAM (MB)",
                title="",
                error=error,
                secondary_map=None,
                secondary_label=None,
                output_path=f"{base}_memory.pdf",
                ell=ell,
                annotate_pairs=None,
                log_scale=False,
                show_legend=False
            )
            save_legend(h, l, f"{base}_legend.pdf")

def run_rq5(df, stem, error, annotate_pairs, log_scale, ell=8):
    req2spawned = req_to_spawned(df)

    h, l = plot_one(
        df_ell=df[df["ell"] == ell],
        x_col="num_threads_requested",
        x_label="Number of threads",
        y_col="time_s",
        y_label="Time (s)",
        title="",
        error=error,
        secondary_map=None,
        secondary_label=None,
        output_path=f"{stem}_time.pdf",
        ell=ell,
        annotate_pairs=annotate_pairs,
        log_scale=log_scale,
        show_legend=False
    )

    plot_one(
        df_ell=df[df["ell"] == ell],
        x_col="num_threads_requested",
        x_label="Number of threads",
        y_col="peak_ram_mb",
        y_label="Peak RAM (MB)",
        title="",
        error=error,
        secondary_map=None,
        secondary_label=None,
        output_path=f"{stem}_memory.pdf",
        ell=ell,
        annotate_pairs=None,
        log_scale=False,
        show_legend=False
    )
    save_legend(h, l, f"{stem}_legend.pdf")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",    required=True)
    parser.add_argument("--rq",       required=True, choices=["rq3", "rq4", "rq5"])
    parser.add_argument("--output",   required=True,
                        help="Output stem (e.g. rq3)")
    parser.add_argument("--error",    action="store_true",
                        help="Show min/max error bands")
    parser.add_argument("--max-rank", nargs="+", type=ell_d_pair, default=None,
                        metavar="ell:d",
                        help="Annotate points with max_rank for given (ell,d) pairs (e.g. --max-rank 5:0 8:3)")
    parser.add_argument("--log-scale", action="store_true",
                        help="Use log scale on y-axis")
    parser.add_argument("--rq5-ell", type=int, default=5,
                        help="Motif length used for the RQ5 thread-scaling run (default: 5)")
    args = parser.parse_args()

    sns.set_theme(style="white", font_scale=1.0)

    df   = load(args.input)
    stem = str(Path(args.output).with_suffix(""))  # strip extension if given

    annotate_pairs = set(args.max_rank) if args.max_rank is not None else None

    if args.rq == "rq3":
        run_rq3(df, stem, args.error, annotate_pairs, args.log_scale)
    elif args.rq == "rq4":
        run_rq4(df, stem, args.error, annotate_pairs, args.log_scale)
    elif args.rq == "rq5":
        run_rq5(df, stem, args.error, annotate_pairs, args.log_scale, args.rq5_ell)

if __name__ == "__main__":
    main()
