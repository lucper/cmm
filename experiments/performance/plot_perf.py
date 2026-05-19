#!/usr/bin/env python3
"""
Plot CMM performance results.

One figure per (edge_density, ell).  Each figure has:
  - y-axis : time in milliseconds (log scale)
  - x-axis : instances sorted by V, with dual tick labels showing V and
             total_label_len beneath it
  - legend : one line per d value

Usage:
    python3 plot_perf.py results.tsv [--out-dir plots/]
"""

import argparse
import os
from itertools import cycle

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd

# One visually distinct colour + marker per d value (up to d=5)
COLORS  = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
MARKERS = ["o", "s", "^", "D", "v", "P"]

def load(tsv_path: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")
    # Ensure numeric types
    df["edge_density"]    = df["edge_density"].astype(float)
    df["V"]               = df["V"].astype(int)
    df["total_label_len"] = df["total_label_len"].astype(int)
    df["ell"]             = df["ell"].astype(int)
    df["d"]               = df["d"].astype(int)
    df["time_ms"]         = df["time_ms"].astype(float)
    return df

def make_plots(df: pd.DataFrame, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)

    for edge_density, df_ed in df.groupby("edge_density"):
        for (ell, df_ell) in df_ed.groupby("ell"):
            fig, ax = plt.subplots(figsize=(10, 5))

            # x positions: one per distinct (V, total_label_len) pair, sorted by V
            instance_keys = (
                df_ell[["V", "total_label_len"]]
                .drop_duplicates()
                .sort_values("V")
                .reset_index(drop=True)
            )
            x_positions = range(len(instance_keys))
            v_to_x = {row.V: i for i, row in instance_keys.iterrows()}

            for d, df_d in df_ell.groupby("d"):
                df_d = df_d.sort_values("V")
                xs = [v_to_x[v] for v in df_d["V"]]
                color  = COLORS[d % len(COLORS)]
                marker = MARKERS[d % len(MARKERS)]
                ax.plot(xs, df_d["time_ms"] / 1000,
                        color=color, marker=marker,
                        linewidth=1.8, markersize=6,
                        label=f"($\\ell={ell}$, $d$={d})")

            # Primary x ticks: V values
            ax.set_xticks(list(x_positions))
            ax.set_xticklabels(
                [str(row.V) for _, row in instance_keys.iterrows()],
                fontsize=9,
            )

            # Secondary x ticks beneath: total_label_len
            ax2 = ax.twiny()
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(list(x_positions))
            ax2.set_xticklabels(
                [f"{row.total_label_len:,}" for _, row in instance_keys.iterrows()],
                fontsize=7.5, color="#555555",
            )
            ax2.xaxis.set_ticks_position("bottom")
            ax2.xaxis.set_label_position("bottom")
            ax2.spines["bottom"].set_position(("outward", 25))
            ax2.spines["bottom"].set_color(("#555555"))
            ax2.tick_params(axis="x", colors="#555555")
            ax2.set_xlabel("")

            ax.set_xlabel("")
            fig.canvas.draw()
            ax.annotate("$|V|$\t",
                xy=(0.5, 0), xycoords="axes fraction",
                xytext=(0, -48), textcoords="offset points",
                fontsize=10, color="black",
                ha="right", va="top", annotation_clip=False,
            )
            ax.annotate("$N$",
                xy=(0.5, 0), xycoords="axes fraction",
                xytext=(0, -48), textcoords="offset points",
                fontsize=10, color="#555555",
                ha="left", va="top", annotation_clip=False,
            )

            ax.set_ylabel("Time (s)", fontsize=10)
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(
                lambda x, _: f"{int(x):,}" if x >= 1 else f"{x:.2g}"
            ))
            ax.set_title(
                f"Edge density = {edge_density:.2f}",
                fontsize=11,
            )
            ax.legend(title="combination", fontsize=9, title_fontsize=9,
                      loc="upper left", framealpha=0.85)
            ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)

            fname = f"perf_ed{int(round(edge_density * 100)):02d}_ell{ell}.png"
            fpath = os.path.join(out_dir, fname)
            fig.tight_layout()
            fig.savefig(fpath, bbox_inches="tight", dpi=300)
            plt.close(fig)
            print(f"Saved: {fpath}")

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tsv",     help="TSV file produced by cmm_perf")
    parser.add_argument("--out-dir", default="plots", help="Output directory for PDFs")
    args = parser.parse_args()

    df = load(args.tsv)
    make_plots(df, args.out_dir)


if __name__ == "__main__":
    main()
