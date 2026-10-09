#!/usr/bin/env python3

import argparse
import numpy as np
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

MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]

# Appearance of the fitted polynomial reference curve.
FIT_STYLE = dict(
    linestyle=(0, (6, 4)),  # long dash; gap widened to stay legible when thick
    linewidth=2.2,
    color="0.25",           # dark gray, stays neutral against the tab10 palette
    alpha=0.9,
    zorder=4,               # above the data lines (~2) and their markers (~3)
)

# Appearance of the exponent label for the fitted line.
FIT_LABEL_XY = (0.97, 0.04)
FIT_LABEL_STYLE = dict(
    fontsize=14,
    color="0.25",
    fontfamily="monospace",
    ha="right",
    va="bottom",
)

# Appearance of the per-point speedup labels (RQ5).
SPEEDUP_LABEL_STYLE = dict(
    fontsize=13,
    color="gray",
    fontfamily="monospace",
    ha="center",
    va="bottom",
)

def fit_polynomial(x, y, degree=2, relative=True):
    """Fit y = c[0]*x^degree + ... + c[degree] (coeffs highest power first).

    relative=True weights each point by 1/y, so points spanning several
    decades contribute evenly (minimises relative rather than absolute
    error). Without it, plain polyfit is dominated by the largest y values
    and fits small inputs poorly.

    Returns the coefficient array, or None if the data cannot support a fit
    (fewer than degree+1 distinct x values, or non-positive y under relative
    weighting).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    ok = np.isfinite(x) & np.isfinite(y)
    if relative:
        ok &= (y > 0)
    x, y = x[ok], y[ok]
    if len(x) < degree + 1 or len(np.unique(x)) < degree + 1:
        return None

    w = 1.0 / y if relative else None
    return np.polyfit(x, y, degree, w=w)

def _fmt_coeff(v):
    """Format a coefficient magnitude legibly.

    Uses plain decimals for moderate magnitudes and scientific notation
    (mantissa x 10^exp) only when the value is very small or very large. The
    sign is handled by the caller; this returns the magnitude string only.
    """
    a = abs(v)
    if a == 0:
        return "0"
    if 1e-3 <= a < 1e5:
        # plain decimal, trimming trailing zeros; more decimals for small values
        if a >= 100:
            s = f"{a:.0f}"
        elif a >= 1:
            s = f"{a:.2f}"
        else:
            s = f"{a:.4f}"
        return s.rstrip("0").rstrip(".") if "." in s else s
    # scientific, rendered as m\times10^{e} for mathtext
    exp = int(np.floor(np.log10(a)))
    mant = a / 10**exp
    mant_s = f"{mant:.2f}".rstrip("0").rstrip(".")
    return f"{mant_s}\\times10^{{{exp}}}"

def _format_poly(coeffs, var="V"):
    """Render a coefficient array as a legible 'aV^2 + bV + c' mathtext string.

    Signs are pulled out so terms read '... - 1.81V + 120' rather than
    '... + -1.81V'. Magnitudes use plain decimals when readable and fall back
    to m x 10^e only for extreme values. `var` is the axis variable name.
    """
    deg = len(coeffs) - 1
    out = ""
    for i, c in enumerate(coeffs):
        p = deg - i
        if c == 0:
            continue
        sign = "-" if c < 0 else "+"
        mag = _fmt_coeff(c)
        if p == 0:
            term = mag
        elif p == 1:
            term = f"{mag}{var}"
        else:
            term = f"{mag}{var}^{{{p}}}"
        if out == "":
            out = (f"-{term}" if sign == "-" else term)   # leading sign only if negative
        else:
            out += f" {sign} {term}"
    return out if out else "0"


def load(path, mem_per_thread=False):
    df = pd.read_csv(path, sep="\t", header=None, names=COLS)
    df["time_s"]      = df["time_ms"] / 1000.0
    df["peak_ram_mb"] = df["peak_ram_kb"] / 1024.0
    # Peak RSS of the whole process divided by the threads actually spawned,
    # which can be fewer than requested.
    if mem_per_thread:
        df["peak_ram_mb"] /= df["num_threads_spawned"]
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
    """Save a standalone legend PDF on a fixed-size canvas.

    The fixed canvas keeps legends across plots at a uniform scale. Uses
    min(len(labels), max_cols) columns; the legend centers in the canvas, so
    fewer entries cluster in the middle with surrounding whitespace.
    """
    if not labels:
        return
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
             show_legend=True, speedup=False,
             fit_poly_pairs=None, poly_degree=2, poly_label=True):
    """Plot one (ell, x, y) view with one line per d value.

    Args:
        annotate_pairs: Set of (ell, d) tuples to annotate with max_rank,
            or None.
        fit_poly_pairs: Set of (ell, d) tuples to overlay a polynomial fit
            (a*x^deg + ... + c) on, or None. Fit with 1/y weighting so it is
            not dominated by large-y points; evaluated directly at each
            point's true x (no log-axis interpolation), then plotted at that
            point's position.
        poly_degree: Degree of the polynomial fit (default 2).
        poly_label: If True, annotate the polynomial curve with its equation
            in the bottom-right corner. The curve is drawn regardless; this
            only controls the text label.
        speedup: If True, label each point with its speedup relative to the
            smallest x value in its own series (so the leftmost point reads
            1.0x). Only meaningful when x is a thread count and y is a time.
    """
    ds      = sorted(df_ell["d"].unique())
    palette = sns.color_palette("tab10", n_colors=len(ds))

    x_vals      = sorted(df_ell[x_col].unique())
    x_to_pos    = {v: i for i, v in enumerate(x_vals)}
    x_positions = list(range(len(x_vals)))

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)

    # Track the extent of the plotted data (incl. error bands) so an
    # extrapolated fit line cannot rescale the y-axis.
    y_data_lim = None

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

        seen_lo = float(min(sub["vmin"].min(), sub["median"].min()) if error
                        else sub["median"].min())
        seen_hi = float(max(sub["vmax"].max(), sub["median"].max()) if error
                        else sub["median"].max())
        y_data_lim = (seen_lo, seen_hi) if y_data_lim is None else (
            min(y_data_lim[0], seen_lo), max(y_data_lim[1], seen_hi))

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

        if speedup:
            med = np.asarray(sub["median"], dtype=float)
            # Baseline is this series' own smallest x (sub is sorted by x_col).
            base = med[0]
            if base > 0 and np.all(med > 0):
                for x, y, sp in zip(xs, med, base / med):
                    ax.annotate(
                        f"{sp:.1f}x",
                        xy=(x, y),
                        xytext=(0, 8),
                        textcoords="offset points",
                        **SPEEDUP_LABEL_STYLE
                    )
            else:
                print(f"  speedup skipped (ell={ell}, d={d}): non-positive times")

        if fit_poly_pairs is not None and ell is not None and (ell, d) in fit_poly_pairs:
            coeffs = fit_polynomial(sub[x_col], sub["median"], degree=poly_degree)
            if coeffs is None:
                print(f"  poly-fit [{y_col}] skipped (ell={ell}, d={d}): "
                      f"needs >={poly_degree + 1} distinct positive points")
            else:
                # Evaluate at each point's OWN true x, placed at its OWN
                # position -- no interpolation across the (possibly non-
                # geometric) x spacing, so the curve tracks the fit exactly.
                pos_known = np.array([x_to_pos[v] for v in sub[x_col]], dtype=float)
                x_known   = np.asarray(sub[x_col], dtype=float)
                y_poly    = np.polyval(coeffs, x_known)

                # Extend to the borders by evaluating the polynomial at the
                # x implied by the end segments' spacing.
                pos_edges = np.array([-0.5, len(x_vals) - 0.5], dtype=float)
                dx_lo = (x_known[1]  - x_known[0])  / (pos_known[1]  - pos_known[0])
                dx_hi = (x_known[-1] - x_known[-2]) / (pos_known[-1] - pos_known[-2])
                x_edges = np.array([
                    x_known[0]  + dx_lo * (pos_edges[0]  - pos_known[0]),
                    x_known[-1] + dx_hi * (pos_edges[1]  - pos_known[-1]),
                ])
                pos_line = np.concatenate([[pos_edges[0]], pos_known, [pos_edges[1]]])
                x_line   = np.concatenate([[x_edges[0]],   x_known,   [x_edges[1]]])
                y_line   = np.polyval(coeffs, x_line)

                var_sym = {"V": "x", "edge_density_pct": "\\rho",
                           "num_threads_requested": "t"}.get(x_col, "x")
                ax.plot(pos_line, y_line, label="_nolegend_", **FIT_STYLE)
                if poly_label:
                    ax.annotate(
                        f"${_format_poly(coeffs, var=var_sym)}$",
                        xy=FIT_LABEL_XY,
                        xycoords="axes fraction",
                        zorder=5,
                        **FIT_LABEL_STYLE
                    )
                # Keep the extrapolated curve from rescaling the y-axis.
                lo = float(min(y_poly.min(), sub["median"].min()))
                hi = float(max(y_poly.max(), sub["median"].max()))
                y_data_lim = (lo, hi) if y_data_lim is None else (
                    min(y_data_lim[0], lo), max(y_data_lim[1], hi))
                coeff_str = ", ".join(f"{c:.4e}" for c in coeffs)
                print(f"  poly-fit [{y_col}] (ell={ell}, d={d}, deg={poly_degree}): "
                      f"[{coeff_str}]")

    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(v) for v in x_vals], fontsize=18)
    ax.set_xlim(-0.5, len(x_vals) - 0.5)

    if log_scale:
        ax.set_yscale("log")

    # Fits are extrapolated to the plot borders, which would otherwise stretch
    # the y-axis and squash the data. Clip the view back to the data extent.
    if fit_poly_pairs and y_data_lim is not None:
        lo_y, hi_y = y_data_lim
        if log_scale:
            pad = (hi_y / lo_y) ** 0.05
            ax.set_ylim(lo_y / pad, hi_y * pad)
        else:
            pad = 0.05 * (hi_y - lo_y)
            ax.set_ylim(lo_y - pad, hi_y + pad)

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

def run_rq3(df, stem, error, annotate_pairs, log_scale,
            log_scale_mem=False,
            fit_poly_pairs=None, fit_poly_pairs_mem=None,
            poly_degree_time=2, poly_degree_mem=2,
            poly_label=True):
    v2n = v_to_n(df)
    for ell in sorted(df["ell"].unique()):
        for group_name, ds in d_groups(ell):
            sub  = df[(df["ell"] == ell) & (df["d"].isin(ds))]
            if sub.empty:
                continue
            base = f"{stem}_ell{ell}_{group_name}"

            h, l = plot_one(
                df_ell=sub,
                x_col="V",
                x_label="$|V|$",
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
                show_legend=False,
                fit_poly_pairs=fit_poly_pairs,
                poly_degree=poly_degree_time,
                poly_label=poly_label
            )
            plot_one(
                df_ell=sub,
                x_col="V",
                x_label="$|V|$",
                y_col="peak_ram_mb",
                y_label="Peak RAM (MB)",
                title="",
                error=error,
                secondary_map=v2n,
                secondary_label="$N$",
                output_path=f"{base}_memory.pdf",
                ell=ell,
                annotate_pairs=None,
                log_scale=log_scale_mem,
                show_legend=False,
                fit_poly_pairs=fit_poly_pairs_mem,
                poly_degree=poly_degree_mem,
                poly_label=poly_label
            )
            save_legend(h, l, f"{base}_legend.pdf")

def run_rq4(df, stem, error, annotate_pairs, log_scale,
            log_scale_mem=False):
    df["edge_density_pct"] = (df["edge_density"] * 100).round().astype(int)
    for ell in sorted(df["ell"].unique()):
        for group_name, ds in d_groups(ell):
            sub  = df[(df["ell"] == ell) & (df["d"].isin(ds))]
            if sub.empty:
                continue
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
                log_scale=log_scale_mem,
                show_legend=False
            )
            save_legend(h, l, f"{base}_legend.pdf")

def run_rq5(df, stem, error, annotate_pairs, log_scale, ell=8,
            log_scale_mem=False, speedup=False):
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
        show_legend=False,
        speedup=speedup
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
        log_scale=log_scale_mem,
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
    parser.add_argument("--log-scale-time", action="store_true",
                        help="Use log scale on the y-axis of the time plots")
    parser.add_argument("--log-scale-mem", action="store_true",
                        help="Use log scale on the y-axis of the memory plots")
    parser.add_argument("--fit-poly-time", nargs="+", type=ell_d_pair, default=None,
                        metavar="ell:d",
                        help="Overlay a polynomial fit (a*x^deg + ... + c, degree set by "
                             "--poly-degree) on the time plots for the given (ell,d) pairs. "
                             "Fitted with 1/y weighting so it is not dominated by large-y "
                             "points. Coefficients are printed to stdout. (e.g. --fit-poly-time 8:3)")
    parser.add_argument("--fit-poly-mem", nargs="+", type=ell_d_pair, default=None,
                        metavar="ell:d",
                        help="As --fit-poly-time, but for the memory plots.")
    parser.add_argument("--poly-degree", type=int, default=2,
                        help="Default polynomial fit degree used for both axes when "
                             "--poly-degree-time / --poly-degree-mem are not given "
                             "(default: 2, i.e. a*x^2+b*x+c).")
    parser.add_argument("--poly-degree-time", type=int, default=None,
                        help="Polynomial fit degree for the time plots. Overrides "
                             "--poly-degree for time only (e.g. 2 for quadratic).")
    parser.add_argument("--poly-degree-mem", type=int, default=None,
                        help="Polynomial fit degree for the memory plots. Overrides "
                             "--poly-degree for memory only (e.g. 1 for linear).")
    parser.add_argument("--no-poly-label", action="store_true",
                        help="Draw the polynomial fit curve but omit the equation label "
                             "in the bottom-right corner of the plot.")
    parser.add_argument("--speedup", action="store_true",
                        help="Label each point on the RQ5 time plot with its speedup "
                             "relative to the fewest-threads run of the same series")
    parser.add_argument("--rq5-ell", type=int, default=5,
                        help="Motif length used for the RQ5 thread-scaling run (default: 5)")
    parser.add_argument("--mem-per-thread", action="store_true",
                        help="Plot peak RAM divided by the number of threads spawned "
                             "instead of the total peak RAM of the process")
    args = parser.parse_args()

    sns.set_theme(style="white", font_scale=1.0)

    df   = load(args.input, args.mem_per_thread)
    stem = str(Path(args.output).with_suffix(""))  # strip extension if given

    annotate_pairs     = set(args.max_rank) if args.max_rank is not None else None
    fit_poly_pairs     = set(args.fit_poly_time) if args.fit_poly_time is not None else None
    fit_poly_pairs_mem = set(args.fit_poly_mem) if args.fit_poly_mem is not None else None
    poly_label         = not args.no_poly_label
    # Per-axis degree, each falling back to the shared --poly-degree default.
    poly_degree_time   = args.poly_degree_time if args.poly_degree_time is not None else args.poly_degree
    poly_degree_mem    = args.poly_degree_mem  if args.poly_degree_mem  is not None else args.poly_degree

    if args.rq == "rq3":
        run_rq3(df, stem, args.error, annotate_pairs, args.log_scale_time,
                args.log_scale_mem,
                fit_poly_pairs, fit_poly_pairs_mem,
                poly_degree_time, poly_degree_mem,
                poly_label)
    elif args.rq == "rq4":
        run_rq4(df, stem, args.error, annotate_pairs, args.log_scale_time,
                args.log_scale_mem)
    elif args.rq == "rq5":
        run_rq5(df, stem, args.error, annotate_pairs, args.log_scale_time, args.rq5_ell,
                args.log_scale_mem, args.speedup)

if __name__ == "__main__":
    main()
