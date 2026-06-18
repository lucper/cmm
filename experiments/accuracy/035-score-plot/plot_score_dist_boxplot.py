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
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns
from pathlib import Path


def setup_publication_style(use_tex=True, base_fontsize=11, tex_preamble=None):
    """Configure matplotlib for LaTeX-quality output that matches a Times
    (newtx-style) document.

    use_tex=True routes ALL text through the system LaTeX install so body text
    AND math (axis numbers, the percent sign) render in Times exactly as in the
    document. use_tex=False falls back to matplotlib's own mathtext with a
    Times-like serif family -- no TeX needed, close but not identical.

    Font package (tex_preamble): the default is 'mathptmx', which uses the
    standard PostScript Times shipped with every TeX distribution. matplotlib's
    usetex path reads the Type 1 fonts the TeX engine references and embeds
    them; mathptmx's fonts (the classic Adobe 35) always resolve, whereas
    newtx references fonts (ntx-Regular-*) that frequently fail to resolve
    through pdftex.map unless the TeX Gyre Termes fonts and their map entries
    are fully wired up -- the cause of the 'could not be found for TeX font'
    error. mathptmx Times is visually indistinguishable from newtx Times for a
    plot's labels, so the document can still use newtx. To force true newtx
    (after fixing the font map), pass
    tex_preamble=r'\\usepackage{newtxtext}\\usepackage{newtxmath}'.

    Fonts are embedded as TrueType (Type 42), never Type 3, which some
    publishers reject and which renders poorly.

    base_fontsize is the size (pt) text should APPEAR at in the document; size
    the figure to its final width (no scaling in LaTeX) so this is faithful.
    """
    rc = {
        "font.family": "serif",
        "font.size": base_fontsize,
        "axes.labelsize": base_fontsize + 1,
        "axes.titlesize": base_fontsize + 1,
        "xtick.labelsize": base_fontsize,
        "ytick.labelsize": base_fontsize,
        "legend.fontsize": base_fontsize - 1,
        # Embed real outline fonts, not Type 3.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        # Thin, consistent rules for print.
        "axes.linewidth": 0.8,
    }
    if use_tex:
        # mathptmx = standard PostScript Times (text + math); resolves on any
        # TeX install. Override via tex_preamble for newtx etc.
        preamble = tex_preamble if tex_preamble is not None \
            else r"\usepackage{mathptmx}"
        rc.update({
            "text.usetex": True,
            "text.latex.preamble": preamble,
        })
    else:
        rc.update({
            "text.usetex": False,
            "mathtext.fontset": "stix",   # Times-like math without TeX
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman",
                           "DejaVu Serif"],
        })
    matplotlib.rcParams.update(rc)

def tex_escape_percent(s, use_tex):
    r"""Under usetex, a literal '%' starts a LaTeX comment and would swallow the
    rest of the line. Escape it to '\%'. No-op when not using TeX."""
    if use_tex and "%" in s and r"\%" not in s:
        return s.replace("%", r"\%")
    return s

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
             reference_line=None, reference_label=None, use_tex=True):
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
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(bottom=min(ymin, 0),
                    top=max(ymax, reference_line * 1.05))
        ax.yaxis.set_major_formatter(
            PercentFormatter(xmax=100, decimals=0, symbol=""))

    ax.set_xlabel("")
    ax.set_ylabel(tex_escape_percent(score_label, use_tex))
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)

    if reference_line is not None and reference_label is not None:
        handles = [Patch(facecolor=palette[i], edgecolor="black",
                         label=tex_escape_percent(m, use_tex))
                   for i, m in enumerate(methods)]
        handles.append(Line2D([0], [0], color="red", linestyle="--",
                              linewidth=2.0,
                              label=tex_escape_percent(reference_label, use_tex)))
        ax.legend(handles=handles, ncol=len(handles),
                  loc="lower center", bbox_to_anchor=(0.5, 1.02),
                  frameon=False)
    elif n_methods > 1:
        ax.legend(ncol=n_methods,
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
    parser.add_argument("--no-tex", action="store_true",
                        help="Disable usetex; use matplotlib mathtext with a "
                             "Times-like font instead. No TeX install needed, "
                             "close but not an exact newtx match.")
    parser.add_argument("--font-size", type=float, default=11,
                        help="Base font size (pt) as it should APPEAR in the "
                             "document. Size the figure to its final width so "
                             "this stays faithful (no scaling in LaTeX).")
    parser.add_argument("--tex-preamble", default=None,
                        help="Override the usetex font package. Default is "
                             "'\\usepackage{mathptmx}' (standard PostScript "
                             "Times, resolves on any TeX install). For true "
                             "newtx pass "
                             "'\\usepackage{newtxtext}\\usepackage{newtxmath}' "
                             "(requires a correctly wired TeX Gyre Termes font "
                             "map).")
    args = parser.parse_args()

    use_tex = not args.no_tex

    sns.set_style("white")
    setup_publication_style(use_tex=use_tex, base_fontsize=args.font_size,
                            tex_preamble=args.tex_preamble)

    specs = [parse_input_spec(s) for s in args.input]
    df, datasets, methods = load_long_dataframe(specs)
    if df.empty:
        raise SystemExit("Error: no data read from any input file")

    reference_line = None if args.no_reference_line else 100.0

    stem = str(Path(args.output).with_suffix(""))
    try:
        plot_box(df, datasets, methods,
                 output_path=f"{stem}.pdf",
                 score_label=args.score_label,
                 reference_line=reference_line,
                 reference_label=args.reference_label,
                 use_tex=use_tex)
    except RuntimeError as e:
        if use_tex:
            raise SystemExit(
                "Error: LaTeX rendering failed. A working TeX install with the "
                "newtx package is required for --usetex (the default). Install "
                "it (TeX Live: 'newtx' + dependencies), or re-run with --no-tex "
                f"for a TeX-free Times-like fallback.\n\nUnderlying error:\n{e}")
        raise

if __name__ == "__main__":
    main()
