#!/usr/bin/env python3
"""
Grouped box plot of per-(dataset, method) normalized score distributions.

X-axis: datasets. Within each dataset, one box per method (color-coded).

Input is one or more 'dataset:method:path' triples, where each path points at a
030-score output file whose column 3 is the per-rank ratio
(slider_score / exact_score * 100).

This version can emit a PGF figure (--format pgf) for \\input into a LaTeX
document so that ALL text in the figure is typeset by the document's own LaTeX
run -- fonts, math, and the percent sign match the surrounding document
exactly, with nothing embedded. --format pdf reproduces the old standalone
behaviour.
"""

import argparse
import matplotlib

# IMPORTANT: choose the backend BEFORE importing pyplot. For PGF output we must
# be on the 'pgf' backend; do it lazily in main() once we know --format, but
# set a safe default here.
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns
from pathlib import Path


# Default font preamble. mathptmx = standard PostScript Times (text + math),
# resolves on any TeX install. To match a newtx document exactly, pass
# --tex-preamble '\usepackage{newtxtext}\usepackage{newtxmath}'.
DEFAULT_PREAMBLE = r"\usepackage{mathptmx}"


def setup_pgf_style(base_fontsize=11, tex_preamble=None, tex_engine="pdflatex"):
    r"""Configure matplotlib's PGF backend so the saved .pgf inherits the
    document's fonts.

    The key differences from a usetex-into-PDF setup:

    * backend = 'pgf' and we save to a .pgf file. matplotlib writes LaTeX
      drawing commands; text is left as \pgftext nodes that YOUR document
      typesets when you \input the file. Nothing is rasterised and no font is
      embedded into a standalone PDF.

    * pgf.rcfonts = False tells the backend NOT to inject its own font setup,
      so font.family / font.serif here are ignored in favour of whatever the
      preamble (and ultimately your document) defines. This is what makes the
      figure text match the document.

    * pgf.preamble must declare the SAME font packages as your document, since
      matplotlib uses these metrics to lay out (place/rotate/size) the text
      boxes. If the preamble's fonts differ from the document's, spacing can
      drift slightly. Keep them identical.

    * pgf.texsystem selects the engine matplotlib shells out to for metric
      measurement. Use the same engine you compile your document with
      (pdflatex / xelatex / lualatex). newtx works under all three; xelatex /
      lualatex are needed only if your document loads system OpenType fonts via
      fontspec.

    base_fontsize is the size text should APPEAR at in the document. Because a
    \input'd PGF figure is NOT scaled by LaTeX, set the figure to its final
    width and this size is faithful.
    """
    preamble = tex_preamble if tex_preamble is not None else DEFAULT_PREAMBLE

    matplotlib.use("pgf")  # must precede the first real figure
    rc = {
        "pgf.texsystem": tex_engine,
        "pgf.rcfonts": False,      # inherit document fonts; ignore mpl font.*
        "pgf.preamble": preamble,  # same font packages as the document
        "font.family": "serif",
        "font.size": base_fontsize,
        "axes.labelsize": base_fontsize + 1,
        "axes.titlesize": base_fontsize + 1,
        "xtick.labelsize": base_fontsize,
        "ytick.labelsize": base_fontsize,
        "legend.fontsize": base_fontsize - 1,
        "axes.linewidth": 0.8,
    }
    matplotlib.rcParams.update(rc)
    # Re-assert AFTER the bulk update in case anything (seaborn theming, etc.)
    # tried to flip font.family to sans-serif. Under pgf.rcfonts=False this is
    # what maps to \rmfamily; if it becomes 'sans-serif' the backend emits
    # \sffamily and the figure won't match a serif document.
    matplotlib.rcParams["font.family"] = "serif"


def setup_pdf_style(use_tex=True, base_fontsize=11, tex_preamble=None):
    """Standalone-PDF style (the original behaviour). Kept for --format pdf."""
    rc = {
        "font.family": "serif",
        "font.size": base_fontsize,
        "axes.labelsize": base_fontsize + 1,
        "axes.titlesize": base_fontsize + 1,
        "xtick.labelsize": base_fontsize,
        "ytick.labelsize": base_fontsize,
        "legend.fontsize": base_fontsize - 1,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    }
    if use_tex:
        preamble = tex_preamble if tex_preamble is not None else DEFAULT_PREAMBLE
        rc.update({"text.usetex": True, "text.latex.preamble": preamble})
    else:
        rc.update({
            "text.usetex": False,
            "mathtext.fontset": "stix",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman",
                           "DejaVu Serif"],
        })
    matplotlib.rcParams.update(rc)


def tex_escape_percent(s, render_via_tex):
    r"""Under TeX (usetex or pgf), a literal '%' starts a comment. Escape it."""
    if render_via_tex and "%" in s and r"\%" not in s:
        return s.replace("%", r"\%")
    return s


def load_tsv_dataframe(path):
    """Read the 4-column TSV emitted by glob.sh
    (order_idx, dataset, method, score). Dataset x-axis order follows column 1
    (the AWK's xi[] index = datasets.csv order); method/legend order is
    first-seen. Non-numeric rows (e.g. a header) are dropped. Returns
    (df, datasets_order, methods_order)."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["idx", "dataset", "method", "score"])
    df["idx"] = pd.to_numeric(df["idx"], errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["idx", "score"])
    if df.empty:
        return df, [], []

    datasets_order = (df.sort_values("idx")["dataset"]
                        .drop_duplicates().tolist())   # by column-1 index
    methods_order = df["method"].drop_duplicates().tolist()  # first-seen

    df["dataset"] = pd.Categorical(df["dataset"], datasets_order, ordered=True)
    df["method"] = pd.Categorical(df["method"], methods_order, ordered=True)

    print(f"Loaded {len(df)} rows from {path} "
          f"({len(datasets_order)} datasets x {len(methods_order)} methods)")
    return df, datasets_order, methods_order


def plot_box(df, datasets, methods, output_path, score_label,
             reference_line=None, reference_label=None, render_via_tex=True,
             figsize=(7, 4.5), tight_bbox=True, show_legend=True):
    n_methods = len(methods)
    palette = sns.color_palette("tab10", n_colors=max(n_methods, 1))

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

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
    ax.set_ylabel(tex_escape_percent(score_label, render_via_tex))
    ax.tick_params(axis="x", rotation=30)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)

    if not show_legend:
        # Panel figures share one legend drawn once in LaTeX; suppress the
        # per-panel legend so half-width panels don't each carry (and clip) a
        # duplicate. seaborn's boxplot adds a legend by default, so remove it.
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()
    elif reference_line is not None and reference_label is not None:
        handles = [Patch(facecolor=palette[i], edgecolor="black",
                         label=tex_escape_percent(m, render_via_tex))
                   for i, m in enumerate(methods)]
        handles.append(Line2D([0], [0], color="red", linestyle="--",
                              linewidth=2.0,
                              label=tex_escape_percent(reference_label,
                                                       render_via_tex)))
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

    # Guard: under the pgf backend, a sans-serif family silently produces
    # \sffamily labels that won't match a serif document. Catch it here rather
    # than discovering it after a LaTeX build.
    if matplotlib.get_backend().lower() == "pgf":
        fam = matplotlib.rcParams["font.family"]
        if any("sans" in f for f in fam):
            print("Warning: font.family is sans-serif under the pgf backend; "
                  "figure text will use \\sffamily and won't match a serif "
                  "document. Something reset it after setup.")

    # bbox_inches="tight" trims to the artists, so the final box is usually a
    # bit WIDER than figsize (the y-label / tick numbers protrude). That breaks
    # an exact \textwidth fit. With tight_bbox=False the box is exactly figsize;
    # constrained_layout still prevents label clipping by shrinking the axes
    # inward, so nothing is cut off.
    if tight_bbox:
        fig.savefig(output_path, bbox_inches="tight")
    else:
        fig.savefig(output_path)
    plt.close(fig)
    print(f"Saved: {output_path}")


def make_legend_strip(methods, output_path, reference_label=None,
                      render_via_tex=True, width_in=4.98, n_colors=None,
                      tight_bbox=True):
    """Render ONLY the legend as its own PGF, to \\input once above a row of
    side-by-side panels. Uses the same tab10 palette indices as plot_box so the
    swatch colors match the panels exactly. Height auto-sizes to one row."""
    n = n_colors if n_colors is not None else len(methods)
    palette = sns.color_palette("tab10", n_colors=max(n, 1))

    fig, ax = plt.subplots(figsize=(width_in, 0.4), constrained_layout=True)
    ax.axis("off")

    handles = [Patch(facecolor=palette[i], edgecolor="black",
                     label=tex_escape_percent(m, render_via_tex))
               for i, m in enumerate(methods)]
    if reference_label is not None:
        handles.append(Line2D([0], [0], color="red", linestyle="--",
                              linewidth=2.0,
                              label=tex_escape_percent(reference_label,
                                                       render_via_tex)))
    ax.legend(handles=handles, ncol=len(handles), loc="center",
              frameon=False, borderaxespad=0)

    if tight_bbox:
        fig.savefig(output_path, bbox_inches="tight")
    else:
        fig.savefig(output_path)
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True,
                        help="Path to the 4-column TSV from glob.sh "
                             "(order_idx, dataset, method, score). Dataset "
                             "order follows column 1; method order is "
                             "first-seen. Use '-' to read from stdin.")
    parser.add_argument("--output", required=True,
                        help="Output stem (e.g. score_box)")
    parser.add_argument("--format", choices=["pgf", "pdf"], default="pgf",
                        help="pgf: emit .pgf (+ accompanying .pdf images) to "
                             "\\input into LaTeX, text typeset by the document. "
                             "pdf: standalone PDF (old behaviour).")
    parser.add_argument("--score-label", default="Score / exact (%)",
                        help="Y-axis label.")
    parser.add_argument("--reference-label", default="Exact algorithm",
                        help="Legend label for the red 100%% reference line.")
    parser.add_argument("--no-reference-line", action="store_true",
                        help="Suppress the red 100%% reference line and percent "
                             "y-axis formatting.")
    parser.add_argument("--no-tex", action="store_true",
                        help="(pdf format only) Disable usetex; use mathtext "
                             "with a Times-like font. Ignored for pgf.")
    parser.add_argument("--tex-engine", default="pdflatex",
                        choices=["pdflatex", "xelatex", "lualatex"],
                        help="(pgf format) Engine matplotlib uses for text "
                             "metrics. Match your document's compiler.")
    parser.add_argument("--font-size", type=float, default=11,
                        help="Base font size (pt) as it should APPEAR in the "
                             "document. Size the figure to final width.")
    parser.add_argument("--width-in", type=float, default=7.0,
                        help="Figure width in inches (set to your \\linewidth).")
    parser.add_argument("--height-in", type=float, default=4.5,
                        help="Figure height in inches.")
    parser.add_argument("--legend-only", action="store_true",
                        help="Render ONLY the shared legend strip (for "
                             "side-by-side panels). Combine with --no-legend "
                             "on the panels themselves.")
    parser.add_argument("--no-legend", action="store_true",
                        help="Suppress the per-panel legend. Use for "
                             "side-by-side panels that share one legend drawn "
                             "once in the LaTeX figure.")
    parser.add_argument("--no-tight-bbox", action="store_true",
                        help="Save at exactly --width-in/--height-in instead "
                             "of trimming to artists. Use for an exact "
                             "\\textwidth fit; constrained_layout still keeps "
                             "labels from clipping.")
    parser.add_argument("--tex-preamble", default=None,
                        help="Override the font package. Default "
                             "'\\usepackage{mathptmx}'. For newtx pass "
                             "'\\usepackage{newtxtext}\\usepackage{newtxmath}'.")
    args = parser.parse_args()

    # seaborn's set_style/set_theme resets font.family to 'sans-serif'. Apply
    # it FIRST so our serif setting below wins -- otherwise the PGF backend
    # emits every label with \sffamily and the figure never matches a serif
    # document, regardless of --tex-preamble.
    sns.set_style("white")

    if args.format == "pgf":
        setup_pgf_style(base_fontsize=args.font_size,
                        tex_preamble=args.tex_preamble,
                        tex_engine=args.tex_engine)
        render_via_tex = True
        ext = "pgf"
    else:
        use_tex = not args.no_tex
        setup_pdf_style(use_tex=use_tex, base_fontsize=args.font_size,
                        tex_preamble=args.tex_preamble)
        render_via_tex = use_tex
        ext = "pdf"

    if args.input == "-":
        import sys
        import tempfile
        data = sys.stdin.read()
        with tempfile.NamedTemporaryFile("w", suffix=".tsv",
                                         delete=False) as tmp:
            tmp.write(data)
            tsv_path = tmp.name
    else:
        tsv_path = args.input
        if not Path(tsv_path).is_file():
            raise SystemExit(
                f"Error: input TSV not found: {tsv_path}\n"
                "--input now takes the 4-column TSV from glob.sh (a path or "
                "'-' for stdin), not 'dataset:method:path' triples.")
    df, datasets, methods = load_tsv_dataframe(tsv_path)
    if df.empty:
        raise SystemExit(f"Error: no data read from {args.input}")

    reference_line = None if args.no_reference_line else 100.0
    stem = str(Path(args.output).with_suffix(""))

    if args.legend_only:
        ref_label = None if args.no_reference_line else args.reference_label
        try:
            make_legend_strip(methods, output_path=f"{stem}.{ext}",
                              reference_label=ref_label,
                              render_via_tex=render_via_tex,
                              width_in=args.width_in,
                              tight_bbox=not args.no_tight_bbox)
        except RuntimeError as e:
            raise SystemExit(
                "Error: LaTeX/PGF rendering failed.\n\nUnderlying error:\n"
                f"{e}")
        return

    try:
        plot_box(df, datasets, methods,
                 output_path=f"{stem}.{ext}",
                 score_label=args.score_label,
                 reference_line=reference_line,
                 reference_label=args.reference_label,
                 render_via_tex=render_via_tex,
                 figsize=(args.width_in, args.height_in),
                 tight_bbox=not args.no_tight_bbox,
                 show_legend=not args.no_legend)
    except RuntimeError as e:
        raise SystemExit(
            "Error: LaTeX/PGF rendering failed. A working TeX install with the "
            "requested font package is required.\n\nUnderlying error:\n"
            f"{e}")


if __name__ == "__main__":
    main()
