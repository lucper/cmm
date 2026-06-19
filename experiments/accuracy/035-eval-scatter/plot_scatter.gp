#!/usr/bin/env gnuplot
# Score-vs-gap scatter, with a separate legend-only mode for shared legends.
#
# Usage:
#   Panel (no key):  gnuplot -c plot_scatter.gp plot   <methods.csv> <out> <exact-name> <width_cm> [height_cm] <f1> <f2>
#   Legend only:     gnuplot -c plot_scatter.gp legend <methods.csv> <out> <exact-name> <width_cm> [height_cm]
#
#   methods.csv : lines of  key,"Display Name"  (e.g.  m_slider,"M-SLIDER")
#   out         : output basename (may include dir, e.g. img/scatter_ecoli)
#   exact       : the exact method's KEY in methods.csv (e.g. cmm)
#   f1, f2      : the two heuristics' b2a TSVs (cmm_eval output)
#                 cols: idx <tab> selfX <tab> selfY <tab> selfScore
#                       <tab> otherX <tab> otherY <tab> otherScore <tab> similarity
#   width/height in cm. Produces <out>.tex + <out>.pdf.
#
# In 'plot' mode the per-panel key is OFF (legend is shared, drawn separately).
# In 'legend' mode only the key renders: no axes, border, data, or labels.
#
# The two heuristic display names are taken from methods.csv, in file order,
# skipping the exact key. The plot assumes exactly two heuristics.

mode  = ARG1
csv   = ARG2
out   = ARG3
exact = ARG4
W     = ARG5 + 0

# Optional height. We detect its presence by argument COUNT, not by probing the
# value (probing breaks because a trailing filename is non-numeric):
#   legend mode:  ... <W> [H]                  -> ARGC 5 or 6
#   plot  mode:   ... <W> [H] <f1> <f2>         -> ARGC 7 (no H) or 8 (with H)
has_h = 0
f1 = ""
f2 = ""
if (mode eq "legend") {
    has_h = (ARGC >= 6) ? 1 : 0
} else {
    has_h = (ARGC >= 8) ? 1 : 0
    if (has_h) { f1 = ARG7; f2 = ARG8 } else { f1 = ARG6; f2 = ARG7 }
}
H = has_h ? (ARG6 + 0) : (W * 0.85)

eval sprintf("set terminal cairolatex pdf input color size %gcm,%gcm font ',9'", W, H)
set output out . '.tex'
set datafile separator "\t"

# ---- shared styling (identical in both modes so legend matches panels) ----
# Two heuristics: distinct hue AND distinct shape so they separate under overplot.
H1COL = "#1f6fb3"   # blue   -> heuristic 1
H2COL = "#c14a3d"   # red    -> heuristic 2
H1PT  = 7           # filled circle
H2PT  = 9           # filled triangle
H1PS  = 0.45
H2PS  = 0.5         # triangles read smaller; nudge up

# Pull the two heuristic display names from the CSV, in file order, skipping the
# exact KEY. Emits one display name per line.
heur_names = system( \
  "awk -F, -v ex=\"".exact."\" '$1!=ex {n=$2; gsub(/\"/,\"\",n); print n}' ".csv)
N1 = word(heur_names, 1)
N2 = word(heur_names, 2)

# Exact display name (for the x-axis label in plot mode).
exact_name = system( \
  "awk -F, -v ex=\"".exact."\" '$1==ex {n=$2; gsub(/\"/,\"\",n); print n}' ".csv)

if (mode eq "legend") {
    # --- legend-only: hide everything but the key ---
    unset border
    unset xtics
    unset ytics
    unset xlabel
    unset ylabel
    set yrange [0:1]
    set xrange [0:1]
    set key center center horizontal maxrows 1 samplen 2 width 4 spacing 1.5 noenhanced
    # plot OUT OF RANGE (x=2, beyond [0:1]) so only the key samples render
    plot \
      '+' using (2):(2) with points pt H1PT ps H1PS lc rgb H1COL title N1, \
      '+' using (2):(2) with points pt H2PT ps H2PS lc rgb H2COL title N2
} else {
    # --- panel: real scatter, key OFF (shared legend lives elsewhere) ---
    set xlabel '\normalsize '.exact_name.' score $f_{\chi^2}$'
    set ylabel '\normalsize dist.\ to nearest pair $1-s$'

    set yrange [0:1.04]
    set ytics nomirror 0.5
    set border 3
    set xtics nomirror
    # Scores run to ~30000; raw labels collide. Compact SI format (16000->"16k")
    # plus an explicit coarse interval so spacing is deterministic across datasets.
    set format x '%.0s%c'
    set xtics 5000
    unset key

    plot \
      f1 using 4:(1-$8) with points pt H1PT ps H1PS lc rgb H1COL notitle, \
      f2 using 4:(1-$8) with points pt H2PT ps H2PS lc rgb H2COL notitle
}
