#!/usr/bin/env gnuplot
# Kernel density estimate of distance-to-nearest-heuristic-pair (1 - similarity)
# over S_B, with a separate legend-only mode for a shared legend atop a grid.
#
# Usage:
#   Panel (no key):  gnuplot -c plot_density.gp plot   <methods.csv> <out> <exact-key> <width_cm> <height_cm> <f1> <f2> [thr]
#   Legend only:     gnuplot -c plot_density.gp legend <methods.csv> <out> <exact-key> <width_cm> [height_cm]
#
#   methods.csv : lines of  key,"Display Name"
#   out         : output basename (may include dir, e.g. img/dens_ecoli)
#   exact       : the exact method's KEY in methods.csv (e.g. cmm) -- used only
#                 to pick the two heuristic display names for the legend.
#   f1, f2      : the two heuristics' b2a TSVs (cmm_eval output); similarity is
#                 column 8, so distance is (1 - $8).
#   thr         : (optional) threshold in (0,1]. Draws a black dashed vertical
#                 line there with the position as a percentage above it. Omit for none.
#   Produces <out>.tex + <out>.pdf.
#
# y-axis is a NORMALIZED density (gnuplot 'smooth kdensity', default bandwidth):
# the area under each curve integrates to 1, so panels with different total
# counts are directly comparable in SHAPE. The two heuristics are NOT comparable
# in magnitude (each integrates to 1 independently). Domain is [0,1]; with the
# default kernel the curve near the boundaries is a smoothed estimate.

mode  = ARG1
csv   = ARG2
out   = ARG3
exact = ARG4
W     = ARG5 + 0

# Argument layout (height REQUIRED in plot mode so the only optional is trailing):
#   legend:  mode csv out exact W [H]                 -> ARGC 5 or 6
#   plot:    mode csv out exact W H f1 f2 [thr]        -> ARGC 8 (no thr) or 9
f1 = "" ; f2 = "" ; thr = 0
if (mode eq "legend") {
    has_h = (ARGC >= 6) ? 1 : 0
    H = has_h ? (ARG6 + 0) : (W * 0.85)
} else {
    H  = ARG6 + 0
    f1 = ARG7
    f2 = ARG8
    if (ARGC >= 9) { thr = ARG9 + 0 }   # 0 => no threshold line
}

eval sprintf("set terminal cairolatex pdf input color size %gcm,%gcm font ',9'", W, H)
set output out . '.tex'
set datafile separator "\t"

# ---- shared styling (identical in both modes so legend matches panels) ----
# Each series: an opaque coloured line for the curve, with a faint translucent
# fill of the same hue beneath it. Cairo blends the translucent fills on overlap.
H1FILL = "#401f6fb3"   # blue, faint fill   -> heuristic 1 area
H2FILL = "#40c14a3d"   # red,  faint fill   -> heuristic 2 area
H1COL  = "#1f6fb3"     # blue, opaque       -> heuristic 1 line
H2COL  = "#c14a3d"     # red,  opaque       -> heuristic 2 line
LW = 2.5               # curve line width

# Heuristic display names from the CSV (file order, skipping the exact key).
heur_names = system( \
  "awk -F, -v ex=\"".exact."\" '$1!=ex {n=$2; gsub(/\"/,\"\",n); print n}' ".csv)
N1 = word(heur_names, 1)
N2 = word(heur_names, 2)

set style fill solid 1.0 noborder

if (mode eq "legend") {
    # --- legend-only: hide everything but the key ---
    unset border ; unset xtics ; unset ytics ; unset xlabel ; unset ylabel
    set yrange [0:1] ; set xrange [0:1]
    set key center center horizontal maxrows 1 samplen 2 width 6 spacing 1.5 noenhanced
    # lines drawn out of range (x=2) so only the key samples show
    plot \
      '+' using (2):(1) with lines lw LW lc rgb H1COL title N1, \
      '+' using (2):(1) with lines lw LW lc rgb H2COL title N2
} else {
    # --- panel: normalized KDE of distance, two overlaid series ---
    set xlabel '\normalsize dist.\ to closest'
    set ylabel '\normalsize density'

    set xrange [0:1]
    set xtics nomirror ("0" 0, "0.5" 0.5, "1" 1)
    set border 3
    set yrange [0:*]          # autoscale; densities are comparable across panels
    set ytics nomirror 1 scale 0.5   # sparse labeled ticks, short marks, no minors
    unset mytics                      # no minor tick marks cluttering the axis
    unset key

    # --- threshold line + percentage label (only if thr in (0,1]) ---
    if (thr > 0) {
        pct = sprintf("\\footnotesize %d\\%%", int(thr*100 + 0.5))
        set arrow from thr, graph 0 to thr, graph 1 nohead dt 2 lw 1.5 lc rgb "black" front
        set label 1 pct at thr, graph 1.06 center front
    }

    # Per-series sample size N, so each KDE can be weighted by 1/N -> area = 1
    # (a normalized probability density, comparable across panels regardless of
    # count). 'smooth kdensity' uses gnuplot's default bandwidth.
    stats f1 using (1-$8) nooutput ; N1c = STATS_records
    stats f2 using (1-$8) nooutput ; N2c = STATS_records

    # Each series: only the faint filled area under the normalized KDE, NO border
    # line (matches the borderless histogram style). Fill alpha comes from the
    # colour's leading byte (#40..) so the two areas blend where they overlap.
    plot \
      f1 using (1-$8):(1.0/N1c) smooth kdensity with filledcurves y1=0 lc rgb H1FILL notitle, \
      f2 using (1-$8):(1.0/N2c) smooth kdensity with filledcurves y1=0 lc rgb H2FILL notitle
}
