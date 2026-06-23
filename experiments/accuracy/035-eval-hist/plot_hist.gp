#!/usr/bin/env gnuplot
# Histogram of distance-to-nearest-heuristic-pair (1 - similarity) over S_B,
# with a separate legend-only mode for a shared legend atop a grid.
#
# Usage:
#   Panel (no key):  gnuplot -c plot_hist.gp plot   <methods.csv> <out> <exact-key> <width_cm> <height_cm> <f1> <f2> [thr]
#   Legend only:     gnuplot -c plot_hist.gp legend <methods.csv> <out> <exact-key> <width_cm> [height_cm]
#
#   methods.csv : lines of  key,"Display Name"
#   out         : output basename (may include dir, e.g. img/hist_ecoli)
#   exact       : the exact method's KEY in methods.csv (e.g. cmm) -- used only
#                 to pick the two heuristic display names for the legend.
#   f1, f2      : the two heuristics' b2a TSVs (cmm_eval output); similarity is
#                 column 8, so distance is (1 - $8).
#   thr         : (optional) threshold in {0.1,0.2,...,1.0}. If given, draws a
#                 black dashed vertical line there, labels it with the position
#                 as a percentage centered above the line, and prints+annotates
#                 the count of pairs with distance >= thr for each series (in the
#                 series colour, to the right of the line). Omit for no line.
#   Produces <out>.tex + <out>.pdf.
#
# Fixed 10 bins over the [0,1] distance domain:
#   [0,0.1) [0.1,0.2) ... [0.9,1.0]   (last CLOSED at 1.0)
# y-axis is raw COUNT (no normalization); each panel autoscales its own y.
# Two heuristics OVERLAID at the same bin centres with solid translucent fills
# (cairo blends solid alpha; patterns are not used so transparency works).

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
# Solid translucent fills (alpha byte ~ 55%); cairo BLENDS solid alpha, so the
# two overlaid series mix where they coincide. Opaque variants used for text.
H1FILL = "#8c1f6fb3"   # blue, translucent  -> heuristic 1 bars
H2FILL = "#8cc14a3d"   # red,  translucent  -> heuristic 2 bars
H1COL  = "#1f6fb3"     # blue, opaque       -> heuristic 1 count label
H2COL  = "#c14a3d"     # red,  opaque       -> heuristic 2 count label

# Ten fixed bins of width 0.1, both series at the SAME bin centre (overlaid).
BINW = 0.1
BOXW = 0.09           # slightly under bin width for a thin gap between bins
# Quantize distance d in [0,1] to bin centre (0.05,0.15,...,0.95); d==1.0 -> 0.95
# (last bin closed) instead of spilling into a phantom 11th bin.
binc(d) = (d >= 1.0 ? 0.95 : (floor(d / BINW) + 0.5) * BINW)

# Heuristic display names from the CSV (file order, skipping the exact key).
heur_names = system( \
  "awk -F, -v ex=\"".exact."\" '$1!=ex {n=$2; gsub(/\"/,\"\",n); print n}' ".csv)
N1 = word(heur_names, 1)
N2 = word(heur_names, 2)

set style fill solid 1.0 border lt -1

if (mode eq "legend") {
    # --- legend-only: hide everything but the key ---
    unset border ; unset xtics ; unset ytics ; unset xlabel ; unset ylabel
    set yrange [0:1] ; set xrange [0:1]
    set key center center horizontal maxrows 1 samplen 2 width 4 spacing 1.5 noenhanced
    # boxes drawn out of range (x=2) so only the key swatches show
    set boxwidth BOXW
    plot \
      '+' using (2):(1) with boxes fs solid 1.0 lc rgb H1FILL title N1, \
      '+' using (2):(1) with boxes fs solid 1.0 lc rgb H2FILL title N2
} else {
    # --- panel: 10-bin count histogram, two translucent overlaid series ---
    set xlabel '\normalsize dist.\ to closest'
    set ylabel '\normalsize \# motif pairs'

    set xrange [0:1]
    set xtics nomirror ("0" 0, "0.5" 0.5, "1" 1)
    set border 3
    set boxwidth BOXW
    set format y '%.0s%c'    # y counts as 1k, 10k, ... when >= 1000
    unset key

    # Bin each series to a temp table (bin-centre, count): used for (a) nice
    # y-ticks and (b) the tail counts past the threshold.
    tmp1 = out . '.b1.dat'
    tmp2 = out . '.b2.dat'
    set table tmp1
      plot f1 using (binc(1-$8)):(1) smooth freq
    unset table
    set table tmp2
      plot f2 using (binc(1-$8)):(1) smooth freq
    unset table

    # Max bar height across both series -> "nice" y-tick interval (~3 ticks).
    # Temp tables are whitespace-separated; switch separator to read them.
    set datafile separator whitespace
    stats tmp1 using 2 nooutput ; m1 = STATS_max
    stats tmp2 using 2 nooutput ; m2 = STATS_max
    set datafile separator "\t"
    ymax = (m1 > m2) ? m1 : m2
    yraw = ymax / 3.0
    yp10 = 10.0 ** floor(log10(yraw))
    yn   = yraw / yp10
    ynice = (yn <= 1) ? 1 : (yn <= 2) ? 2 : (yn <= 5) ? 5 : 10
    ystep = ynice * yp10
    set yrange [0 : ymax * 1.20]   # headroom so the percentage label clears the top
    set ytics nomirror ystep

    # --- threshold line + percentage label (only if thr in (0,1]) ---
    if (thr > 0) {
        pct = sprintf("\\footnotesize %d\\%%", int(thr*100 + 0.5))
        # vertical dashed black line at thr (no head)
        set arrow from thr, graph 0 to thr, graph 1 nohead dt 2 lw 1.5 lc rgb "black" front
        # percentage centred ABOVE the line, lifted clear of the plot top so it
        # never overlaps the line or the top bars (headroom added to yrange below)
        set label 1 pct at thr, graph 1.10 center front
    }

    # Overlaid translucent bars at the SAME bin centre.
    plot \
      f1 using (binc(1-$8)):(1) smooth freq with boxes \
         fs solid 1.0 lc rgb H1FILL notitle, \
      f2 using (binc(1-$8)):(1) smooth freq with boxes \
         fs solid 1.0 lc rgb H2FILL notitle
}
