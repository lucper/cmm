#!/usr/bin/env gnuplot
# Histogram of distance-to-nearest-heuristic-pair (1 - similarity) over S_B,
# with a separate legend-only mode for a shared legend atop a grid.
#
# Usage:
#   Panel (no key):  gnuplot -c plot_hist.gp plot   <methods.csv> <out> <exact-key> <width_cm> [height_cm] <f1> <f2>
#   Legend only:     gnuplot -c plot_hist.gp legend <methods.csv> <out> <exact-key> <width_cm> [height_cm]
#
#   methods.csv : lines of  key,"Display Name"
#   out         : output basename (may include dir, e.g. img/hist_ecoli)
#   exact       : the exact method's KEY in methods.csv (e.g. cmm) -- used only
#                 to pick the two heuristic display names for the legend.
#   f1, f2      : the two heuristics' b2a TSVs (cmm_eval output); similarity is
#                 column 8, so distance is (1 - $8).
#   Produces <out>.tex + <out>.pdf.
#
# Fixed 5 bins over the [0,1] distance domain:
#   [0,0.2) [0.2,0.4) [0.4,0.6) [0.6,0.8) [0.8,1.0]   (last CLOSED at 1.0)
# y-axis is raw COUNT (no normalization); each panel autosciales its own y.
# Two heuristics overlaid: colored, semi-transparent, different fill patterns.

mode  = ARG1
csv   = ARG2
out   = ARG3
exact = ARG4
W     = ARG5 + 0

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
# Opaque colors: cairo renders pattern fills opaque regardless of alpha, so
# transparency does not blend here. Instead we OFFSET the two series within each
# bin so both are always visible, and give them clearly different patterns.
H1COL = "#1f6fb3"   # blue -> heuristic 1
H2COL = "#c14a3d"   # red  -> heuristic 2
H1PAT = 1           # fillstyle pattern: dense diagonal /////
H2PAT = 6           # fillstyle pattern: opposite-direction / cross feel \\\\\

# Five fixed bins of width 0.2. Each series sits in a half-width box, the two
# offset left/right of the bin center so they sit side by side, not stacked.
BINW   = 0.2
BOXW   = 0.085
OFF    = 0.045      # half-offset: series1 at center-OFF, series2 at center+OFF
# Quantize distance d in [0,1] to bin center (0.1,0.3,0.5,0.7,0.9); d==1.0 -> 0.9
# (last bin closed). The +/-OFF shift separates the two overlaid series.
binc(d)  = (d >= 1.0 ? 0.9 : (floor(d / BINW) + 0.5) * BINW)
binL(d)  = binc(d) - OFF
binR(d)  = binc(d) + OFF

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
      '+' using (2):(1) with boxes fs pattern H1PAT lc rgb H1COL title N1, \
      '+' using (2):(1) with boxes fs pattern H2PAT lc rgb H2COL title N2
} else {
    # --- panel: 5-bin count histogram, two opaque hatched series, side by side ---
    set xlabel '\normalsize dist.\ to nearest'
    set ylabel '\normalsize \# motif pairs'

    set xrange [0:1]
    set xtics nomirror ("0" 0, "0.5" 0.5, "1" 1)
    set border 3
    set boxwidth BOXW
    unset key

    # Bin each series to a temp table (bin-center, count) so we can (a) read the
    # true max bar height for nice y-ticks and (b) plot from stable data.
    tmp1 = out . '.b1.dat'
    tmp2 = out . '.b2.dat'
    set table tmp1
      plot f1 using (binc(1-$8)):(1) smooth freq
    unset table
    set table tmp2
      plot f2 using (binc(1-$8)):(1) smooth freq
    unset table

    # Max count across both series -> choose a "nice" y-tick interval (~3 ticks).
    # The temp tables are space-separated (gnuplot 'set table' format), so switch
    # the datafile separator away from tab before statting them.
    set datafile separator whitespace
    stats tmp1 using 2 nooutput ; m1 = STATS_max
    stats tmp2 using 2 nooutput ; m2 = STATS_max
    set datafile separator "\t"   # restore for the data files below
    ymax = (m1 > m2) ? m1 : m2
    yraw = ymax / 3.0
    yp10 = 10.0 ** floor(log10(yraw))
    yn   = yraw / yp10
    ynice = (yn <= 1) ? 1 : (yn <= 2) ? 2 : (yn <= 5) ? 5 : 10
    ystep = ynice * yp10
    set yrange [0 : ymax * 1.08]
    set ytics nomirror ystep

    # Re-bin with the L/R offset for the actual side-by-side draw.
    plot \
      f1 using (binL(1-$8)):(1) smooth freq with boxes \
         fs pattern H1PAT lc rgb H1COL notitle, \
      f2 using (binR(1-$8)):(1) smooth freq with boxes \
         fs pattern H2PAT lc rgb H2COL notitle
}
