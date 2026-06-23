#!/usr/bin/env gnuplot
# Quality-profile curves, with a separate legend-only mode for shared legends.
#
# Usage:
#   Panel (no key):    gnuplot -c curve.gp plot   <data.tsv> <out> <exact-name> <width_cm> [height_cm]
#   Legend only:       gnuplot -c curve.gp legend <methods.csv> <out> <exact-name> <width_cm> [height_cm]
#
#   data.tsv : cols  rank <tab> proper_method <tab> normscore
#   out      : output basename (may include dir, e.g. img/curve_ecoli)
#   exact    : the exact method's PROPER name, exactly as in the data
#   width/height in cm. Produces <out>.tex + <out>.pdf.
#
# In 'plot' mode the per-panel key is OFF (legend is shared, drawn separately).
# In 'legend' mode only the key renders: no axes, border, data, or labels.

mode  = ARG1
data  = ARG2
out   = ARG3
exact = ARG4
W     = ARG5 + 0
H     = (ARGC >= 6) ? (ARG6 + 0) : (W * 0.72)

eval sprintf("set terminal cairolatex pdf input color size %gcm,%gcm font ',9'", W, H)
set output out . '.tex'
set datafile separator "\t"

# ---- shared styling (identical in both modes so legend matches panels) ----
# Colours keyed by METHOD NAME (not position) so they match the histogram:
#   M-SLIDER -> blue, SEQ-SLIDER -> red. Exact baseline -> green (was red).
# Unknown names fall back to a neutral grey.
H1COL   = "#1f6fb3"   # blue  -> M-SLIDER   (matches histogram)
H2COL   = "#c14a3d"   # red   -> SEQ-SLIDER (matches histogram)
EXCOL   = "#2c8a3d"   # green -> exact baseline
FALLBK  = "#7f7f7f"   # grey  -> any other method name

# Line widths: thicker so the dashed heuristics are visible without zooming.
LW_HEUR  = 3.5
LW_EXACT = 5.0

# Map a method's display name to its colour.
methcol(name) = (name eq "M-SLIDER")   ? H1COL : \
                (name eq "SEQ-SLIDER") ? H2COL : FALLBK

if (mode eq "legend") {
    # --- legend-only: hide everything but the key ---
    # method names come from methods.csv (col 2, quoted), excluding the exact name
    methods = system("awk -F, -v ex=\"".exact."\" '{n=$2; gsub(/\"/,\"\",n)} n!=ex && !(n in s){s[n]; printf \"%s\\n\", n}' ".data)
    NM = words(methods)
    unset border
    unset xtics
    unset ytics
    unset xlabel
    unset ylabel
    set yrange [0:1]
    set xrange [0:1]
    # key centered, horizontal row; this is the only thing that renders
    set key center center horizontal maxrows 1 samplen 3 width 6 spacing 1.5
    # plot the series OUT OF RANGE (x=2, beyond [0:1]) so only key samples show
    plot \
      for [i=1:NM] '+' using (2):(2) with lines lw LW_HEUR dashtype 2 \
         lc rgb methcol(word(methods,i)) title word(methods,i), \
      '+' using (2):(2) with lines lw LW_EXACT lc rgb EXCOL title exact
} else {
    # --- panel: real plot, key OFF (shared legend lives elsewhere) ---
    # method names come from the TSV (col 2, tab), excluding the exact name
    methods = system("awk -F'\t' -v ex=\"".exact."\" '\$2!=ex && !(\$2 in s){s[\$2]; printf \"%s\\n\", \$2}' ".data)
    NM = words(methods)
    set xlabel '\normalsize rank'
    set ylabel '\normalsize norm. $f_{\chi^2}$'
    stats data using 1 nooutput
    set yrange [0:1.03]
    set xrange [1:STATS_max]
    set border 3
    set xtics nomirror
    set xtics (1, STATS_max)
    set ytics nomirror 0.5
    unset key
    plot \
      for [i=1:NM] data using 1:(strcol(2) eq word(methods,i) ? $3 : 1/0) \
         with lines lw LW_HEUR dashtype 2 lc rgb methcol(word(methods,i)) notitle, \
      data using 1:(strcol(2) eq exact ? $3 : 1/0) \
         with lines lw LW_EXACT lc rgb EXCOL notitle
}
