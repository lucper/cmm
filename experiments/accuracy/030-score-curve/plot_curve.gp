#!/usr/bin/env gnuplot

mode  = ARG1
csv   = ARG2

# Argument layout differs by mode only in whether a data file precedes <out>:
#   plot:    mode csv data out exact W [H]   -> ARGC 7 or 8
#   legend:  mode csv      out exact W [H]   -> ARGC 6 or 7
if (mode eq "legend") {
    data  = ""
    out   = ARG3
    exact = ARG4
    W     = ARG5 + 0
    H     = (ARGC >= 6) ? (ARG6 + 0) : (W * 0.72)
} else {
    data  = ARG3
    out   = ARG4
    exact = ARG5
    W     = ARG6 + 0
    H     = (ARGC >= 7) ? (ARG7 + 0) : (W * 0.72)
}

# Resolve the exact method's DISPLAY name from its key via methods.csv. The
# data TSV's column 2 holds display names, so we match against this resolved
# name (and use it as the exact legend title).
exact_name = system( \
  "awk -F, -v k=\"".exact."\" '$1==k {n=$2; gsub(/\"/,\"\",n); print n; exit}' ".csv)

eval sprintf("set terminal cairolatex pdf input color size %gcm,%gcm font ',9'", W, H)
set output out . '.tex'
set datafile separator "\t"

# ---- shared styling (identical in both modes so legend matches panels) ----
# Colours keyed by METHOD NAME (not position) so they match the histogram:
#   M-SLIDER -> blue, SEQ-SLIDER -> red. Exact baseline -> green.
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
    # heuristic display names come from methods.csv (col 2, quoted), excluding
    # the exact KEY (col 1).
    methods = system("awk -F, -v ex=\"".exact."\" '$1!=ex {n=$2; gsub(/\"/,\"\",n); print n}' ".csv)
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
      '+' using (2):(2) with lines lw LW_EXACT lc rgb EXCOL title exact_name
} else {
    # --- panel: real plot, key OFF (shared legend lives elsewhere) ---
    # method names come from the TSV (col 2, tab), excluding the exact display name
    methods = system("awk -F'\t' -v ex=\"".exact_name."\" '\$2!=ex && !(\$2 in s){s[\$2]; printf \"%s\\n\", \$2}' ".data)
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
      data using 1:(strcol(2) eq exact_name ? $3 : 1/0) \
         with lines lw LW_EXACT lc rgb EXCOL notitle
}
