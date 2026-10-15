#!/usr/bin/env gnuplot
# Time or peak memory of one panel, one line per d, against V with N on the top axis (RQ3), the
# edge density (RQ4), or the number of threads (RQ5).
#
# Usage: gnuplot -c plot_efficiency.gp <mode> <rq> <dat> <out> <W> <H>
#   mode   time, memory, or legend
#   rq     rq3, rq4, or rq5
#   dat    a panel file of make_panels.py
#   W, H   size in cm of the plot area, or of the whole legend in legend mode

mode = ARG1 ; rq = ARG2 ; dat = ARG3 ; out = ARG4 ; W = ARG5 + 0 ; H = ARG6 + 0

top  = (rq eq "rq3")
logy = (rq ne "rq4")
xlab = (rq eq "rq3") ? '$|V|$' : (rq eq "rq4") ? 'Edge density (\%)' : 'Number of threads'

set datafile separator "\t"
set datafile columnheaders
do for [i = 1:5] {
    set style line i lw 2 ps 0.6 pt (word("7 5 9 13 11", i) + 0) lc rgb word("#1f6fb3 #e08214 #2c8a3d #c14a3d #7b3294", i)
}

if (mode eq "legend") {
    eval sprintf("set terminal cairolatex pdf input color size %gcm,%gcm font ',9'", W, H)
    set output out . '.tex'
    stats dat using 1 nooutput
    unset border ; unset xtics ; unset ytics
    set xrange [0:1] ; set yrange [0:1]
    set key center center horizontal maxrows 1 samplen 2 width 4 spacing 1.5
    plot for [i = 0:STATS_blocks - 1] dat index i using (2):(2) with linespoints ls i + 1 title columnheader(1)
    exit
}

# Margins in cm, as in all plots of the experiments; the top one holds the N axis on RQ3.
LM = 1.2 ; RM = 0.3 ; BM = 1.0 ; TM = top ? 1.4 : 0.5
whole_bp(cm) = ceil(cm / 2.54 * 72) / 72.0 * 2.54
PW = W ; PH = H ; W = whole_bp(LM + PW + RM) ; H = whole_bp(BM + PH + TM)
eval sprintf("set terminal cairolatex pdf input color size %gcm,%gcm font ',9'", W, H)
set output out . '.tex'
set lmargin at screen LM / W ; set rmargin at screen (LM + PW) / W
set bmargin at screen BM / H ; set tmargin at screen (BM + PH) / H

set border (top ? 7 : 3)
set offsets 0.5, 0.5, 0, 0
set xtics nomirror
set ytics nomirror
if (top) { set x2tics nomirror rotate by 30 left ; set link x2 ; set x2label '\normalsize $N$' }
if (logy) { set logscale y ; set format y '$10^{%T}$' }
set xlabel '\normalsize '.xlab
set ylabel (mode eq "time") ? '\normalsize Time (s)' : '\normalsize Peak RAM (MB)'
unset key
col = (mode eq "time") ? 3 : 4

plot for [i = 0:*] dat index i using 0:col:xtic(1):x2tic(2) with linespoints ls i + 1, \
     dat using 0:col:(mode eq "time" && strcol(5) ne "-" ? sprintf('\\scriptsize\\texttt{%s}', strcol(5)) : "") \
        with labels offset 0, 0.7 tc rgb "#7f7f7f"
