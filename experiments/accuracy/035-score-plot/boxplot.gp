#!/usr/bin/env gnuplot

set terminal pdfcairo size 9,6 enhanced font "Helvetica,12"
set output 'boxplot.pdf'
set datafile separator "\t"

set style data boxplot
#set style boxplot outliers pointtype 7
set style boxplot range 1e6 nooutliers
set style fill solid 0.4 border -1
set boxwidth 0.25

set ylabel "normalised score (%)"
set yrange [0:105]
set border 2
set ytics nomirror
unset xtics
set xtics nomirror scale 0
set xtics rotate by 30 right
set key outside top center horizontal maxcols 3 samplen 3 width 6

stats ARG1 using 1 nooutput
ND = STATS_max

methods = system("awk -F'\t' '!(\$3 in s){s[\$3]; printf \"%s \", \$3}' ".ARG1)
NM = words(methods)              # expected: 2

bw = 0.30
sep = 0.36                       # centre-to-centre gap between the NM boxes
offset(j) = (j - (NM+1)/2.0) * sep

array COL[2]
COL[1] = "#4e79a7"
COL[2] = "#e15759"

set for [i=1:ND] xtics ( \
  system(sprintf("awk -F'\t' '$1==%d{print $2; exit}' %s", i, ARG1)) i )

plot \
  for [j=1:NM] for [i=1:ND] ARG1 \
     using (i+offset(j)):($1==i && stringcolumn(3) eq word(methods,j) ? $4 : 1/0):(bw) \
     with boxplot lc rgb COL[j] notitle, \
  for [j=1:NM] NaN with boxes lc rgb COL[j] title word(methods,j), \
  100 with lines lc rgb "red" dashtype 3 lw 2 title "Exact"
