#!/usr/bin/env gnuplot

set terminal cairolatex pdf input color size 14cm,9cm font ",11"
set output ARG2 . '.tex'
set datafile separator "\t"

set style data boxplot
set style boxplot nooutliers
set style fill solid 0.4 border -1
set boxwidth 0.25

set ylabel '\normalsize normalised score (\%)'
set yrange [0:105]
set border 2
set ytics nomirror
unset xtics
set xtics nomirror scale 0 rotate by 30 right
#set key outside top center horizontal maxcols 3 samplen 3 width 6
set key outside top center horizontal maxrows 1 samplen 2 width 5 spacing 1

set lmargin 8
set bmargin 6
set tmargin 4
set rmargin 2

stats ARG1 using 1 nooutput
ND = STATS_max
methods = system("awk -F'\t' '!(\$3 in s){s[\$3]; printf \"%s \", \$3}' ".ARG1)
NM = words(methods)

bw = 0.30
sep = 0.36
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
