#!/bin/bash

awk '
FILENAME ~ /datasets\.csv$/ {
    split($0,a,","); gsub(/"/,"",a[2])
    dname[a[1]]=a[2]
    xi[a[1]] = ++nx
    next
}
FILENAME ~ /methods\.csv$/  {
    split($0,a,","); gsub(/"/,"",a[2])
    mname[a[1]]=a[2]
    next
}
FNR==1 {
    split(FILENAME, p, "/"); split(p[length(p)], f, ".")
    rawds = f[1]; method = mname[f[2]]
}
{ OFS="\t"; if ($3 != "") print xi[rawds], dname[rawds], method, $3 }
' "$@"
