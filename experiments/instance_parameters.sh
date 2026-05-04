#!/usr/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <nodes.csv> <edges.csv>"
    exit 1
fi

NODE_FILE=$1
EDGE_FILE=$2

awk '
    FILENAME == ARGV[1] {
        if (/^>/) {
            current_id = substr($0, 2)
            sub(/[ \t].*/, "", current_id)  # take first token only
            exists[current_id] = 1
        } else if (current_id != "") {
            total_chars += length($0)
        }
        next
    }

    FILENAME == ARGV[2] {
        if (NF < 2) next
        count[$1]++
        count[$2]++
        edges++

        exists[$1] = 1
        exists[$2] = 1
    }

    END {
        min = 1e18
        max = 0
        sum_deg = 0
        vertex_count = 0

        for (v in exists) {
            vertex_count++
            deg = count[v] + 0
            sum_deg += deg

            if (deg < min) min = deg
            if (deg > max) max = deg
        }

        ed = (1.0 * edges) / ((vertex_count * (vertex_count - 1)) / 2)

        if (vertex_count == 0) {
            print "No valid data found."
        } else {
            printf "Number of nodes (V)      %d\n", vertex_count
            printf "    Min degree           %d\n", min
            printf "    Max degree           %d\n", max
            printf "    Avg degree           %.2f\n", sum_deg / vertex_count
            printf "Number of edges (E)      %d\n", edges
            printf "    Edge density         %.3f\n", ed
            printf "Total label length (N)   %d\n", total_chars
            printf "    Avg label length     %.2f\n", total_chars / vertex_count
        }
    }
' "$NODE_FILE" "$EDGE_FILE"
