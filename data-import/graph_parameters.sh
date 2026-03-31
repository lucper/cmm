#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <edges.csv> <sequences.csv>"
    exit 1
fi

EDGE_FILE=$1
SEQ_FILE=$2

if [ ! -f "$EDGE_FILE" ] || [ ! -f "$SEQ_FILE" ]; then
    echo "Error: One or both files not found."
    exit 1
fi

awk -F';' '
    # --- PASS 1: Process Sequences File ---
    FILENAME == ARGV[2] {
        if (FNR == 1) next # Skip header
        
        total_chars += length($2)
        
        exists[$1] = 1
        next
    }

    # --- PASS 2: Process Edges File ---
    FILENAME == ARGV[1] {
        if (FNR == 1) next # Skip header
        
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
            deg = count[v] + 0 # +0 handles nodes with no edges
            sum_deg += deg

            if (deg < min) min = deg
            if (deg > max) max = deg
        }

        ed = edges / ((vertex_count * (vertex_count - 1)) / 2)

        if (vertex_count == 0) {
            print "No valid data found."
        } else {
            printf "==========================================\n"
            printf "             Graph Statistics             \n"
            printf "==========================================\n"
            printf "Total Vertices (V):      %d\n", vertex_count
            printf "Total Edges (E):         %d\n", edges
            printf "Total Sequence (N):      %d\n", total_chars
            printf "------------------------------------------\n"
            printf "Min Degree:              %d\n", min
            printf "Max Degree:              %d\n", max
            printf "Avg Degree:              %.2f\n", sum_deg / vertex_count
            printf "Edge Density:            %.2f\n", ed
            printf "==========================================\n"
        }
    }
' "$EDGE_FILE" "$SEQ_FILE"
