import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict

def plot_comparison(data_dict, metric_index, ylabel, filename):
    """
    data_dict: {'MethodName': [(ed, prec, rec), ...]}
    metric_index: 1 for Precision, 2 for Recall
    """
    plt.figure(figsize=(10, 6))
    
    # Define some distinct markers/colors for variety
    markers = ['o', 's', '^', 'D', 'v']
    colors = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#ff7f0e']

    for i, (method, points) in enumerate(data_dict.items()):
        points.sort() # sort by ed
        x = [p[0] for p in points]
        y = [p[metric_index] for p in points]
        
        plt.plot(x, y, label=method, marker=markers[i % len(markers)], 
                 color=colors[i % len(colors)], linewidth=2)

    plt.xlabel('Edge density')
    plt.ylabel(ylabel)
    plt.ylim(0, 1.05) # Extra room for clarity
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='lower', frameon=True, edgecolor='black')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

def main():
    method_groups = defaultdict(list)
    
    # parse from stdin: "Method_A 0.05 0.2 0.3"
    for line in sys.stdin:
        parts = line.split()
        if len(parts) == 4:
            name, param, prec, rec = parts
            method_groups[name].append((float(param), float(prec), float(rec)))
        else:
            print(f"Illegal line: {line}. Skipping.", file=sys.stderr)

    if not method_groups:
        raise SystemExit(f"Could not parse any line from {sys.stdin}")

    plot_comparison(method_groups, 1, "Precision", "comparison_precision.png")
    plot_comparison(method_groups, 2, "Recall", "comparison_recall.png")

if __name__ == "__main__":
    main()
