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
    sns.set_theme(style="ticks", context="paper", font_scale=2.0)

    rows = [
        {'Method': method, 'Edge density': p[0], ylabel: p[metric_index]}
        for method, points in data_dict.items()
        for p in points
    ]
    df = pd.DataFrame(rows)

    plt.figure(figsize=(10, 6))

    plot = sns.lineplot(
        data=df,
        x='Edge density',
        y=ylabel,
        hue='Method',
        style='Method',
        markers=True,
        dashes=False,
        linewidth=2.5,
        markersize=10
    )

    plt.ylim(0, 1.05)
    plt.legend(title=None, frameon=True, edgecolor='black')

    sns.despine()

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
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
