#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

FILE_PATHS = {
    3: "/storage/group/izg5139/default/xaris/not_zipf/not_zipf_3mers_final.txt",
    4: "/storage/group/izg5139/default/xaris/not_zipf/not_zipf_4mers_final.txt",
    5: "/storage/group/izg5139/default/xaris/not_zipf/not_zipf_5mers_final.txt",
    6: "/storage/group/izg5139/default/xaris/not_zipf/not_zipf_6mers_final.txt",
    7: "/storage/group/izg5139/default/xaris/not_zipf/not_zipf_7mers_final.txt",
    8: "/storage/group/izg5139/default/xaris/not_zipf/not_zipf_8mers_final.txt",
}

OUTPUT_PLOT = "not_zipf_box_plot.png"

def extract_metrics(file_path):
    r2_vals = []
    with open(file_path) as f:
        for line in f:
            parts = line.strip().split()
            for token in parts:
                if token.startswith('R2='):
                    try:
                        r2_vals.append(float(token.split('=', 1)[1]))
                    except ValueError:
                        pass
    return r2_vals

def custom_boxplot(ax, data, positions):
    for i, vals in enumerate(data):
        vals = np.sort(vals)
        q1 = np.percentile(vals, 25)
        q3 = np.percentile(vals, 75)
        iqr = q3 - q1

        # Calculate standard whisker positions
        standard_lower_whisker = q1 - 1.5 * iqr
        
        # Lower whisker: clip extreme negatives
        lower_whisker = max(standard_lower_whisker, -10.0)
        
        # Upper whisker: extend to maximum value
        upper_whisker = vals.max()

        ax.plot([positions[i]] * 2, [q1, q3], color='orange', linewidth=10)

        median = np.percentile(vals, 50)
        cap_width = 0.05
        ax.plot([positions[i] - cap_width, positions[i] + cap_width], [median] * 2, color='black', linewidth=2)

        ax.plot([positions[i], positions[i]], [lower_whisker, q1], color='black')
        ax.plot([positions[i], positions[i]], [q3, upper_whisker], color='black')

        cap_width = 0.1
        ax.plot([positions[i] - cap_width, positions[i] + cap_width], [lower_whisker] * 2, color='black')
        ax.plot([positions[i] - cap_width, positions[i] + cap_width], [upper_whisker] * 2, color='black')

def main():
    metrics_by_k = {}
    for k, file_path in FILE_PATHS.items():
        if not os.path.exists(file_path):
            print(f"Warning: file for k={k} not found at {file_path}, skipping.")
            continue
        r2_vals = extract_metrics(file_path)
        r2_clean = [v for v in r2_vals if not np.isnan(v)]
        if len(r2_clean) < 2:
            print(f"Not enough data for k={k} (R2 count {len(r2_clean)}), skipping.")
            continue
        metrics_by_k[k] = r2_clean

    if not metrics_by_k:
        print("No metrics to process. Exiting.")
        return

    ks_plot = sorted(metrics_by_k.keys())
    data = [metrics_by_k[k] for k in ks_plot]

    fig, ax = plt.subplots(figsize=(10, 6))
    custom_boxplot(ax, data, ks_plot)
    ax.axhline(0, linestyle='--', color='gray')
    ax.set_xticks(ks_plot)
    ax.set_xlabel('k-mer length', fontsize=16)
    ax.set_ylabel(r'$R^2$', fontsize=16)
    ax.tick_params(labelsize=16)

    # Set y-limits to show the full range
    ax.set_ylim(bottom=-10, top=1.05)

    yticks = [-10, -5, -2.5, -1, 0, 1.0]
    ax.set_yticks(yticks)

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"Custom box plot saved to {OUTPUT_PLOT}")

if __name__ == '__main__':
    main()