#!/usr/bin/env python3
"""
Computes the average R2 and standard deviation for k = 3,4,5,6 for each distribution,
Prints the count of values used, how many were negative, and plots a grouped bar chart of average R2 vs k with ±1 std error bars.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

files = {
    'Truncated': 'Results_truncated.xlsx',
    'Zipf-Mandelbrot': 'Results_Zipf_Mandelbrot.xlsx'
}
output_plot = 'Average_R2.png'
s_values = [3, 4, 5, 6]

avg_r2 = {name: [] for name in files}
std_r2 = {name: [] for name in files}
counts = {name: [] for name in files}

for name, path in files.items():
    for k in s_values:
        sheet = f'k{k}'
        try:
            df = pd.read_excel(path, sheet_name=sheet)
        except Exception as e:
            print(f"Error reading {sheet} from {path}: {e}")
            avg_r2[name].append(np.nan)
            std_r2[name].append(np.nan)
            counts[name].append(0)
            continue

        r2_vals = df['R2'].dropna().astype(float)
        total_count = len(r2_vals)
        neg_count = (r2_vals < 0).sum()

        nonneg_r2 = r2_vals[r2_vals >= 0]
        used_count = len(nonneg_r2)
        mean_r2 = nonneg_r2.mean() if used_count > 0 else np.nan
        std = nonneg_r2.std() if used_count > 0 else np.nan

        counts[name].append(used_count)
        avg_r2[name].append(mean_r2)
        std_r2[name].append(std)

        print(
            f"{name}, k={k}: {used_count} values used ({neg_count} negative),"
            f" avg R2 = {mean_r2:.4f}, std = {std:.4f}"
        )

x = np.arange(len(s_values)) 
width = 0.35

fig, ax = plt.subplots()
for i, name in enumerate(files):
    offsets = x + (i - 0.5) * width
    ax.bar(offsets, avg_r2[name], width,
           yerr=std_r2[name], capsize=5,
           label=name)

ax.set_xlabel('k-mer length', fontsize=16)
ax.set_ylabel('Average R²', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels([str(k) for k in s_values], fontsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.grid(axis='y', linestyle='--', alpha=0.7)

ax.set_ylim(0, 1)
ax.set_yticks(np.arange(0.05, 1.0, 0.1))

ax.legend(
    loc='lower center',      
    bbox_to_anchor=(0.5, -0.35), 
    ncol=len(files),            
    frameon=False,              
    fontsize=14                 
)

plt.tight_layout()
plt.savefig(output_plot, bbox_inches='tight')
print(f"Plot saved to {output_plot}")