#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

Taxonomies = True   # Set to False for single overall plot, True for per-taxonomy plots + combined
ks        = [3, 4, 5, 6, 7, 8]
base_path = '/storage/group/izg5139/default/xaris'

overall_template = 'avg_{0}mers.txt'
overall_output   = 'gini_kmers.png'

taxonomies = ['archaea', 'bacteria', 'eukaryote', 'viral']
color_map = {
    'archaea': 'cyan',
    'bacteria': 'red',
    'eukaryote': 'green',
    'viral': 'purple'
}

def compute_gini(counts):
    """counts: positive 1D array -> Gini in [0,1], NaN if no data"""
    x = np.sort(counts)
    n = x.size
    if n == 0:
        return np.nan
    return (2.0 * np.sum(np.arange(1, n+1) * x) / (n * x.sum())) - (n + 1) / n

def gather_gini(file_template):
    gini_vals = []
    for k in ks:
        fp = os.path.join(base_path, file_template.format(k))
        if not os.path.isfile(fp):
            print(f"Warning: not found {fp}")
            gini_vals.append(np.nan)
            continue
        df = pd.read_csv(fp, delim_whitespace=True, header=None, names=['avg','std'])
        df['avg'] = pd.to_numeric(df['avg'], errors='coerce')
        df = df.dropna(subset=['avg'])
        gini_vals.append(compute_gini(df['avg'].values))
    return gini_vals

def plot_gini(ks, gini_vals, color, label, output_png=None):
    plt.figure(figsize=(8,5))
    plt.plot(ks, gini_vals, marker='o', linewidth=2, color=color, label=label)
    for xi, yi in zip(ks, gini_vals):
        if not np.isnan(yi):
            plt.text(xi, yi+0.01, f"{yi:.2f}", ha='center', fontsize=10)
    plt.xlabel('k-mer length', fontsize=16)
    plt.ylabel('Gini Coefficient', fontsize=16)
    plt.xticks(ks, fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim(0,1)
    plt.grid(True, linestyle='--', alpha=0.5)
    if label:
        plt.legend()
    plt.tight_layout()
    if output_png:
        plt.savefig(output_png, dpi=300)
        plt.close()
        print(f"Saved plot to {output_png}")
    else:
        plt.show()

if __name__ == '__main__':
    if not Taxonomies:
        # Single overall plot
        gini_vals = gather_gini(overall_template)
        plot_gini(ks, gini_vals, 'black', None, overall_output)
    else:
        # Per-taxonomy plots
        all_ginis = {}
        for tax in taxonomies:
            template = f"{tax}_avg_{{0}}mers.txt"
            outp = f"gini_kmers_{tax}.png"
            vals = gather_gini(template)
            all_ginis[tax] = vals
            plot_gini(ks, vals, color_map[tax], None, outp)
        # Combined plot
        combined_out = 'gini_kmers_all_taxonomies.png'
        plt.figure(figsize=(8,5))
        for tax, vals in all_ginis.items():
            plt.plot(ks, vals, marker='o', linewidth=2, color=color_map[tax], label=tax)
        plt.xlabel('k-mer length', fontsize=16)
        plt.ylabel('Gini Coefficient', fontsize=16)
        plt.xticks(ks, fontsize=13)
        plt.yticks(fontsize=13)
        plt.ylim(0,1)
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.legend(title='Taxonomy', fontsize=13)
        plt.tight_layout()
        plt.savefig(combined_out, dpi=300)
        plt.close()
        print(f"Saved combined plot to {combined_out}")