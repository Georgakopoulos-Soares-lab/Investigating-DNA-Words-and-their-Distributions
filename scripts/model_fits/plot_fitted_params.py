#!/usr/bin/env python3
"""
Computes mean and std for k = 3-6, excluding rows with R2 < 0,
and generates five 6-panel figures:
  1. All data
  2. Archaea
  3. Bacteria
  4. Viral
  5. Eukaryote
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

files = {
    'Truncated': 'Results_truncated.xlsx',
    'Zipf-Mandelbrot': 'Results_Zipf_Mandelbrot.xlsx'
}
k_values = [3, 4, 5, 6]
distribution_params = {
    'Truncated': ['alpha', 'lambda', 'scale'],
    'Zipf-Mandelbrot': ['alpha', 'beta', 'scale']
}

taxonomies = {
    'All': None,
    'Archaea': ['archaea'],
    'Bacteria': ['bacteria'],
    'Viral': ['viral'],
    'Eukaryote': ['protozoa', 'vertebrate', 'vertebrate(other)', 'fungi', 'plant', 'invertebrate']
}

colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
indices = np.arange(len(k_values))

for tax_name, tax_list in taxonomies.items():
    stats = {
        dist: {param: {'mean': [], 'std': []} for param in params}
        for dist, params in distribution_params.items()
    }
    print(f"\nTaxonomy: {tax_name}")
    for dist, path in files.items():
        for k in k_values:
            df = pd.read_excel(path, sheet_name=f'k{k}')
            if tax_list is not None:
                df = df[df['Taxonomy'].str.lower().isin(tax_list)]
            neg_count = (df['R2'] < 0).sum()
            print(f"  {dist}, k={k}: excluded {neg_count} negative R2 values")
            df = df[df['R2'] >= 0]
            for param in distribution_params[dist]:
                vals = pd.to_numeric(df[param], errors='coerce').dropna()
                stats[dist][param]['mean'].append(vals.mean() if not vals.empty else np.nan)
                stats[dist][param]['std'].append(vals.std() if not vals.empty else np.nan)

    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 12), sharex=True)
    for row_idx, (param_trunc, param_zipf) in enumerate(zip(distribution_params['Truncated'],
                                                               distribution_params['Zipf-Mandelbrot'])):
        # Left: Truncated
        ax = axes[row_idx, 0]
        means = stats['Truncated'][param_trunc]['mean']
        stds = stats['Truncated'][param_trunc]['std']
        c = colors[row_idx]
        ax.plot(indices, means, marker='o', color=c)
        ax.fill_between(indices,
                        np.array(means) - np.array(stds),
                        np.array(means) + np.array(stds),
                        color=c, alpha=0.3)
        ax.set_title(f'Truncated: {param_trunc}', fontsize=14)
        ax.grid(axis='y', linestyle='--', alpha=0.5)
        ax.tick_params(axis='both', labelsize=14)
        if row_idx == 2:
            ax.set_xticks(indices)
            ax.set_xticklabels([str(k) for k in k_values], fontsize=14)

        # Right: Zipf-Mandelbrot
        ax2 = axes[row_idx, 1]
        means2 = stats['Zipf-Mandelbrot'][param_zipf]['mean']
        stds2 = stats['Zipf-Mandelbrot'][param_zipf]['std']
        c2 = colors[row_idx + 3]
        ax2.plot(indices, means2, marker='o', color=c2)
        ax2.fill_between(indices,
                         np.array(means2) - np.array(stds2),
                         np.array(means2) + np.array(stds2),
                         color=c2, alpha=0.3)
        ax2.set_title(f'Zipf-Mandelbrot: {param_zipf}', fontsize=14)
        ax2.grid(axis='y', linestyle='--', alpha=0.5)
        ax2.tick_params(axis='both', labelsize=14)
        if row_idx == 2:
            ax2.set_xticks(indices)
            ax2.set_xticklabels([str(k) for k in k_values], fontsize=14)

    fig.text(0.5, 0.04, 'k-mer length', ha='center', fontsize=16)
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    outname = f'Parameters_by_k_{tax_name}.png'.replace(' ', '_')
    plt.savefig(outname)
    print(f"Figure saved to {outname}\n")