#!/usr/bin/env python3
"""
Computes the average R2 (excluding negative values) for k = 3, 4, 5, 6 across four taxonomic groups:
  - Archaea
  - Bacteria
  - Viral
  - Eukaryote
Generates two separate heatmap figures (one per distribution) showing taxonomies on the y-axis and k-mer lengths on the x-axis,
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

distributions = {
    'Truncated': 'Results_truncated.xlsx',
    'Zipf-Mandelbrot': 'Results_Zipf_Mandelbrot.xlsx'
}
k_values = [3, 4, 5, 6]
taxonomies = {
    'Archaea': ['archaea'],
    'Bacteria': ['bacteria'],
    'Viral': ['viral'],
    'Eukaryote': ['protozoa', 'vertebrate', 'vertebrate_other', 'fungi', 'plant', 'invertebrate', 'vertebrate_mammalian']
}

def filter_taxonomy(df, tax_list):
    return df if tax_list is None else df[df['Taxonomy'].str.lower().isin(tax_list)]

tax_keys = list(taxonomies.keys())

avg_r2 = {dist: np.zeros((len(tax_keys), len(k_values))) for dist in distributions}
for i, tax in enumerate(tax_keys):
    tax_list = taxonomies[tax]
    for dist, path in distributions.items():
        for j, k in enumerate(k_values):
            df = pd.read_excel(path, sheet_name=f'k{k}')
            df = filter_taxonomy(df, tax_list)
            vals = df['R2'][df['R2'] >= 0].dropna().astype(float)
            avg_r2[dist][i, j] = vals.mean() if not vals.empty else np.nan

all_vals = np.concatenate([avg_r2[d].flatten() for d in distributions])
vmin, vmax = np.nanmin(all_vals), np.nanmax(all_vals)

for dist, matrix in avg_r2.items():
    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(matrix, aspect='auto', vmin=vmin, vmax=vmax, cmap='viridis')
    ax.set_xticks(np.arange(len(k_values)))
    ax.set_xticklabels(k_values, fontsize=14)
    ax.set_yticks(np.arange(len(tax_keys)))
    ax.set_yticklabels(tax_keys, fontsize=14)
    ax.set_xlabel('k-mer length', fontsize=14)

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            if not np.isnan(val):
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', color='white', fontsize=12)

    cbar = fig.colorbar(im, ax=ax, orientation='vertical', pad=0.04)
    cbar.ax.tick_params(labelsize=14)

    plt.tight_layout()
    outname = f'R2_heatmap_{dist.replace(" ", "_")}.png'
    plt.savefig(outname)
    print(f'Heatmap saved to {outname}')
