#!/usr/bin/env python3
"""
Generates eight final hexbin plots (one per distribution per k for k=3,4,5,6) of R² vs. genome size (bp, log scale).
The plot shows data density and:
  - A colorbar indicating point density (log scale).
  - Major & minor grid lines.
  - Spearman correlation analysis (ρ) with a p-value in a textbox.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, rankdata
from matplotlib.colors import LogNorm

files = {
    'Truncated': '/storage/group/izg5139/default/xaris/Investigating_Dna_Words/scripts/model_fits/Results_truncated.xlsx',
    'Zipf-Mandelbrot': '/storage/group/izg5139/default/xaris/Investigating_Dna_Words/scripts/model_fits/Results_Zipf_Mandelbrot.xlsx'
}
k_values = [3, 4, 5, 6]
TICK_FONTSIZE = 16

def interpret_rho(rho):
    abs_r = abs(rho)
    if abs_r < 0.1:
        return 'Very weak'
    elif abs_r < 0.3:
        return 'Weak'
    elif abs_r < 0.5:
        return 'Moderate'
    elif abs_r < 0.7:
        return 'Strong'
    else:
        return 'Very strong'

for dist, path in files.items():
    for k in k_values:
        df = pd.read_excel(path, sheet_name=f'k{k}')
        df = df[df['R2'] >= 0]
        genome_sizes = pd.to_numeric(df['Genome Size (bp)'], errors='coerce').dropna()
        r2_values = df.loc[genome_sizes.index, 'R2'].astype(float)
        if genome_sizes.empty:
            print(f'No data for {dist} k={k} after filtering.')
            continue

        rho, pval = spearmanr(genome_sizes, r2_values)
        
        # If p-value is very small, display as "p = 0"
        if pval < 0.0001:
            p_string = 'p = 0'
        else:
            p_string = f'p = {pval:.4f}' # Otherwise, format to 4 decimal places

        interp = interpret_rho(rho)
        direction = 'increase' if rho > 0 else ('decrease' if rho < 0 else 'no change')
        textstr = (
            f'ρ = {rho:.4f}\n'
            f'{p_string}\n'  
            f'{interp} {direction}'
        )

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.minorticks_on()
        ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.7)
        ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)

        hb = ax.hexbin(
            genome_sizes,
            r2_values,
            gridsize=75,
            cmap='viridis',
            xscale='log',
            mincnt=1,
            norm=LogNorm()
        )
        
        cbar = fig.colorbar(hb, ax=ax)
        cbar.ax.tick_params(labelsize=TICK_FONTSIZE)
        
        ax.set_xlabel('Genome size (bp)', fontsize=16)
        ax.set_ylabel('R²', fontsize=16)

        ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)

        props = dict(boxstyle='round,pad=1.2', facecolor='wheat', alpha=0.8)
        ax.text(
            0.95, 0.05, textstr,
            transform=ax.transAxes,
            fontsize=14,
            verticalalignment='bottom',
            horizontalalignment='right',
            bbox=props
        )

        plt.tight_layout()
        outname = f'R2_vs_genome_size_{dist.replace(" ","_")}_k{k}.png'
        plt.savefig(outname, dpi=300, bbox_inches='tight')
        print(f'Saved plot: {outname}')
        plt.close(fig)