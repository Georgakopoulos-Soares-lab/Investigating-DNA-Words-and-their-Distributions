#!/usr/bin/env python3
"""
Generates eight hexbin plots (one per distribution per k for k=3,4,5,6) of R² vs. GC Content (%).
Each plot shows data density and:
  - A colorbar indicating point density (log scale).
  - Major & minor grid lines.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

files = {
    'Truncated': '/storage/group/izg5139/default/xaris/Investigating_Dna_Words/scripts/model_fits/Results_truncated.xlsx',
    'Zipf-Mandelbrot': '/storage/group/izg5139/default/xaris/Investigating_Dna_Words/scripts/model_fits/Results_Zipf_Mandelbrot.xlsx'
}

k_values = [3, 4, 5, 6]
TICK_FONTSIZE = 16

for dist, path in files.items():
    for k in k_values:
        df = pd.read_excel(path, sheet_name=f'k{k}')
        df = df[df['R2'] >= 0]
        
        gc_content = pd.to_numeric(df['GC Content (%)'], errors='coerce').dropna()
        r2_values = df.loc[gc_content.index, 'R2'].astype(float)
        
        if gc_content.empty:
            print(f'No data for {dist} k={k} after filtering.')
            continue

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.minorticks_on()
        ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.7)
        ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)

        hb = ax.hexbin(
            gc_content,
            r2_values,
            gridsize=75,
            cmap='viridis', 
            mincnt=1,
            norm=LogNorm()
        )
        
        # Create the colorbar
        cbar = fig.colorbar(hb, ax=ax)
        cbar.ax.tick_params(labelsize=TICK_FONTSIZE)
        
        ax.set_xlabel('GC Content (%)', fontsize=16)
        ax.set_ylabel('R²', fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)

        plt.tight_layout()
        outname = f'R2_vs_GC_content_{dist.replace(" ","_")}_k{k}.png'
        plt.savefig(outname, dpi=300, bbox_inches='tight')
        print(f'Saved plot: {outname}')
        plt.close(fig)