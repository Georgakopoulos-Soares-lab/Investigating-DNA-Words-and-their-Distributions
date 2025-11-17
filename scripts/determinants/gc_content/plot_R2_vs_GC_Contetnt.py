#!/usr/bin/env python3
"""
Generates eight hexbin plots (one per distribution per k for k=3,4,5,6) of R² vs. GC Content (%).
Each plot shows data density and:
  - A colorbar indicating point density (log scale).
  - Major & minor grid lines.

Additionally, computes Spearman correlations between GC Content (%) and R²
within GC-content quartiles for each (distribution, k) and prints them.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import spearmanr

files = {
    'Truncated': '/work/10906/hariskil/vista/zipf/xaris/Investigating_Dna_Words/scripts/model_fits/Results_truncated.xlsx',
    'Zipf-Mandelbrot': '/work/10906/hariskil/vista/zipf/xaris/Investigating_Dna_Words/scripts/model_fits/Results_Zipf_Mandelbrot.xlsx'
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

        # ---------- New: Spearman per GC-content quartile ----------
        # Compute quartile boundaries
        q1, q2, q3 = np.quantile(gc_content, [0.25, 0.5, 0.75])
        bounds = [gc_content.min(), q1, q2, q3, gc_content.max()]

        print(f'\n[Spearman GC vs R²] {dist}, k={k}')
        for q in range(1, 5):
            lo = bounds[q - 1]
            hi = bounds[q]
            # Last bin includes the upper bound, others are half-open
            if q < 4:
                mask = (gc_content >= lo) & (gc_content < hi)
            else:
                mask = (gc_content >= lo) & (gc_content <= hi)

            gc_q = gc_content[mask]
            r2_q = r2_values[mask]

            if len(gc_q) < 2 or gc_q.nunique() < 2 or r2_q.nunique() < 2:
                print(f'  Q{q}: N={len(gc_q)} -> not enough variation for Spearman.')
                continue

            rho, pval = spearmanr(gc_q, r2_q)
            print(
                f'  Q{q}: [{lo:.2f}, {hi:.2f}] N={len(gc_q)} '
                f'rho={rho:.3f}, p={pval:.3e}'
            )

        # fig, ax = plt.subplots(figsize=(12, 8))
        # ax.minorticks_on()
        # ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.7)
        # ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)

        # hb = ax.hexbin(
        #     gc_content,
        #     r2_values,
        #     gridsize=75,
        #     cmap='viridis',
        #     mincnt=1,
        #     norm=LogNorm()
        # )
        
        # # Create the colorbar
        # cbar = fig.colorbar(hb, ax=ax)
        # cbar.ax.tick_params(labelsize=TICK_FONTSIZE)
        
        # ax.set_xlabel('GC Content (%)', fontsize=16)
        # ax.set_ylabel('R²', fontsize=16)
        # ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)

        # plt.tight_layout()
        # outname = f'R2_vs_GC_content_{dist.replace(" ","_")}_k{k}.png'
        # plt.savefig(outname, dpi=300, bbox_inches='tight')
        # print(f'Saved plot: {outname}')
        # plt.close(fig)