#!/usr/bin/env python3
"""
Loads a txt file that contains all pre calculated genic percentages.
Generates eight hexbin plots (one per distribution per k for k=3,4,5,6) of R² vs. Genic Percentage (%).
The plot shows data density and:
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

# Path to the genic percentage data file
genic_file_path = '/storage/group/izg5139/default/xaris/genic_percentage.txt'

def load_genic_data(path):
    """
    Reads a text file with "filename : value" format.
    Strips suffixes and stores data in a dictionary for fast lookup.
    Returns a dictionary mapping base filename to genic percentage.
    """
    genic_map = {}
    try:
        with open(path, 'r') as f:
            for line in f:
                if ' : ' in line:
                    parts = line.split(' : ')
                    filename = parts[0].strip()
                    value = float(parts[1].strip())
                    
                    base_name = filename.replace('.gff.gz', '')
                    genic_map[base_name] = value
    except FileNotFoundError:
        print(f"Error: Genic data file not found at '{path}'")
        print("Please ensure the file exists and the path is correct.")
        exit()
        
    print(f"Successfully loaded {len(genic_map)} entries from genic data file.")
    return genic_map

genic_data_map = load_genic_data(genic_file_path)

for dist, path in files.items():
    for k in k_values:
        df = pd.read_excel(path, sheet_name=f'k{k}')
        df = df[df['R2'] >= 0].copy() 

        def get_base_name(excel_filename):
            return excel_filename.replace('.fna.gz', '')

        df['Genic Percentage'] = df['Filename'].apply(
            lambda x: genic_data_map.get(get_base_name(x), np.nan)
        )

        df.dropna(subset=['Genic Percentage', 'R2'], inplace=True)
        
        initial_count = len(df)
        df = df[(df['Genic Percentage'] >= 0) & (df['Genic Percentage'] <= 100)]
        final_count = len(df)
        
        print(f"--- For {dist} k={k}: Found {final_count} files with valid genic percentages (filtered from {initial_count}). ---")
        
        genic_percentage = df['Genic Percentage']
        r2_values = df['R2'].astype(float)

        if genic_percentage.empty:
            print(f'No data to plot for {dist} k={k}. Skipping.')
            continue

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.minorticks_on()
        ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.7)
        ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)

        hb = ax.hexbin(
            genic_percentage,
            r2_values,
            gridsize=75,
            cmap='viridis',
            mincnt=1,
            norm=LogNorm()
        )
        
        cbar = fig.colorbar(hb, ax=ax)
        cbar.ax.tick_params(labelsize=TICK_FONTSIZE)
        ax.set_xlabel('Genic Percentage (%)', fontsize=16)
        ax.set_ylabel('R²', fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)

        plt.tight_layout()
        
        outname = f'R2_vs_Genic_Percentage_{dist.replace(" ","_")}_k{k}.png'
        plt.savefig(outname, dpi=300, bbox_inches='tight')
        print(f'Saved plot: {outname}\n')
        plt.close(fig)