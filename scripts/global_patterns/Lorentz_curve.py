import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ks = [3, 4, 5, 6, 7, 8]
base_path = '/storage/group/izg5139/default/xaris'    
file_template = 'eukaryote_avg_{0}mers.txt'     
output_file = 'lorenz_curve_eukaryote.png'      

custom_colors = {
    3: '#1f77b4',  # blue
    4: '#ff7f0e',  # orange
    5: '#2ca02c',  # green
    6: '#d62728',  # red
    7: '#9467bd',  # purple
    8: '#8c564b'   # brown
}

plt.figure(figsize=(8, 6))
plt.plot([0, 1], [0, 1], linestyle='--', color='gray')

for k in ks:
    filepath = os.path.join(base_path, file_template.format(k))

    df = pd.read_csv(filepath,
                     sep=r'\s+',
                     header=None,
                     names=['avg', 'std'])
    
    counts = np.sort(df['avg'].values)
    probs = counts / counts.sum()
    
    cum_kmers = np.arange(1, len(probs) + 1) / len(probs)
    cum_counts = np.cumsum(probs)
    
    plt.plot(cum_kmers, cum_counts, color=custom_colors[k])

plt.xlabel('Cumulative share of k-mers', fontsize=16)
plt.ylabel('Cumulative share of total counts', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.tight_layout()

plt.savefig(output_file, dpi=300, bbox_inches='tight')

print(f"Lorenz curve plot saved as: {output_file}")