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
    'Eukaryote': ['protozoa', 'vertebrate', 'vertebrate_other', 'fungi', 'plant', 'invertebrate']
}

colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
indices = np.arange(len(k_values))

def _fmt(x):
    try:
        if pd.isna(x):
            return "nan"
        return f"{float(x):.6g}"
    except Exception:
        return "nan"

for tax_name, tax_list in taxonomies.items():
    stats = {
        dist: {param: {'mean': [], 'std': []} for param in params}
        for dist, params in distribution_params.items()
    }
    extra = {
        dist: {param: {'ci_low': [], 'ci_high': []}
            for param in params}
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
                stats[dist][param]['std'].append(vals.std(ddof=1) if len(vals) > 1 else (0.0 if len(vals)==1 else np.nan))

                # CI
                n = int(vals.shape[0])
                mean = stats[dist][param]['mean'][-1]
                sd = stats[dist][param]['std'][-1]
                if n and n > 1 and not pd.isna(mean) and not pd.isna(sd):
                    half_width = 1.96 * (sd / np.sqrt(n))
                    ci_low = mean - half_width
                    ci_high = mean + half_width
                else:
                    ci_low = np.nan
                    ci_high = np.nan

                extra[dist][param]['ci_low'].append(ci_low)
                extra[dist][param]['ci_high'].append(ci_high)


    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 12), sharex=True)
    for row_idx, (param_trunc, param_zipf) in enumerate(
        zip(distribution_params['Truncated'], distribution_params['Zipf-Mandelbrot'])
    ):
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

    def write_table(dist_name: str):
        params = distribution_params[dist_name]
        fname = f"Table_{dist_name}_{tax_name}.txt".replace(' ', '_')
        with open(fname, "w", encoding="utf-8") as f:
            f.write(f"# Distribution: {dist_name}\n")
            f.write(f"# Taxonomy: {tax_name}\n")
            f.write("# Rows with R2 < 0 excluded.\n")
            f.write("# Columns: for each parameter, CI_low, CI_high, CI_width (95% normal approx).\n")

            cols = ["k"]
            for p in params:
                cols += [f"{p}_ci_low", f"{p}_ci_high", f"{p}_ci_width"]
            f.write("\t".join(cols) + "\n")

            for i, k in enumerate(k_values):
                row = [str(k)]
                for p in params:
                    ci_low = extra[dist_name][p]['ci_low'][i]
                    ci_high = extra[dist_name][p]['ci_high'][i]
                    ci_width = (ci_high - ci_low) if (not pd.isna(ci_low) and not pd.isna(ci_high)) else np.nan
                    row += [ _fmt(ci_low), _fmt(ci_high), _fmt(ci_width) ]
                f.write("\t".join(row) + "\n")
        print(f"Wrote table: {fname}")

    write_table('Truncated')
    write_table('Zipf-Mandelbrot')