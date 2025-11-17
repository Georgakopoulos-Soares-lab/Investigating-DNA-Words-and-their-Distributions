#!/usr/bin/env python3
"""
Computes the average AIC and standard deviation for k = 3,4,5,6 for each distribution,
prints the count of values used, how many were negative, and plots a grouped bar chart
with a broken y-axis so k=6 doesn't dominate.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

files = {
    'Truncated': 'Results_truncated.xlsx',
    'Zipf-Mandelbrot': 'Results_Zipf_Mandelbrot.xlsx'
}
output_plot = 'Average_AIC.png'
s_values = [3, 4, 5, 6]

avg_aic = {name: [] for name in files}
std_aic = {name: [] for name in files}
counts = {name: [] for name in files}

for name, path in files.items():
    for k in s_values:
        sheet = f'k{k}'
        try:
            df = pd.read_excel(path, sheet_name=sheet)
        except Exception as e:
            print(f"Error reading {sheet} from {path}: {e}")
            avg_aic[name].append(np.nan)
            std_aic[name].append(np.nan)
            counts[name].append(0)
            continue

        aic_vals = df['AIC'].dropna().astype(float)
        total_count = len(aic_vals)
        neg_count = (aic_vals < 0).sum()

        used_count = total_count
        mean_aic = aic_vals.mean() if used_count > 0 else np.nan
        std = aic_vals.std(ddof=1) if used_count > 0 else np.nan

        counts[name].append(used_count)
        avg_aic[name].append(mean_aic)
        std_aic[name].append(std)

        mean_str = f"{mean_aic:.4f}" if pd.notna(mean_aic) else "nan"
        std_str = f"{std:.4f}" if pd.notna(std) else "nan"
        print(
            f"{name}, k={k}: {used_count} values used ({neg_count} negative), "
            f"avg AIC = {mean_str}, std = {std_str}"
        )

#Broken y-axis so k=6 doesn't crush the scale
x = np.arange(len(s_values))
width = 0.35

non_k6_idx = [i for i, k in enumerate(s_values) if k != 6]
k6_idx = s_values.index(6)

def _lims_from(indices):
    uppers, lowers = [], []
    for name in files:
        a = np.array(avg_aic[name], dtype=float)
        s = np.array(std_aic[name], dtype=float)
        s = np.where(np.isfinite(s), s, 0.0)
        a_sel = a[indices]
        s_sel = s[indices]
        upp = a_sel + s_sel
        low = a_sel - s_sel
        uppers.extend(list(upp[np.isfinite(upp)]))
        lowers.extend(list(low[np.isfinite(low)]))
    if not uppers or not lowers:
        return (-1, 1)
    top = max(uppers)
    bottom = min(lowers)
    pad = 0.08 * max(1.0, abs(top - bottom))
    return (bottom - pad, top + pad)

ylim_top = _lims_from(non_k6_idx)     
ylim_bot = _lims_from([k6_idx])      

fig, (ax_top, ax_bot) = plt.subplots(
    2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 2]}
)

def draw_bars(ax):
    for i, name in enumerate(files):
        offsets = x + (i - 0.5) * width
        ax.bar(
            offsets, avg_aic[name], width,
            yerr=std_aic[name], capsize=5  
        )
    ax.grid(axis='y', linestyle='--', alpha=0.7)

draw_bars(ax_top)
draw_bars(ax_bot)

ax_top.set_ylim(ylim_top)
ax_bot.set_ylim(ylim_bot)

ax_top.spines['bottom'].set_visible(False)
ax_bot.spines['top'].set_visible(False)
ax_top.tick_params(labeltop=False)
ax_bot.xaxis.tick_bottom()

kwargs = dict(color='k', clip_on=False, linewidth=1)
ax_top.plot((-0.015, +0.015), (-0.02, +0.02), transform=ax_top.transAxes, **kwargs)
ax_top.plot((0.985, 1.015), (-0.02, +0.02), transform=ax_top.transAxes, **kwargs)
ax_bot.plot((-0.015, +0.015), (1 - 0.02, 1 + 0.02), transform=ax_bot.transAxes, **kwargs)
ax_bot.plot((0.985, 1.015), (1 - 0.02, 1 + 0.02), transform=ax_bot.transAxes, **kwargs)

ax_bot.set_xlabel('k-mer length', fontsize=16)
for ax in (ax_top, ax_bot):
    ax.tick_params(axis='y', labelsize=14)

try:
    fig.supylabel('Average AIC', fontsize=16)
except AttributeError:
    fig.text(0.02, 0.5, 'Average AIC', va='center', rotation='vertical', fontsize=16)

ax_bot.set_xticks(x)
ax_bot.set_xticklabels([str(k) for k in s_values], fontsize=14)

plt.tight_layout()
plt.savefig(output_plot, bbox_inches='tight')
print(f"Plot saved to {output_plot}")