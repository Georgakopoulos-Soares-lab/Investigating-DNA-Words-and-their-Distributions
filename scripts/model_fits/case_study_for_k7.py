import pandas as pd
import matplotlib.pyplot as plt

TRUNC_FILE       = 'Results_truncated.xlsx'
ZIPF_FILE        = 'Results_Zipf_Mandelbrot.xlsx'
OUTPUT_COUNT_BAR = 'Organisms_with_positive_negative_R2_k7.png'
OUTPUT_AVG_R2    = 'Avg_R2_k7.png'
OUTPUT_TRUNC_POS = 'Trunc_positive_params.png'
OUTPUT_ZIPF_POS  = 'Zipf_Mandelbrot_positive_params.png'

df_trunc = pd.read_excel(TRUNC_FILE, sheet_name='k7')
df_zipf  = pd.read_excel(ZIPF_FILE,  sheet_name='k7')

count_data = {
    'Truncated Power Law': {
        'Positive': (df_trunc['R2'] > 0).sum(),
        'Negative': (df_trunc['R2'] < 0).sum(),
    },
    'Zipf Mandelbrot': {
        'Positive': (df_zipf['R2'] > 0).sum(),
        'Negative': (df_zipf['R2'] < 0).sum(),
    }
}
counts_df = pd.DataFrame(count_data).T

avg_r2 = {
    'Truncated Power Law': df_trunc.loc[df_trunc['R2'] > 0, 'R2'].mean(),
    'Zipf Mandelbrot':     df_zipf.loc[df_zipf['R2'] > 0, 'R2'].mean(),
}
avg_r2_series = pd.Series(avg_r2)

fig1, ax1 = plt.subplots(figsize=(6, 5))
counts_df.plot(kind='bar', ax=ax1)
ax1.set_xlabel('') 
ax1.set_xticklabels(counts_df.index, rotation=0)
ax1.set_ylabel('Number of Organisms')
ax1.xaxis.label.set_fontsize(16)
ax1.yaxis.label.set_fontsize(16)
ax1.tick_params(axis='both', labelsize=16)
ax1.legend(title='R² Sign', loc='upper left')

def annotate_bars(ax):
    for p in ax.patches:
        ax.annotate(
            format(int(p.get_height()), ','),
            (p.get_x() + p.get_width() / 2., p.get_height()),
            ha='center', va='bottom', fontsize=14
        )

annotate_bars(ax1)
plt.tight_layout()
fig1.savefig(OUTPUT_COUNT_BAR, dpi=300)
plt.close(fig1)

fig2, ax2 = plt.subplots(figsize=(6, 5))
avg_r2_series.plot(kind='bar', ax=ax2, color=['C0', 'C1'])
ax2.set_xlabel('')
ax2.set_xticklabels(avg_r2_series.index, rotation=0)
ax2.set_ylabel('Average R²')
ax2.xaxis.label.set_fontsize(16)
ax2.yaxis.label.set_fontsize(16)
ax2.tick_params(axis='both', labelsize=16)

for p in ax2.patches:
    ax2.annotate(
        f"{p.get_height():.3f}",
        (p.get_x() + p.get_width() / 2., p.get_height()),
        ha='center', va='bottom', fontsize=14
    )

plt.tight_layout()
fig2.savefig(OUTPUT_AVG_R2, dpi=300)
plt.close(fig2)

trunc_pos_df = df_trunc[df_trunc['R2'] > 0][['alpha', 'lambda', 'scale']]
zipf_pos_df  = df_zipf [df_zipf['R2']  > 0][['alpha', 'beta',  'scale']]
mean_trunc   = trunc_pos_df.mean()
std_trunc    = trunc_pos_df.std()
mean_zipf    = zipf_pos_df.mean()
std_zipf     = zipf_pos_df.std()

def save_table_trunc(mean_series, std_series, filename, title):
    rows = [f"{mean_series[idx]:.3f} ± {std_series[idx]:.3f}" for idx in mean_series.index]
    df = pd.DataFrame(rows, index=mean_series.index, columns=['Mean ± Std'])
    fig, ax = plt.subplots(figsize=(4, len(df)*0.5 + 1))
    ax.axis('off')
    ax.table(
        cellText=df.values,
        colLabels=df.columns,
        rowLabels=df.index,
        cellLoc='center',
        loc='center'
    )
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close(fig)

def save_table_zipf(mean_series, std_series, filename, title):
    rows = []
    for idx in mean_series.index:
        m, s = mean_series[idx], std_series[idx]
        if idx == 'scale':
            rows.append(f"{m:.3e} ± {s:.3e}")
        else:
            rows.append(f"{m:.3f} ± {s:.3f}")
    df = pd.DataFrame(rows, index=mean_series.index, columns=['Mean ± Std'])
    fig, ax = plt.subplots(figsize=(4, len(df)*0.5 + 1))
    ax.axis('off')
    ax.table(
        cellText=df.values,
        colLabels=df.columns,
        rowLabels=df.index,
        cellLoc='center',
        loc='center'
    )
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close(fig)

save_table_trunc(
    mean_trunc, std_trunc,
    OUTPUT_TRUNC_POS,
    ''
)
save_table_zipf(
    mean_zipf, std_zipf,
    OUTPUT_ZIPF_POS,
    ''
)