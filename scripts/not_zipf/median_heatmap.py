import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

kmer3_txt = '/storage/group/izg5139/default/xaris/not_zipf/not_zipf_3mers_final.txt'
kmer4_txt = '/storage/group/izg5139/default/xaris/not_zipf/not_zipf_4mers_final.txt'
kmer5_txt = '/storage/group/izg5139/default/xaris/not_zipf/not_zipf_5mers_final.txt'
kmer6_txt = '/storage/group/izg5139/default/xaris/not_zipf/not_zipf_6mers_final.txt'
kmer7_txt = '/storage/group/izg5139/default/xaris/not_zipf/not_zipf_7mers_final.txt'
kmer8_txt = '/storage/group/izg5139/default/xaris/not_zipf/not_zipf_8mers_final.txt'

txt_paths = {
    '3': kmer3_txt,
    '4': kmer4_txt,
    '5': kmer5_txt,
    '6': kmer6_txt,
    '7': kmer7_txt,
    '8': kmer8_txt
}

taxonomy_csv = '/storage/group/izg5139/default/xaris/taxonomies.csv'

tax_df = pd.read_csv(
    taxonomy_csv,
    usecols=['Accession (GCF)', 'Assembly Name', 'Genome Type and Domain'],
    dtype=str
)

tax_df = tax_df.set_index(['Accession (GCF)', 'Assembly Name'], drop=False)
tax_df.sort_index(inplace=True)

eukaryote_subtypes = {
    "protozoa",
    "vertebrate",
    "vertebrate(other)",
    "fungi",
    "plant",
    "invertebrate"
}

domains = ['viral', 'bacteria', 'archaea', 'eukaryote']
collected = { k: { dom: [] for dom in domains } for k in txt_paths.keys() }

pattern = re.compile(r'^(?P<acc>GCA_[^_]+|GCF_[^_]+)_(?P<asm>[^_]+)_')

for k, path in txt_paths.items():
    matched_count = 0
    unmatched_list = []

    with open(path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            try:
                fname, metrics = line.split(':', 1)
            except ValueError:
                unmatched_list.append(f"{fname if 'fname' in locals() else line}")
                continue

            # Extract R² from metrics
            r2_match = re.search(r'R2=([-\d\.eE\+]+)', metrics)
            if not r2_match:
                unmatched_list.append(fname)
                continue
            r2_val = float(r2_match.group(1))

            # Extract accession & assembly from filename
            m = pattern.match(fname)
            if not m:
                unmatched_list.append(fname)
                continue
            acc = m.group('acc')
            asm = m.group('asm')

            try:
                raw_tax_field = tax_df.loc[(acc, asm), 'Genome Type and Domain']
            except KeyError:
                unmatched_list.append(fname)
                continue

            if isinstance(raw_tax_field, pd.Series):
                raw_tax = raw_tax_field.iloc[0].strip().lower()
            else:
                raw_tax = raw_tax_field.strip().lower()

            if raw_tax in eukaryote_subtypes:
                domain = 'eukaryote'
            elif raw_tax == 'viral':
                domain = 'viral'
            elif raw_tax == 'bacteria':
                domain = 'bacteria'
            elif raw_tax == 'archaea':
                domain = 'archaea'
            else:
                unmatched_list.append(fname)
                continue

            collected[k][domain].append(r2_val)
            matched_count += 1

    print(f"\nSummary for file (k = {k}):\n  Path: {path}")
    print(f" Successfully matched taxonomy for {matched_count} filenames.")

median_r2 = pd.DataFrame(index=domains, columns=sorted(txt_paths.keys(), key=int), dtype=float)
for k in sorted(collected.keys(), key=int):
    for dom in domains:
        vals = collected[k][dom]
        median_r2.loc[dom, k] = np.nan if not vals else np.median(vals)
median_r2 = median_r2[['3','4','5','6','7','8']]

plt.figure(figsize=(8, 6))
im = plt.imshow(
    median_r2.values,
    aspect='auto',
    origin='lower',
    interpolation='nearest'
)

plt.xticks(
    ticks=np.arange(len(median_r2.columns)),
    labels=median_r2.columns.astype(int),
    fontsize=16
)

plt.yticks(
    ticks=np.arange(len(median_r2.index)),
    labels=median_r2.index.str.capitalize(),
    fontsize=16
)

cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize=16)

plt.xlabel('k-mer length', fontsize=16)
plt.ylabel('Domain', fontsize=16)

for i, dom in enumerate(median_r2.index):
    for j, k in enumerate(median_r2.columns):
        val = median_r2.loc[dom, k]
        if not np.isnan(val):
            text = f"{val:.2f}"
            plt.text(
                j, i, text,
                ha='center', va='center',
                fontsize=16,
                color='black'
            )

plt.tight_layout()
plt.savefig('median_heatmap', dpi=300)