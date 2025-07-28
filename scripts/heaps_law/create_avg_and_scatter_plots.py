import os
import csv
import math
import re
import matplotlib.pyplot as plt

PLOT_SCATTER = True
k_mer_length = 12

MERGED_FILE = '/scratch/cpk5664/heap_results.txt'
TAX_CSV     = '/storage/group/izg5139/default/xaris/taxonomies.csv'

DOMAINS = ['viral', 'bacteria', 'archaea','eukaryote']
COLORS  = {
    'archaea':   'cyan',  
    'bacteria':  'red',
    'viral':     'purple',
    'eukaryote': 'green',
}

EUK_LIST = ['protozoa', 'vertebrate', 'vertebrate(other)', 'fungi', 'plant', 'invertebrate']
KS_ALL = list(range(6, 16))

def extract_ids(fname):
    base  = os.path.basename(fname)
    parts = base.split('_')
    acc   = "_".join(parts[:2])
    asm   = parts[2]
    return acc, asm

def detect_domain(raw):
    val = raw.strip().lower()
    if 'archaea' in val:   return 'archaea'
    if 'bacteria' in val:  return 'bacteria'
    if 'viral' in val:     return 'viral'
    for e in EUK_LIST:
        if e in val:        return 'eukaryote'
    return None

def build_needed_keys(path):
    needed = set()
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            fname = line.split('\t',1)[0]
            needed.add(extract_ids(fname))
    return needed

def build_tax_map(csv_path, needed_keys):
    taxmap = {}
    with open(csv_path, newline='') as csvf:
        reader = csv.DictReader(csvf)
        for row in reader:
            key = (row['Accession (GCF)'], row['Assembly Name'])
            if key in needed_keys:
                dom = detect_domain(row['Genome Type and Domain'])
                if dom:
                    taxmap[key] = dom
    return taxmap

def main():
    needed = build_needed_keys(MERGED_FILE)
    print(f"Total unique assemblies in merged file: {len(needed)}")

    taxmap = build_tax_map(TAX_CSV, needed)
    matched   = len(taxmap)
    unmatched = len(needed) - matched
    print(f"→ {matched} assemblies matched a domain")
    print(f"→ {unmatched} assemblies left unmatched and will be skipped")

    domain_counts = {dom: 0 for dom in DOMAINS}
    for dom in taxmap.values():
        if dom in domain_counts:
            domain_counts[dom] += 1
    for dom in DOMAINS:
        print(f"{domain_counts[dom]} files successfully matched into {dom}")

    stats = {
        dom: {
            k: {'sum_K':0.0, 'sum_K2':0.0, 'sum_beta':0.0, 'sum_beta2':0.0, 'count':0}
            for k in KS_ALL
        }
        for dom in DOMAINS
    }
    data_points = []  # list of (k, K, β, domain)

    k_re = re.compile(r'K:([^ \t]+)')
    b_re = re.compile(r'β:([^ \t]+)')
    with open(MERGED_FILE) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            acc_asm = extract_ids(parts[0])
            dom = taxmap.get(acc_asm)
            if not dom:
                continue
            for k in KS_ALL:
                idx = k - 5
                if idx < 0 or idx >= len(parts):
                    continue
                field = parts[idx]
                mK = k_re.search(field)
                mB = b_re.search(field)
                if not (mK and mB):
                    continue
                try:
                    Kv = float(mK.group(1))
                    Bv = float(mB.group(1))
                except ValueError:
                    continue
                rec = stats[dom][k]
                rec['sum_K']     += Kv
                rec['sum_K2']    += Kv*Kv
                rec['sum_beta']  += Bv
                rec['sum_beta2'] += Bv*Bv
                rec['count']     += 1
                data_points.append((k, Kv, Bv, dom))

    for dom in DOMAINS:
        meanK, stdK = [], []
        meanB, stdB = [], []
        for k in KS_ALL:
            rec = stats[dom][k]
            cnt = rec['count']
            if cnt > 0:
                mK = rec['sum_K'] / cnt
                vK = rec['sum_K2']/cnt - mK*mK
                mB = rec['sum_beta'] / cnt
                vB = rec['sum_beta2']/cnt - mB*mB
                meanK.append(mK)
                stdK.append(math.sqrt(vK) if vK>0 else 0.0)
                meanB.append(mB)
                stdB.append(math.sqrt(vB) if vB>0 else 0.0)
            else:
                meanK.append(math.nan); stdK.append(0.0)
                meanB.append(math.nan); stdB.append(0.0)

        color = COLORS[dom]
        fig, (axK, axB) = plt.subplots(1, 2, figsize=(12,4), sharex=True)

        # Avg K with ±1σ
        axK.plot(KS_ALL, meanK, marker='o', color=color)
        lowerK = [max(m-s, 0) for m,s in zip(meanK,stdK)]
        upperK = [m+s for m,s in zip(meanK,stdK)]
        axK.fill_between(KS_ALL, lowerK, upperK, alpha=0.2, color=color)
        axK.set_xlabel('k-mer length', fontsize=16)
        axK.set_ylabel('Average K', fontsize=16)
        axK.set_xticks(KS_ALL)
        axK.set_ylim(bottom=0)
        axK.tick_params(axis='both', labelsize=14)
        axK.grid(True, linestyle='--', alpha=0.5)

        # Avg β with ±1σ
        axB.plot(KS_ALL, meanB, marker='o', color=color)
        lowerB = [max(m-s, 0) for m,s in zip(meanB,stdB)]
        upperB = [m+s for m,s in zip(meanB,stdB)]
        axB.fill_between(KS_ALL, lowerB, upperB, alpha=0.2, color=color)
        axB.set_xlabel('k-mer length', fontsize=16)
        axB.set_ylabel('Average β', fontsize=16)
        axB.set_xticks(KS_ALL)
        axB.set_ylim(bottom=0)
        axB.tick_params(axis='both', labelsize=14)
        axB.grid(True, linestyle='--', alpha=0.5)

        fig.tight_layout()
        out_png = f'{dom}_avgK_beta.png'
        fig.savefig(out_png)
        plt.close(fig)
        print(f"Wrote plot for {dom}: {out_png}")

    if PLOT_SCATTER:
        pts = data_points
        if k_mer_length is not None:
            pts = [pt for pt in pts if pt[0] == k_mer_length]
        if not pts:
            print("No data points collected for scatter plot.")
        else:
            fig, ax = plt.subplots(figsize=(6,6))
            for dom in DOMAINS:
                dom_pts = [(K, B) for (k, K, B, d) in pts if d == dom]
                if dom_pts:
                    Ks = [pt[0] for pt in dom_pts]
                    Bs = [pt[1] for pt in dom_pts]
                    MARKERS = {'archaea':'o','bacteria':'^','viral':'s','eukaryote':'D'}
                    ax.scatter(Bs, Ks, s=20, alpha=0.4, label=dom, color=COLORS[dom])
            ax.grid(True, linestyle='--', alpha=0.3)
            ax.set_xlabel('β', fontsize=16)
            ax.set_yscale('log')
            ax.set_ylabel('K', fontsize=16)
            ax.tick_params(axis='both', labelsize=14)
            ax.legend()
            fig.tight_layout()
            out_scatter = 'scatter_K_beta_12.png'
            fig.savefig(out_scatter)
            plt.close(fig)
            print(f"Wrote scatter plot: {out_scatter}")

    print("All plots generated.")

if __name__ == '__main__':
    main()