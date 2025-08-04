#!/usr/bin/env python3
import os
import re
import tarfile
import numpy as np
import pandas as pd

def collapse_domain(dom):
    """
    Collapse the full CSV 'Genome Type and Domain' into one of:
      'archaea', 'bacteria', 'viral', 'eukaryote'
    """
    dom = dom.lower()
    if dom in ('archaea', 'bacteria', 'viral'):
        return dom
    euk_groups = {
        'fungi',
        'invertebrate',
        'plant',
        'protozoa',
        'vertebrate_mammalian',
        'vertebrate_other',
    }
    if dom in euk_groups:
        return 'eukaryote'
    raise KeyError(f"Unknown domain '{dom}'")


def update_accumulator(acc, data):
    """
    Update running count, sum, and sum of squares for online mean/std.
    """
    count = acc['count']
    sum1 = acc['sum1']
    sum2 = acc['sum2']
    M = data.shape[0]

    if sum1 is None:
        sum1 = np.zeros(M, dtype=float)
        sum2 = np.zeros(M, dtype=float)
    elif sum1.size < M:
        new_sum1 = np.zeros(M, dtype=float)
        new_sum1[:sum1.size] = sum1
        sum1 = new_sum1
        new_sum2 = np.zeros(M, dtype=float)
        new_sum2[:sum2.size] = sum2
        sum2 = new_sum2

    sum1[:M] += data
    sum2[:M] += data**2
    count += 1

    acc['count'] = count
    acc['sum1'] = sum1
    acc['sum2'] = sum2

def main():
    tar_dir     = "/scratch/cpk5664"
    mapping_csv = "/storage/group/izg5139/default/xaris/taxonomies.csv"
    flat_dir    = "/storage/group/izg5139/default/xaris/sorted_4mers"
    output_dir  = "/storage/group/izg5139/default/xaris"
    ks          = [3,4,5,6,7,8]

    # load taxonomy mappings
    df = pd.read_csv(mapping_csv, low_memory=False)
    df.columns = df.columns.str.strip()
    acc2dom = df.set_index('Accession (GCF)')['Genome Type and Domain'].to_dict()
    asm2dom = df.set_index('Assembly Name')['Genome Type and Domain'].to_dict()

    # regex to parse filenames
    pat = re.compile(r'(GCA_\d+\.\d+)_([^_/]+)_genomic')
    os.makedirs(output_dir, exist_ok=True)

    domains = ('archaea', 'bacteria', 'viral', 'eukaryote')
    for k in ks:
        accumulators = {
            dom: {'count': 0, 'sum1': None, 'sum2': None}
            for dom in domains
        }

        # scan tarfiles
        for fn in os.listdir(tar_dir):
            if not (fn.endswith('.tar') or fn.endswith('.tar.gz')):
                continue
            if f"sorted_{k}mers" not in fn:
                continue
            tarpath = os.path.join(tar_dir, fn)
            with tarfile.open(tarpath, 'r:*') as tar:
                for member in tar.getmembers():
                    if not member.name.endswith(f"_kmers_{k}_sorted.txt"):
                        continue
                    m = pat.search(member.name)
                    if not m:
                        print("SKIP (no regex match):", member.name)
                        continue
                    acc, asm = m.group(1), m.group(2)
                    dom_full = acc2dom.get(acc) or asm2dom.get(asm)
                    if dom_full is None:
                        print("WARN: no taxonomy for", acc, asm)
                        continue
                    dom0 = collapse_domain(dom_full)
                    with tar.extractfile(member) as fh:
                        data = np.loadtxt(fh, usecols=1)
                    update_accumulator(accumulators[dom0], data)

        # scan normal dir
        for fn in os.listdir(flat_dir):
            if not fn.endswith(f"_kmers_{k}_sorted.txt"):
                continue
            m = pat.search(fn)
            if not m:
                print("SKIP (no regex match):", fn)
                continue
            acc, asm = m.group(1), m.group(2)
            dom_full = acc2dom.get(acc) or asm2dom.get(asm)
            if dom_full is None:
                print("WARN: no taxonomy for", acc, asm)
                continue
            dom0 = collapse_domain(dom_full)
            path_txt = os.path.join(flat_dir, fn)
            data = np.loadtxt(path_txt, usecols=1)
            update_accumulator(accumulators[dom0], data)

        print(f"k={k} summary:")
        for dom0, acc in accumulators.items():
            print(f"  {dom0}: {acc['count']} files processed")

        for dom0, acc in accumulators.items():
            count = acc['count']
            if count == 0:
                print(f"WARNING: no data for k={k}, domain={dom0}")
                continue
            sum1, sum2 = acc['sum1'], acc['sum2']
            mean = sum1 / count
            var = (sum2 / count) - mean**2
            std = np.sqrt(np.maximum(var, 0))
            out = np.vstack([mean, std]).T

            outpath = os.path.join(output_dir, f"{dom0}_avg_{k}mers.txt")
            np.savetxt(outpath, out, fmt="%.6f", header="mean\tstd", comments="")
            print(f"WROTE: {outpath}\n")

if __name__ == "__main__":
    main()
