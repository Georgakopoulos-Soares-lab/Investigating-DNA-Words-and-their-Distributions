#!/usr/bin/env python3

import tarfile
import matplotlib.pyplot as plt

tarfiles = [
    "/scratch/cpk5664/heap_6mers.tar",
    "/scratch/cpk5664/heap_7mers.tar",
    "/scratch/cpk5664/heap_8mers.tar",
    "/scratch/cpk5664/heap_9mers.tar",
    "/scratch/cpk5664/heap_10mers.tar",
    "/scratch/cpk5664/heap_11mers.tar",
    "/scratch/cpk5664/heap_12mers.tar",
    "/scratch/cpk5664/heap_13mers.tar",
    "/scratch/cpk5664/heap_14mers.tar",
    "/scratch/cpk5664/heap_15mers.tar"
]

# E.coli: "GCA_000599625.1_ASM59962v1_genomic.fna.gz.txt"
# homosapiens : "GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz.txt"

target_basename = "GCA_000599625.1_ASM59962v1_genomic.fna.gz.txt"
output_png      = "E.coli_heap_all_k.png"

plt.figure(figsize=(8,6),constrained_layout=True)

for tarpath in tarfiles:
    try:
        k = int(tarpath.split('_')[1].replace('mers.tar',''))
    except Exception:
        k = tarpath

    with tarfile.open(tarpath, 'r') as tar:
        members = [m for m in tar.getmembers() if m.name.endswith(target_basename)]
        if not members:
            print(f" {target_basename!r} not found in {tarpath}, skipping.")
            continue

        f = tar.extractfile(members[0])
        if not f:
            print(f" Failed to extract {members[0].name} from {tarpath}, skipping.")
            continue

        lines = f.read().decode('utf-8').strip().splitlines()

    xs, ys = [], []
    for line in lines:
        if ':' not in line:
            continue
        start, end = line.split(':', 1)
        try:
            xs.append(int(end))
            ys.append(int(start))
        except ValueError:
            continue

    plt.plot(xs, ys, marker='o', linestyle='-', label=f'k={k}')

plt.xlabel('Total k-mers', fontsize=17)
plt.ylabel('Distinct k-mers', fontsize=17)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.grid(True)
plt.xscale('log')
plt.yscale('log')

plt.savefig(output_png, dpi=300)
print(f"Plot saved as {output_png}")