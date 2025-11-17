#!/usr/bin/env python3

import sys
import json
import os
import gzip
from Bio import SeqIO

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATCG", "TAGC")
    return seq.translate(complement_map)[::-1]

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <bucket_id> <k> <scheduler.json>")
        sys.exit(1)

    bucket_id = sys.argv[1]
    k = int(sys.argv[2])
    scheduler_path = sys.argv[3]

    with open(scheduler_path, 'r') as f:
        scheduler = json.load(f)
    if bucket_id not in scheduler:
        print(f"Bucket ID '{bucket_id}' not found in scheduler.")
        sys.exit(1)

    output_dir = "/scratch/cpk5664/extra_3_tru_mers"
    os.makedirs(output_dir, exist_ok=True)

    for fasta_path in scheduler[bucket_id]:
        base_name = os.path.basename(fasta_path)
        output_filename = os.path.join(output_dir, f"{base_name}_kmers_{k}.txt")
        kmer_counts = {}

        if fasta_path.endswith(".gz"):
            handle = gzip.open(fasta_path, "rt")
        else:
            handle = open(fasta_path, "r")

        with handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_str = str(record.seq).upper()
                seq_length = len(seq_str)
                if seq_length < k:
                    continue
                for i in range(seq_length - k + 1):
                    kmer = seq_str[i:i+k]
                    if any(base not in "ATCG" for base in kmer):
                        continue
                    rc_kmer = reverse_complement(kmer)
                    canonical_kmer = kmer if kmer < rc_kmer else rc_kmer
                    kmer_counts[canonical_kmer] = kmer_counts.get(canonical_kmer, 0) + 1

        with open(output_filename, 'w') as out_f:
            out_f.write("#kmer\tcount\n")
            # Write kmers sorted by count in descending order
            for kmer, count in sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True):
                out_f.write(f"{kmer}\t{count}\n")

        print(f"Processed {fasta_path}, results written to: {output_filename}")
        kmer_counts.clear()

if __name__ == "__main__":
    main()