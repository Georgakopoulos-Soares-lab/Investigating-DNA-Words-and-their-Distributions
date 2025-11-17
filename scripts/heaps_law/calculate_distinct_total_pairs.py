import json
import argparse
import gzip
import os
from pathlib import Path
from Bio import SeqIO
import sys

_COMPLEMENT = str.maketrans({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'})

def extract_genome_name(file_path):
    return Path(file_path).name

def open_fasta(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    else:
        return open(path, 'r')

def count_total_windows(path, k):
    total = 0
    with open_fasta(path) as handle:
        for rec in SeqIO.parse(handle, 'fasta'):
            seq = str(rec.seq).upper()
            run_len = 0
            for c in seq:
                if c in 'ATCG':
                    run_len += 1
                else:
                    if run_len >= k:
                        total += run_len - k + 1
                    run_len = 0
            if run_len >= k:
                total += run_len - k + 1
    return total

def canonical_kmer(kmer):
    """Return lexicographically smaller of kmer vs its reverse-complement."""
    rc = kmer.translate(_COMPLEMENT)[::-1]
    return rc if rc < kmer else kmer

def generate_kmers_progressive_streaming(path, k, output_file, genome_name):
    # figure out how many windows there are
    tot_windows = count_total_windows(path, k)
    if tot_windows < 1:
        print(f"Error: No valid {k}-mers in {genome_name}", file=sys.stderr)
        return

    seen = set()
    total_count = 0
    done = set()              # which %'s we’ve already recorded
    pending = set(range(1, 101))

    with open(output_file, 'w') as out:
        # slide over each valid run
        proc = 0
        with open_fasta(path) as handle:
            for rec in SeqIO.parse(handle, 'fasta'):
                seq = str(rec.seq).upper()
                n = len(seq)
                i = 0
                while i <= n - k:
                    if seq[i] not in 'ATCG':
                        i += 1
                        continue
                    # find the end of this ATCG run
                    j = i
                    while j < n and seq[j] in 'ATCG':
                        j += 1
                    run_len = j - i
                    if run_len >= k:
                        # slide across this clean run
                        for r in range(i, j - k + 1):
                            kmer = seq[r:r+k]
                            can = canonical_kmer(kmer)
                            seen.add(can)
                            total_count += 1
                            proc += 1

                            pct = int(proc * 100 / tot_windows)
                            if pct in pending:
                                out.write(f"{len(seen)}:{total_count}\n")
                                pending.remove(pct)
                    i = j + 1

        # ensure we wrote 100% if it never hit exactly
        if 100 in pending:
            out.write(f"{len(seen)}:{total_count}\n")

def process_genome_file(path, k, output_dir):
    name = extract_genome_name(path)
    outp = os.path.join(output_dir, f"{name}.txt")
    if os.path.exists(outp):
        print(f"Skipping {name}: output exists", file=sys.stderr)
        return
    try:
        generate_kmers_progressive_streaming(path, k, outp, name)
        print(f"Done {name}")
    except Exception as e:
        print(f"Error {name}: {e}", file=sys.stderr)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('scheduler_file')
    p.add_argument('k', type=int)
    p.add_argument('output_dir')
    p.add_argument('--bucket-id', type=str, default=None,
                  help="Bucket id to process; if omitted, uses SLURM_PROCID.")
    args = p.parse_args()

    if not os.path.exists(args.scheduler_file):
        sys.exit(f"Scheduler not found: {args.scheduler_file}")
    if args.k < 1:
        sys.exit(f"k must be ≥1, got {args.k}")

    if args.bucket_id is not None:
        bucket_id = args.bucket_id
    else:
        bucket_id = os.environ.get("SLURM_PROCID")
        if bucket_id is None:
            sys.exit("No --bucket-id provided and SLURM_PROCID is not set.")

    os.makedirs(args.output_dir, exist_ok=True)
    with open(args.scheduler_file) as f:
        sched = json.load(f)

    if bucket_id not in sched:
        print(f"Bucket '{bucket_id}' not in scheduler; nothing to do.", file=sys.stderr)
        return

    for fasta in sched[bucket_id]:
        if not os.path.exists(fasta):
            print(f"Missing file, skipping: {fasta}", file=sys.stderr)
            continue
        process_genome_file(fasta, args.k, args.output_dir)

if __name__ == "__main__":
    main()