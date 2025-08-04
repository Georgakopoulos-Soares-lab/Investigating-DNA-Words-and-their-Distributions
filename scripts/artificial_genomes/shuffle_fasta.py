#!/usr/bin/env python3
import os
import sys
import json
import gzip
import random
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def shuffle_sequence(seq):
    """Return a new sequence string with nucleotides shuffled."""
    seq_list = list(str(seq))
    random.shuffle(seq_list)
    return "".join(seq_list)

def process_file(inpath, outpath):
    """Read inpath .fna.gz, shuffle each record, and write to outpath .fna.gz."""
    with gzip.open(inpath, "rt") as hin, gzip.open(outpath, "wt") as hout:
        records = []
        for rec in SeqIO.parse(hin, "fasta"):
            shuffled = shuffle_sequence(rec.seq)
            records.append(SeqRecord(Seq(shuffled),
                                     id=rec.id,
                                     description=rec.description))
        SeqIO.write(records, hout, "fasta")

def make_shuffled_name(filename: str) -> str:
    suffix = ".fna.gz"
    if filename.endswith(suffix):
        base = filename[:-len(suffix)]
        return f"{base}_shuffled{suffix}"
    else:
        return f"{filename}_shuffled"

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python shuffle.py <scheduler.json> <outdir>")

    scheduler_path = sys.argv[1]
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        with open(scheduler_path) as f:
            sched = json.load(f)
    except Exception as e:
        sys.exit(f"Error reading scheduler file: {e}")

    for bucket_id, files in sched.items():
        print(f"Processing bucket {bucket_id} with {len(files)} files…")
        for filepath in files:
            infile = Path(filepath)
            if not infile.exists():
                print(f"  [!] WARNING: {infile} not found, skipping.")
                continue

            new_name = make_shuffled_name(infile.name)
            outfile = outdir / new_name

            print(f"  → Shuffling {infile.name} → {new_name}")
            process_file(str(infile), str(outfile))

if __name__ == "__main__":
    main()