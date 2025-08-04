#!/usr/bin/env python3
import sys
import os
import tarfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


USE_TAR = True

# When USE_TAR is True
TAR_PATH = '/scratch/cpk5664/sorted_5mers.tar'
TXT_MEMBER = 'sorted_5mers/GCA_002853715.1_ASM285371v1_genomic.fna.gz_kmers_5_sorted.txt' #E.coli

# When USE_TAR is False
PLAINTXT_PATH = '/path/to/counts.txt'

def extract_lines_from_tar(tar_path, member_path):
    try:
        with tarfile.open(tar_path, 'r') as tar:
            member = tar.getmember(member_path)
            f = tar.extractfile(member)
            if f is None:
                print(f"Error: Could not extract '{member_path}' from '{tar_path}'")
                sys.exit(1)
            content = f.read().decode('utf-8')
            return content.strip().split('\n')
    except FileNotFoundError:
        print(f"Error: Tar file '{tar_path}' not found.")
        sys.exit(1)
    except KeyError:
        print(f"Error: Member '{member_path}' not found in '{tar_path}'.")
        sys.exit(1)


def main():
    # Load lines either from tar or from a plain file
    if USE_TAR:
        lines = extract_lines_from_tar(TAR_PATH, TXT_MEMBER)
        base_name = os.path.splitext(os.path.basename(TXT_MEMBER))[0]
    else:
        try:
            with open(PLAINTXT_PATH, 'r') as f:
                lines = f.read().strip().split('\n')
        except FileNotFoundError:
            print(f"Error: File '{PLAINTXT_PATH}' not found.")
            sys.exit(1)
        base_name = os.path.splitext(os.path.basename(PLAINTXT_PATH))[0]

    # Parse k-mer counts
    kmers = []
    for line in lines:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        kmer, cnt_str = parts[0], parts[1]
        try:
            count = int(cnt_str)
        except ValueError:
            continue
        kmers.append((kmer, count))
    
    if not kmers:
        print("No data provided in the input file.")
        sys.exit(0)
    
    max_rank = min(2000, len(kmers))
    
    # Extract counts and compute standard Zipf (based on rank-1 count)
    rank1_count = kmers[0][1]
    ranks = list(range(1, max_rank + 1))
    real_counts = [count for (_, count) in kmers[:max_rank]]
    theory_counts = [rank1_count / r for r in ranks]

    fig, ax = plt.subplots(figsize=(7, 5))

    # Log-Log Zipf plot
    ax.loglog(ranks, real_counts, marker='o', linestyle='none', label='Real Data')
    ax.loglog(ranks, theory_counts, marker='x', linestyle='-', label='Standard Zipf')
    ax.set_xlabel('Rank', fontsize =16)
    ax.set_ylabel('Count', fontsize =16)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)

    ax.legend(loc='best')

    # # Linear Zipf plot
    # ax2.plot(ranks, real_counts, marker='o', linestyle='none', label='Real Data')
    # ax2.plot(ranks, theory_counts, marker='x', linestyle='-', label='Standard Zipf')
    # ax2.set_xlabel('Rank')
    # ax2.set_ylabel('Count')
    # ax2.set_title('Real vs Theoretical Zipf Plot for Escherichia coli (k=6) (Linear scale)')
    # ax2.legend(loc='best')

    plt.tight_layout()
    plt.savefig("E.coli_5_not_zipf.png", dpi=150)
    plt.close()

if __name__ == "__main__":
    main()