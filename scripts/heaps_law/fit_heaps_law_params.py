#!/usr/bin/env python3
import json
import os
import argparse
import tarfile
import numpy as np
from scipy.optimize import curve_fit

TAR_FILE = "/scratch/cpk5664/heap_15mers.tar"
OUTPUT_DIR = "/scratch/cpk5664"

def heaps_model(n, K, beta):
    return K * n**beta

def fit_heaps_scipy(V, N):
    popt, _ = curve_fit(heaps_model, N, V, p0=[1.0, 0.5], maxfev=10_000)
    return popt 

def main():
    parser = argparse.ArgumentParser(
        description="Fit Heaps' law (V = K n^β) for each bucket"
    )
    parser.add_argument("scheduler", help="path to scheduler")
    parser.add_argument("bucket", help="bucket ID to process")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    with open(args.scheduler) as f:
        sched = json.load(f)

    if args.bucket not in sched:
        print(f"Error: bucket '{args.bucket}' not found in scheduler.")
        return

    files = sched[args.bucket]
    
    if not os.path.exists(TAR_FILE):
        print(f"Error: '{TAR_FILE}' not found")
        return

    print(f"Processing bucket {args.bucket} with {len(files)} files")
    print(f"Opening tar file: {TAR_FILE}")
    
    results = []
    
    with tarfile.open(TAR_FILE, "r") as tar:
        for relpath in files:
            try:
                member = tar.getmember(relpath)
                f_txt = tar.extractfile(member)
                lines = f_txt.read().decode("utf-8").splitlines()
                
                # Parse V:N pairs
                V, N = [], []
                for line in lines:
                    if ":" not in line:
                        continue
                    v_str, n_str = line.split(":", 1)
                    try:
                        V.append(int(v_str))
                        N.append(int(n_str))
                    except ValueError:
                        continue

                if len(V) < 2:
                    print(f"Warning: insufficient data in '{relpath}'")
                    continue

                # Fit K and β
                try:
                    K_fit, beta_fit = fit_heaps_scipy(np.array(V), np.array(N))
                    results.append(f"{relpath}  K: {K_fit:.4g}  β: {beta_fit:.4g}")
                except Exception as e:
                    print(f"Warning: fit failed for '{relpath}': {e}")
                    continue
                    
            except KeyError:
                print(f"Warning: '{relpath}' not found in tar file")
                continue
            except Exception as e:
                print(f"Warning: error processing '{relpath}': {e}")
                continue

    out_fname = os.path.join(OUTPUT_DIR, f"results_bucket_{args.bucket}.txt")
    with open(out_fname, "w") as out:
        for result in results:
            out.write(result + "\n")


if __name__ == "__main__":
    main()