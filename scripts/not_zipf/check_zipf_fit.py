#!/usr/bin/env python3

import sys
import os
import json
import numpy as np
from scipy.stats import spearmanr

DATA_DIR_PATH = "/storage/group/izg5139/default/xaris/sorted_4mers" 

def fit_statistics(y_true, y_pred, num_params):
    sse = np.sum((y_true - y_pred) ** 2)
    n = len(y_true)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r2 = 1 - sse / ss_tot if ss_tot > 0 else np.nan
    rmse = np.sqrt(sse / n)
    if sse <= 0:
        aic = 2 * num_params
        bic = num_params * np.log(n)
    else:
        aic = n * np.log(sse / n) + 2 * num_params
        bic = n * np.log(sse / n) + num_params * np.log(n)
    return {'R2': r2, 'RMSE': rmse, 'AIC': aic, 'BIC': bic}

def main():
    if len(sys.argv) < 3:
        print("Usage: python kmers_evaluation_zipf_raw.py <bucket_id> <scheduler_file.json>")
        sys.exit(1)

    bucket_id      = sys.argv[1]
    scheduler_file = sys.argv[2]

    with open(scheduler_file) as sf:
        sched = json.load(sf)

    if bucket_id not in sched or not sched[bucket_id]:
        print(f"No files listed for bucket '{bucket_id}'")
        sys.exit(0)

    out_fname = f"not_zipf_law_4mers_{bucket_id}.txt"

    with open(out_fname, 'w') as out_f:
        for file_path in sched[bucket_id]:
            full_path = os.path.join(DATA_DIR_PATH, file_path)
            if not os.path.isfile(full_path):
                out_f.write(f"{os.path.basename(file_path)}: FileNotFound\n")
                continue

            try:
                with open(full_path, 'r') as f:
                    lines = [line.strip() for line in f]
            except Exception:
                out_f.write(f"{os.path.basename(file_path)}: ReadError\n")
                continue

            counts = []
            for line in lines:
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    continue
                try:
                    counts.append(float(parts[1]))
                except ValueError:
                    pass

            if not counts:
                out_f.write(f"{os.path.basename(file_path)}: NoData\n")
                continue

            counts = np.array(counts, dtype=float)
            k = np.arange(1, len(counts) + 1, dtype=float)
            pred = counts[0] / k

            stats = fit_statistics(counts, pred, num_params=1)
            rho, _ = spearmanr(counts, pred)

            out_f.write(
                f"{os.path.basename(file_path)}:"
                f" R2={stats['R2']:.4g}"
                f" Spearman={rho:.4g}"
                f" RMSE={stats['RMSE']:.4g}"
                f" AIC={stats['AIC']:.4g}"
                f" BIC={stats['BIC']:.4g}\n"
            )

    print(f"Done. Results saved to: {out_fname}")

if __name__ == "__main__":
    main()