#!/usr/bin/env python3

import sys
import os
import json
import numpy as np
from scipy.optimize import curve_fit

def truncated_power_law(k, alpha, lambda_, scale):
    """
    Truncated Power-Law model:
        freq(k) = scale * (k ** -alpha) * exp(-lambda * k)
    """
    return scale * (k ** -alpha) * np.exp(-lambda_ * k)

def fit_statistics(y_true, y_pred, num_params):
    
    sse = np.sum((y_true - y_pred) ** 2)
    n = len(y_true)
    ss_total = np.sum((y_true - np.mean(y_true)) ** 2)
    r2 = 1 - (sse / ss_total) if ss_total > 0 else np.nan

    if sse <= 0:
        aic = 2 * num_params
    else:
        aic = n * np.log(sse / n) + 2 * num_params

    return {
        'R2': r2,
        'AIC': aic
    }

def main():
    if len(sys.argv) < 3:
        print("Usage: python fit_truncated_powerlaw.py <bucket_id> <scheduler_file.json>")
        sys.exit(1)

    bucket_id_str = sys.argv[1]
    scheduler_file = sys.argv[2]

    with open(scheduler_file, 'r') as sf:
        data = json.load(sf)

    if bucket_id_str not in data:
        print(f"No entry found in JSON for bucket_id='{bucket_id_str}'")
        sys.exit(0)

    file_list = data[bucket_id_str]
    if not file_list:
        print(f"No files listed for bucket {bucket_id_str}")
        sys.exit(0)

    out_filename = f"truncated_power_law_3mers_{bucket_id_str}.txt"
    with open(out_filename, 'w') as out_f:
        for filepath in file_list:
            try:
                counts = []
                with open(filepath, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) != 2:
                            continue
                        try:
                            counts.append(float(parts[1]))
                        except ValueError:
                            continue
            except FileNotFoundError:
                out_f.write(f"{os.path.basename(filepath)}: FileNotFound\n")
                continue

            if len(counts) == 0:
                continue

            counts = np.array(counts, dtype=float)
            k_array = np.arange(1, len(counts) + 1, dtype=float)
            freq = counts / np.sum(counts)

            initial_guess = [1, 0.1, 1]
            try:
                popt, _ = curve_fit(
                    truncated_power_law,
                    k_array,
                    freq,
                    p0=initial_guess,
                    bounds=(0, np.inf),
                    method='trf',
                    max_nfev=1000000
                )
            except RuntimeError as e:
                out_f.write(f"{os.path.basename(filepath)}: FitError={str(e)}\n")
                continue

            alpha_fit, lambda_fit, scale_fit = popt
            pred = truncated_power_law(k_array, alpha_fit, lambda_fit, scale_fit)

            stats = fit_statistics(freq, pred, len(popt))

            out_f.write(
                f"{os.path.basename(filepath)}:"
                f" alpha={alpha_fit:.6g}"
                f" lambda={lambda_fit:.6g}"
                f" scale={scale_fit:.6g}"
                f" R2={stats['R2']:.4g}"
                f" AIC={stats['AIC']:.4g}\n"
            )

    print(f"Done. Results saved to: {out_filename}")

if __name__ == "__main__":
    main()