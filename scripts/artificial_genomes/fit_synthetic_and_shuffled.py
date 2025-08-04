#!/usr/bin/env python3
import re
from pathlib import Path
import numpy as np
from scipy.optimize import curve_fit

def truncated_power_law(r, alpha, lam, scale):
    return scale * r ** (-alpha) * np.exp(-lam * r)

def zipf_mandelbrot(r, alpha, beta, scale):
    return scale * (r + beta) ** (-alpha)

def fit_and_evaluate(model_func, x, y, p0, bounds, maxfev=1000000):
    try:
        popt, _ = curve_fit(
            model_func, x, y,
            p0=p0, bounds=bounds,
            maxfev=maxfev
        )
    except RuntimeError as e:
        # fit failed
        return None, None, None, None, None

    y_pred = model_func(x, *popt)
    resid = y - y_pred
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res/ss_tot

    n = len(y)
    k = len(popt)
    aic = n * np.log(ss_res/n) + 2 * k
    bic = n * np.log(ss_res/n) + k * np.log(n)
    rmse = np.sqrt(ss_res / n)

    return popt, r2, aic, bic, rmse

def load_counts(path: Path):
    """
    Load a two-column tab file (kmer<TAB>count), return a sorted list of counts.
    """
    lines = path.read_text().strip().splitlines()
    counts = [int(line.split('\t')[1]) for line in lines if line]
    return sorted(counts, reverse=True)

def prepare_rank_freq(counts):
    total = sum(counts)
    freqs = np.array(counts, dtype=float) / total
    x = np.arange(1, len(freqs) + 1)
    return x, freqs

def main(input_dir, output_summary):
    p = re.compile(r'(.+)_([34567])\.txt$')
    files = list(Path(input_dir).glob("*.txt"))

    genomes = {}
    for f in files:
        m = p.match(f.name)
        if not m:
            continue
        prefix, k = m.group(1), int(m.group(2))
        genomes.setdefault(prefix, []).append((k, f))

    with open(output_summary, "w") as out:
        for prefix in sorted(genomes):
            out.write(f"{prefix}\n")

            entries = sorted(genomes[prefix], key=lambda t: t[0])

            # truncated power law
            out.write("Truncated power law:\n")
            for k, path in entries:
                counts = load_counts(path)
                x, y = prepare_rank_freq(counts)

                popt, r2, aic, bic, rmse = fit_and_evaluate(
                    truncated_power_law,
                    x, y,
                    p0=(1.5, 0.1, 1.0),
                    bounds=(0, np.inf),
                    maxfev=1000000 
                )
                if popt is None:
                    out.write(f" {k}: fit failed\n")
                else:
                    alpha, lam, scale = popt
                    out.write(
                        f" {k}: alpha={alpha:.4f} "
                        f"lambda={lam:.6f} "
                        f"scale={scale:.4e} "
                        f"R2={r2:.4f} "
                        f"AIC={aic:.2f} "
                        f"BIC={bic:.2f} "
                        f"RMSE={rmse:.6f}\n"
                    )

            # Zipf–Mandelbrot
            out.write("Zipf–Mandelbrot:\n")
            for k, path in entries:
                counts = load_counts(path)
                x, y = prepare_rank_freq(counts)

                popt, r2, aic, bic, rmse = fit_and_evaluate(
                    zipf_mandelbrot,
                    x, y,
                    p0=(1.5, 1.0, 1.0),
                    bounds=(0, np.inf),
                    maxfev=1000000
                )
                if popt is None:
                    out.write(f" {k}: fit failed\n")
                else:
                    alpha, beta, scale = popt
                    out.write(
                        f" {k}: alpha={alpha:.4f} "
                        f"beta={beta:.4f} "
                        f"scale={scale:.4e} "
                        f"R2={r2:.4f} "
                        f"AIC={aic:.2f} "
                        f"BIC={bic:.2f} "
                        f"RMSE={rmse:.6f}\n"
                    )

            out.write("\n")

if __name__ == "__main__":
    INPUT_DIR      = "/scratch/cpk5664/shuffled_kmers"
    OUTPUT_SUMMARY = "/scratch/cpk5664/results_fit_shuffled.txt"
    main(INPUT_DIR, OUTPUT_SUMMARY)