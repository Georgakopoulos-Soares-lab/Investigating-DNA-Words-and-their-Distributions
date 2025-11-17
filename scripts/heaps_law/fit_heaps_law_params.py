#!/usr/bin/env python3
import json
import os
import argparse
import numpy as np
from scipy.optimize import curve_fit

def heaps_model(N, K, beta):
    """Heaps' law: V = K * N^beta."""
    return K * (N ** beta)

def menzerath_model(N, A, b):
    """Menzerath-like: M = A * N^b, where M = V/N."""
    return A * (N ** b)

def fit_heaps_scipy(V, N):
    """
    Fit Heaps' law V = K * N^β using SciPy curve_fit on linear scale.
    """
    V = np.asarray(V, dtype=float)
    N = np.asarray(N, dtype=float)

    # Keep only strictly positive pairs
    mask = (V > 0) & (N > 0)
    V = V[mask]
    N = N[mask]

    if V.size < 2:
        raise ValueError("Need at least two positive (V, N) points for Heaps fit")

    # Initial guess: K ~ V[0] / N[0]^0.5, beta ~ 0.5
    beta0 = 0.5
    K0 = V[0] / (N[0] ** beta0) if N[0] > 0 else 1.0

    popt, _ = curve_fit(
        heaps_model,
        N,
        V,
        p0=[K0, beta0],
        maxfev=10_000
    )
    K_fit, beta_fit = popt
    return K_fit, beta_fit


def fit_menzerath_from_pairs(V, N):
    """
    From progressive (V, N) pairs, fit the Menzerath-like relation:

        M(N) = V(N) / N = A * N^b

    on linear scale via curve_fit with bounds:
        A in [0, 1e9],  b in [-5, 5]
    """
    V = np.asarray(V, dtype=float)
    N = np.asarray(N, dtype=float)

    # Keep only strictly positive pairs
    mask = (V > 0) & (N > 0)
    V = V[mask]
    N = N[mask]

    if V.size < 2:
        raise ValueError("Need at least two positive (V, N) points for Menzerath fit")

    M = V / N  # vocabulary density

    # Bounds
    A_low, A_high = 0.0, 1e9
    b_low, b_high = -5.0, 5.0

    # Initial guess
    b0 = -0.5
    # Pivot at geometric mean N
    N0 = float(np.exp(np.mean(np.log(N))))
    M0 = float(np.interp(N0, N, M))
    A0 = M0 / (N0 ** b0)

    # Clip initial guess into bounds
    A0 = min(max(A0, A_low), A_high)
    b0 = min(max(b0, b_low), b_high)

    popt, _ = curve_fit(
        menzerath_model,
        N,
        M,
        p0=[A0, b0],
        bounds=([A_low, b_low], [A_high, b_high]),
        maxfev=10_000,
    )
    A_fit, b_fit = popt
    return A_fit, b_fit

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Fit Heaps' law (V = K N^β) and Menzerath (M = V/N = A N^b) "
            "for each genome file in a scheduler bucket."
        )
    )
    parser.add_argument("scheduler", help="Path to scheduler JSON")
    parser.add_argument("bucket", help="Bucket ID (key in scheduler JSON)")
    parser.add_argument("input_dir", help="Directory containing V:N text files")
    parser.add_argument("output_dir", help="Directory to write the results file")
    parser.add_argument(
        "--menzerath-only",
        action="store_true",
        help="Only fit Menzerath (skip Heaps); useful if Heaps already computed."
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load scheduler
    with open(args.scheduler) as f:
        sched = json.load(f)

    if args.bucket not in sched:
        raise SystemExit(f"Error: bucket '{args.bucket}' not found in scheduler.")

    files = sched[args.bucket]
    print(f"Processing bucket {args.bucket} with {len(files)} files")
    print(f"Input directory:  {args.input_dir}")
    print(f"Output directory: {args.output_dir}")

    results_lines = []

    for relpath in files:
        full_path = os.path.join(args.input_dir, relpath)

        try:
            if not os.path.exists(full_path):
                print(f"Warning: '{full_path}' not found on disk, skipping")
                continue

            with open(full_path, "r") as f_txt:
                lines = f_txt.read().splitlines()

            # Parse V:N pairs
            V_vals, N_vals = [], []
            for line in lines:
                if ":" not in line:
                    continue
                v_str, n_str = line.split(":", 1)
                try:
                    V_vals.append(int(v_str))
                    N_vals.append(int(n_str))
                except ValueError:
                    continue

            if len(V_vals) < 2:
                print(f"Warning: insufficient data in '{relpath}', skipping")
                continue

            V_arr = np.array(V_vals, dtype=float)
            N_arr = np.array(N_vals, dtype=float)

            # ---- Heaps fit ----
            if not args.menzerath_only:
                try:
                    K_fit, beta_fit = fit_heaps_scipy(V_arr, N_arr)
                except Exception as e:
                    print(f"Warning: Heaps fit failed for '{relpath}': {e}")
                    K_fit, beta_fit = np.nan, np.nan
            else:
                K_fit, beta_fit = np.nan, np.nan

            # ---- Menzerath fit (M = V/N vs N, bounded) ----
            try:
                A_fit, b_fit = fit_menzerath_from_pairs(V_arr, N_arr)
            except Exception as e:
                print(f"Warning: Menzerath fit failed for '{relpath}': {e}")
                A_fit, b_fit = np.nan, np.nan

            # ---- Compose output line ----
            if args.menzerath_only:
                line_out = (
                    f"{relpath}\t"
                    f"A_M:{A_fit:.4g}\t"
                    f"b_M:{b_fit:.4g}"
                )
            else:
                if np.isnan(beta_fit):
                    b_theory = np.nan
                    delta_b = np.nan
                else:
                    b_theory = beta_fit - 1.0
                    delta_b = b_fit - b_theory if not np.isnan(b_fit) else np.nan

                line_out = (
                    f"{relpath}\t"
                    f"K:{K_fit:.4g}\t"
                    f"β:{beta_fit:.4g}\t"
                    f"A_M:{A_fit:.4g}\t"
                    f"b_M:{b_fit:.4g}\t"
                    f"(b_M-(β-1)):{delta_b:.4g}"
                )

            results_lines.append(line_out)

        except Exception as e:
            print(f"Warning: error processing '{relpath}': {e}")
            continue

    # Write results
    out_fname = os.path.join(args.output_dir, f"results_bucket_{args.bucket}.txt")
    with open(out_fname, "w") as out_f:
        for line in results_lines:
            out_f.write(line + "\n")

    print(f"Wrote results to: {out_fname}")


if __name__ == "__main__":
    main()
