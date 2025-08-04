import re
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

path_shuffled = "/scratch/cpk5664/results_fit_shuffled.txt"
path_synthetic = "/scratch/cpk5664/results_fit_synthetic.txt"
path_real = "/scratch/cpk5664/results_real.txt"

model_headers = {
    "truncated": re.compile(r"^Truncated power law:", re.IGNORECASE),
    "zipf": re.compile(r"^Zipf[\u2013-]Mandelbrot:", re.IGNORECASE),
}

pattern_truncated = re.compile(
    r"^\s*(\d+):\s*alpha=([\d.eE+-]+)\s+lambda=([\d.eE+-]+)\s+scale=([\d.eE+-]+)\s+R2=([\d.eE+-]+)"
)
pattern_zipf = re.compile(
    r"^\s*(\d+):\s*alpha=([\d.eE+-]+)\s+beta=([\d.eE+-]+)\s+scale=([\d.eE+-]+)\s+R2=([\d.eE+-]+)"
)

pattern_truncated_real = re.compile(
    r"^\s*(\d+):\s*alpha=([\d.eE+-]+)\s+lambda=([\d.eE+-]+)\s+scale=([\d.eE+-]+).*?R2=([\d.eE+-]+)"
)
pattern_zipf_real = re.compile(
    r"^\s*(\d+):\s*alpha=([\d.eE+-]+)\s+beta=([\d.eE+-]+)\s+scale=([\d.eE+-]+).*?R2=([\d.eE+-]+)"
)

def parse_file(path, label):
    records = []
    if not os.path.exists(path):
        print(f"[WARN] File not found: {path}")
        return records
    current_genome = None
    current_model = None
    with open(path, 'r', encoding='utf-8') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if all(not hdr.match(line) for hdr in model_headers.values()) and ":" not in line:
                current_genome = line
                current_model = None
                continue
            model_switched = False
            for key, hdr in model_headers.items():
                if hdr.match(line):
                    current_model = key
                    model_switched = True
                    break
            if model_switched:
                continue

            if current_model == "truncated":
                m = pattern_truncated.match(line) or pattern_truncated_real.match(line)
                if m:
                    k, alpha, lam, scale, r2 = m.groups()
                    records.append({
                        "genome": current_genome,
                        "model": "truncated",
                        "k": int(k),
                        "alpha": float(alpha),
                        "lambda": float(lam),
                        "beta": None,
                        "scale": float(scale),
                        "R2": float(r2),
                        "dataset": label
                    })
            elif current_model == "zipf":
                m = pattern_zipf.match(line) or pattern_zipf_real.match(line)
                if m:
                    k, alpha, beta, scale, r2 = m.groups()
                    records.append({
                        "genome": current_genome,
                        "model": "zipf",
                        "k": int(k),
                        "alpha": float(alpha),
                        "lambda": None,
                        "beta": float(beta),
                        "scale": float(scale),
                        "R2": float(r2),
                        "dataset": label
                    })
    return records

data = []
data.extend(parse_file(path_shuffled, "shuffled"))
data.extend(parse_file(path_synthetic, "synthetic"))
data.extend(parse_file(path_real, "real"))

df = pd.DataFrame(data)
if df.empty:
    raise SystemExit("No data parsed. Check file paths and formats.")

print(f"Parsed {len(df)} parameter rows across {df['genome'].nunique()} genomes.")

out_dir = "comparative_plots"
os.makedirs(out_dir, exist_ok=True)

def save_plots(fig, fname):
    path = os.path.join(out_dir, fname)
    fig.savefig(path, dpi=220, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {path}")

colors = {
    "synthetic": "#1f77b4",
    "real": "#ff7f0e",
    "shuffled": "#2ca02c",
}

markers = {
    "synthetic": "o",
    "real": "s",
    "shuffled": "^",
}

# Box plots for all key parameters
def boxplot_param(df, param: str, model: str, fname: str):
    sub = df[(df['model'] == model) & df[param].notnull()]
    if sub.empty:
        print(f"No data for model={model}, param={param}")
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.boxplot(
        data=sub,
        x="k",
        y=param,
        hue="dataset",
        palette=colors,
        ax=ax
    )
    ax.set_xlabel("kmer length", fontsize=16)
    ax.set_ylabel(param, fontsize=16)
    ax.tick_params(axis='both', labelsize=16)
    ax.legend(title="Dataset", fontsize=12, title_fontsize=12)
    fig.tight_layout()
    save_plots(fig, fname)

boxplot_param(df, "alpha", "truncated", "boxplot_alpha_truncated.png")
boxplot_param(df, "lambda", "truncated", "boxplot_lambda_truncated.png")
boxplot_param(df, "alpha", "zipf", "boxplot_alpha_zipf.png")
boxplot_param(df, "beta", "zipf", "boxplot_beta_zipf.png")

# log‐scale β plot
sub = df[(df['model'] == 'zipf') & df['beta'].notnull()]

fig, ax = plt.subplots(figsize=(8, 6))
sns.boxplot(
    data=sub,
    x="k",
    y="beta",
    hue="dataset",
    palette=colors,
    ax=ax,
    showfliers=True
)
ax.set_yscale('log')
ax.set_xlabel("kmer length", fontsize=16)
ax.set_ylabel("beta", fontsize=16)
ax.tick_params(axis='both', labelsize=16)
ax.legend(title="Dataset", fontsize=12, title_fontsize=12)
fig.tight_layout()
save_plots(fig, "boxplot_beta_zipf_log.png")

# log‐scale lambda plot
sub = df[(df['model'] == 'truncated') & df['lambda'].notnull()]

fig, ax = plt.subplots(figsize=(8, 6))
sns.boxplot(
    data=sub,
    x="k",
    y="lambda",
    hue="dataset",
    palette=colors,
    ax=ax,
    showfliers=True
)
ax.set_yscale('log')
ax.set_xlabel("kmer length", fontsize=16)
ax.set_ylabel("lambda", fontsize=16)
ax.tick_params(axis='both', labelsize=16)
ax.legend(title="Dataset", fontsize=12, title_fontsize=12)
fig.tight_layout()
save_plots(fig, "boxplot_lambda_truncated_log.png")