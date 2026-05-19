#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
ALGOS = ["chainx", "chainx-opt", "llchain"]
ALGO_LABELS = {"chainx": "ChainX", "chainx-opt": "ChainX-opt", "llchain": "llchain"}


def load_file(prefix, algo, strip_s=False):
    path = os.path.join(DATA_DIR, f"{prefix}_{algo}")
    anchors, times = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            anchors.append(int(parts[0]))
            t = parts[1].rstrip("s") if strip_s else parts[1]
            times.append(float(t))
    return np.array(anchors), np.array(times)


# Load seeding (has trailing 's') and chaining data for all algorithms
seeding = {a: load_file("seeding_times_human_mem", a, strip_s=True) for a in ALGOS}
chaining = {a: load_file("chaining_times_human_mem", a, strip_s=False) for a in ALGOS}


def avg_seeding(data):
    """Average seeding times across algorithms; use all 3 where available, else 2."""
    min_len = min(len(data[a][0]) for a in ALGOS)
    anchors_full = data["chainx"][0]
    avg_all = np.stack([data[a][1][:min_len] for a in ALGOS], axis=1).mean(axis=1)
    if len(anchors_full) > min_len:
        avg_partial = np.stack(
            [data[a][1][min_len:] for a in ["chainx", "chainx-opt"]], axis=1
        ).mean(axis=1)
        return anchors_full, np.concatenate([avg_all, avg_partial])
    return anchors_full[:min_len], avg_all


def filter_positive(anchors, times):
    mask = (anchors > 0) & (times > 0)
    return anchors[mask], times[mask]


def bin_trend(x, y, n_bins=50):
    """Average into n_bins equal-width bins in log space."""
    edges = np.logspace(np.log10(x.min()), np.log10(x.max()), n_bins + 1)
    idx = np.digitize(x, edges) - 1
    idx = np.clip(idx, 0, n_bins - 1)
    bx, by = [], []
    for i in range(n_bins):
        mask = idx == i
        if mask.sum() > 0:
            bx.append(x[mask].mean())
            by.append(y[mask].mean())
    return np.array(bx), np.array(by)


panels = [
    (
        f"{ALGO_LABELS['chainx']}",
        *filter_positive(*chaining["chainx"]),
        "tab:blue",
        "tab:cyan",
        1.0,
    ),
    (
        f"{ALGO_LABELS['chainx-opt']}",
        *filter_positive(*chaining["chainx-opt"]),
        "tab:orange",
        "tab:red",
        1.0,
    ),
    (
        f"{ALGO_LABELS['llchain']}",
        *filter_positive(*chaining["llchain"]),
        "tab:green",
        "lime",
        1.0,
    ),
    ("Seeding", *filter_positive(*avg_seeding(seeding)), "tab:gray", "black", 0.1),
]

fig, ax = plt.subplots(figsize=(5, 3.5))

for title, x, t, color, c2, s in panels:
    ax.scatter(x, t, s=s, alpha=1.0, color=color, rasterized=True, label=title)
    bx, by = bin_trend(x, t)
    ax.plot(bx, by, color=c2, linewidth=1.5)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("# anchors")
ax.set_ylabel("time [s]")
ax.legend(markerscale=6, loc="upper left")
ax.grid(True, which="major", linestyle="-", linewidth=0.5, alpha=0.7)
# ax.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.5)

plt.tight_layout()
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot.pdf")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Saved to {out_path}")

print(f"\n{'Stage':<20} {'Avg time (s)':>15}")
print("-" * 36)
for title, x, t, color, c2, s in panels:
    print(f"{title:<20} {t.mean():>15.6f}")
