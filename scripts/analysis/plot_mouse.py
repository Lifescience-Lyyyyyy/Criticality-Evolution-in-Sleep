"""Plot every paired mouse session and the state means."""

from pathlib import Path
from statistics import mean
import json

import matplotlib as mpl
import matplotlib.pyplot as plt


config = json.loads(Path("config/paths.json").read_text())
result = Path(config["results_dir"]) / "pnas/mouse.json"
cohort = {}
for row in json.loads(result.read_text())["primary"]:
    cohort.setdefault(row["session"], {})[row["state"]] = row["d2_bits_per_second"]
pairs = [(states["WAKE"], states["NREM"]) for states in cohort.values()]
target = Path(config["figures_dir"]) / "pnas/mouse_wake_nrem"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.7,
    "legend.frameon": False,
})
fig, ax = plt.subplots(figsize=(3.45, 2.65), layout="constrained")
for index, pair in enumerate(pairs):
    ax.plot([0, 1], pair, color="#8e98a4", linewidth=0.85, marker="o",
            markersize=3, alpha=0.76, label="Individual Sessions" if index == 0 else None)
ax.plot([0, 1], [mean(state) for state in zip(*pairs)], color="#126782",
        linewidth=1.8, marker="D", markersize=5, label="Mean Across Sessions")
ax.set_xticks([0, 1], ["Wake", "NREM"])
ax.set(xlim=(-0.25, 1.25), ylim=(0, None), ylabel="d₂ (bits/s)")
ax.margins(y=0.12)
ax.legend(loc="upper left", fontsize=7)
fig.savefig(f"{target}.pdf", facecolor="white")
fig.savefig(f"{target}.svg", facecolor="white")
fig.savefig(f"{target}.png", dpi=600, facecolor="white")
plt.close(fig)
