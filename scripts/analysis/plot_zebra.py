"""Show the paired LMAN d2 estimates at the bird level."""

from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt


paths = json.loads(Path("config/paths.json").read_text())
source = Path(paths["results_dir"]) / "zebra_finch/lman.json"
pairs = {}
for row in json.loads(source.read_text())["primary"]:
    pairs.setdefault(row["bird"], {})[row["state"]] = row["d2_bits_per_second"]
target = Path(paths["figures_dir"]) / "zebra_finch/lman_wake_sws"
target.parent.mkdir(parents=True, exist_ok=True)

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
fig, ax = plt.subplots(figsize=(3.45, 2.55), layout="constrained")
for (bird, states), color, marker in zip(sorted(pairs.items()), ["#126782", "#86543d"], ["o", "D"], strict=True):
    ax.plot([0, 1], [states["WAKE"], states["SWS"]], color=color, marker=marker,
            linewidth=1.4, markersize=5, label=bird)
ax.set_xticks([0, 1], ["Light-on active", "Scored SWS"])
upper = max(value for states in pairs.values() for value in states.values())
ax.set(xlim=(-0.25, 1.25), ylim=(0, 1.1 * upper), ylabel="LMAN d₂ (bits/s)")
ax.legend(loc="lower center", ncol=2, fontsize=7)
fig.savefig(f"{target}.pdf", facecolor="white")
fig.savefig(f"{target}.svg", facecolor="white")
fig.savefig(f"{target}.png", dpi=600, facecolor="white")
plt.close(fig)
