"""Plot the seeded AR(1) validation against its analytical reference."""

from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt


config = json.loads(Path("config/paths.json").read_text())
source = Path(config["results_dir"]) / "pnas/synthetic.json"
cases = json.loads(source.read_text())["ar1"]
target = Path(config["figures_dir"]) / "pnas/synthetic_d2"
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
fig, ax = plt.subplots(figsize=(3.45, 2.42), layout="constrained")
xs = list(range(len(cases)))
references = [row["true_analytic_d2"] for row in cases]
estimated = [row["d2_bits_per_sample"] for row in cases]
ax.scatter([x - 0.07 for x in xs], references, marker="D", s=27,
           facecolor="white", edgecolor="#535b67", label="Analytic reference")
ax.scatter([x + 0.07 for x in xs], estimated, s=29,
           color="#126782", label="MATLAB estimate")
ax.set_xticks(xs, [f"{row['generating_phi']:.2f}" for row in cases])
ax.set(xlim=(-0.35, len(cases) - 0.65), ylim=(0, 0.255),
       xlabel="Generating AR(1) coefficient", ylabel="d₂ (bits/sample)")
ax.legend(loc="upper right", fontsize=7)
fig.savefig(f"{target}.pdf", facecolor="white")
fig.savefig(f"{target}.svg", facecolor="white")
fig.savefig(f"{target}.png", dpi=600, facecolor="white")
plt.close(fig)
