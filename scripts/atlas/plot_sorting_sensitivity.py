"""Plot Phase 11 sorting sensitivities at fixed state intervals."""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from build_atlas_data import ROOT, read
from build_audit_results import OUT, thresholds

sys.path.insert(0, str(ROOT / "scripts/analysis"))
from plot_sleep_dandi import save


def plot():
    cases = read(OUT + "audit_comparisons.json")["cases"]
    picked = []
    for c in cases:
        for region, methods in thresholds(c).items():
            for stem in ("lratio", "contamination"):
                if stem + "_q50" in methods:
                    order = [methods[k] for k in (stem + "_q50", stem + "_q75", "reference")]
                    picked.append((c["session"], region, stem, order))
    assert len(picked) == 5
    fig, axs = plt.subplots(2, 3, figsize=(8.1, 4.45), sharey=True)
    fig.subplots_adjust(left=.095, right=.985, bottom=.17, top=.92, wspace=.32, hspace=.95)
    colors = {"NW": "#a85e2c", "RN": "#2c667d"}
    labels = {"NW": "NREM − Wake", "RN": "REM − NREM"}
    for i, (session, region, stem, methods) in enumerate(picked):
        ax = axs.flat[i]
        ax.axhline(0, color=".55", ls="--", lw=.65)
        units = [next(iter(x["states"].values()))["units"] for x in methods]
        names = [f"q50\n{units[0]} units", f"q75\n{units[1]} units", f"Original\n{units[2]} units"]
        for contrast in ("NW", "RN"):
            vals = [m["effects"][contrast] for m in methods]
            ax.plot(range(3), vals, "-o", label=labels[contrast], color=colors[contrast], linewidth=1, markersize=3.8)
        ax.set(xticks=range(3), xticklabels=names, xlim=(-.17, 2.17),
               title=region + " · " + session)
        ax.tick_params(axis="both", labelsize=7)
        if i % 3 == 0:
            ax.set_ylabel("Δd₂ (bits/s)")
        ax.set_xlabel("L-ratio" if stem == "lratio" else "ContamPct", color=".36", fontsize=7, labelpad=2)
    axs.flat[5].axis("off")
    handles, legend = axs.flat[0].get_legend_handles_labels()
    axs.flat[5].legend(handles, legend, loc="center", frameon=False, fontsize=8)
    destination = ROOT / "figures/phase11/sorting_sensitivity"
    destination.parent.mkdir(parents=True, exist_ok=True)
    save(fig, destination)


if __name__ == "__main__":
    plot()
