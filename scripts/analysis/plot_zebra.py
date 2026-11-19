"""Render the five version-2 zebra finch figures from verified numerical outputs."""

from pathlib import Path
import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from summarize_zebra import PAIRS, STATES, valid


paths = json.loads(Path("config/paths.json").read_text())
result = json.loads((Path(paths["results_dir"]) / "zebra_finch/lman.json").read_text())
assert result["schema_version"] == 2
birds = sorted(result["state_metrics"])
palette = {"ZF1": "#126782", "ZF2": "#86543d"}
markers = {"ZF1": "o", "ZF2": "D"}
labels = ["Light-on\nactive", "Light-on\nnon-active", "Scored\nSWS"]
output = Path(paths["figures_dir"]) / "zebra_finch"
output.mkdir(parents=True, exist_ok=True)
broad = {(x["bird"], x["state"]): x for x in result["diagnostics"] if x["variant"] == "broad"}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "svg.fonttype": "none", "pdf.fonttype": 42,
    "axes.spines.right": False, "axes.spines.top": False, "axes.linewidth": 0.7,
    "legend.frameon": False,
})


def save(fig, name):
    if len(fig.axes) > 1 and os.getenv("NATURE_FIGURE_AUDIT"):
        from audit_panel_alignment import require_matplotlib_panel_alignment
        require_matplotlib_panel_alignment(fig, json_out=f"logs/{name}.alignment.json",
                                           tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    base = output / name
    fig.savefig(f"{base}.pdf", facecolor="white")
    fig.savefig(f"{base}.svg", facecolor="white")
    fig.savefig(f"{base}.png", facecolor="white", dpi=600)
    plt.close(fig)


def by_state(ax, data, ylabel):
    for bird in birds:
        values = data[bird]
        ax.plot(range(3), values, color=palette[bird], marker=markers[bird],
                linewidth=1.35, markersize=4, label=bird)
    ax.set_xticks(range(3), labels, fontsize=7)
    ax.set(ylabel=ylabel, xlim=(-0.22, 2.22))
    ax.margins(y=0.17)


fig, ax = plt.subplots(figsize=(3.50, 2.55), layout="constrained")
by_state(ax, {b: [broad[b, s]["d2"] for s in STATES] for b in birds}, "LMAN d₂ (bits/s)")
fig.legend(loc="outside upper center", ncol=2, fontsize=7, title="Broad units (diagnostic)")
save(fig, "lman_state_comparison")

fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45))
fig.subplots_adjust(left=0.10, right=0.97, bottom=0.205, top=0.81, wspace=0.54)
for ax, key, label in zip(axes, ("firing_rate_hz", "fano", "zero_fraction"),
                          ("Population rate (spikes/s)", "Fano factor", "Zero-count fraction"), strict=True):
    by_state(ax, {b: [broad[b, s][key] for s in STATES] for b in birds}, label)
fig.legend(*axes[0].get_legend_handles_labels(), loc="outside upper center", ncol=2, fontsize=7)
save(fig, "lman_firing_statistics")

fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45))
fig.subplots_adjust(left=0.10, right=0.97, bottom=0.205, top=0.81, wspace=0.54)
for ax, key, label in zip(axes, ("event_rate_per_min", "median_size_spikes", "mean_participation"),
                          ("Bursts per minute", "Median spikes per event", "Mean participation"), strict=True):
    by_state(ax, {b: [result["state_metrics"][b][s]["burst_summary"][key] for s in STATES] for b in birds}, label)
fig.legend(*axes[0].get_legend_handles_labels(), loc="outside upper center", ncol=2, fontsize=7)
save(fig, "lman_burst_statistics")

fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.3), sharex="col", sharey=True)
fig.subplots_adjust(left=0.10, right=0.97, bottom=0.12, top=0.97, wspace=0.27, hspace=0.18)
metrics = (("event_rate_per_min", "Bursts per minute"),
           ("median_size_spikes", "Median spikes per event"),
           ("mean_participation", "Mean participation"))
for bird, row in zip(birds, axes, strict=True):
    windows = [r for r in result["windows"] if r["bird"] == bird and r["state"] == "SWS" and valid(r)]
    row[0].set_ylabel(f"{bird} SWS d₂ (bits/s)")
    for ax, (key, _) in zip(row, metrics, strict=True):
        points = [(r[key], r["d2"]) for r in windows if r.get(key) is not None]
        if points:
            x, y = zip(*points)
            ax.scatter(x, y, s=18, color=palette[bird], marker=markers[bird], alpha=0.85)
for ax, (_, label) in zip(axes[-1], metrics, strict=True):
    ax.set_xlabel(label)
save(fig, "lman_burst_d2")

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.05), sharex=True, sharey=True, layout="constrained")
variants = [("broad", "Broad units", "#126782", "o"),
            ("reach_0p40", "0.40 s AR reach", "#86543d", "D"),
            ("unit_subsampling", "100 unit subsamples", "#426a40", "s"),
            ("rate_thinning", "100 rate thinnings", "#80518a", "^"),
            ("historical", "Historical broad", "#66666a", "x")]
prior = {(r["bird"],r["state"]): r["d2_bits_per_second"] for r in result["historical_reproduction"]["historical_primary"]}
for b, ax in zip(birds, axes, strict=True):
    for i, (hi, lo) in enumerate(PAIRS):
        row = result["contrasts"][b][f"{hi}_minus_{lo}"]
        for j, (kind, name, color, mark) in enumerate(variants):
            if kind == "historical":
                if (hi, lo) != ("SWS", "ACTIVE"):
                    continue
                x = prior[b, "SWS"] - prior[b, "WAKE"]
            else:
                x = row[kind]
            y = i + (j - 2) * 0.125
            legend = name if i == (1 if kind == "historical" else 0) else None
            if isinstance(x, dict):
                if x["n"]:
                    ax.errorbar(x["median"], y, xerr=[[x["median"] - x["low"]], [x["high"] - x["median"]]],
                                marker=mark, color=color, markersize=4, linewidth=1, capsize=2, label=legend)
            elif x is not None:
                ax.plot(x, y, marker=mark, color=color, markersize=4, linestyle="none", label=legend)
    ax.axvline(0, color="#747981", linewidth=0.7)
    ax.set(title=b, xlabel="Paired Δd₂ (bits/s)", ylim=(2.43, -0.43))
    ax.set_yticks(range(3), ["Non-active − active", "SWS − active", "SWS − non-active"], fontsize=7)
strict_note = "Strict single-unit pairs unavailable"
strict = [r for r in result["diagnostics"] if r["variant"] == "strict_unpaired" and valid(r)]
if len(strict) == 1:
    row = strict[0]
    strict_note += f" · {row['bird']} {row['state']} d₂ = {row['d2']:.3f} bits/s"
fig.legend(*axes[0].get_legend_handles_labels(), loc="outside upper center", ncol=3, fontsize=7,
           title=strict_note, title_fontsize=7)
save(fig, "lman_sensitivity")
