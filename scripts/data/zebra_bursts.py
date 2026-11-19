"""Port the released Lorenz LMAN IFR/participation detector to Colab Python."""

import numpy as np
from scipy.signal import convolve, find_peaks, peak_widths


GRID_HZ = 2500
WINDOW_TICKS = 125
STEP_TICKS = 25
SIGMA_SECONDS = 0.015


def participation(spike_times, duration):
    """Evaluate Lorenz's >1 high-IFR sample rule on the 50 ms windows."""
    grid_size = int(np.floor(duration * GRID_HZ)) + 1
    half = WINDOW_TICKS // 2
    centers = np.arange(half, grid_size - half - 1, STEP_TICKS)
    units = np.zeros(len(centers), dtype=np.float64)
    for spikes in spike_times:
        short = np.diff(spikes) <= 0.01
        starts = np.floor(spikes[:-1][short] * GRID_HZ).astype(np.int64) + 1
        stops = np.ceil(spikes[1:][short] * GRID_HZ).astype(np.int64)
        keep = stops > starts
        starts, stops = starts[keep], stops[keep]
        if not len(starts):
            continue
        total = np.r_[0, np.cumsum(stops - starts)]
        edges = np.stack((centers - half, centers + half + 1))
        pos = np.searchsorted(starts, edges, side="left")
        before = total[pos] - np.where(pos, np.maximum(stops[pos - 1] - edges, 0), 0)
        units += before[1] - before[0] > 1
    return centers, units / len(spike_times)


def check_centers(spike_times, centers, observed):
    """Independent, dense evaluation at supplied real or synthetic centers."""
    half = WINDOW_TICKS // 2
    for center, value in zip(centers, observed, strict=True):
        grid = np.arange(center - half, center + half + 1) / GRID_HZ
        active = 0
        for spikes in spike_times:
            idx = np.searchsorted(spikes, grid, side="right")
            inside = (idx > 0) & (idx < len(spikes))
            left, right = spikes[idx[inside] - 1], spikes[idx[inside]]
            positive = (grid[inside] > left) & (grid[inside] < right) & ((right - left) <= 0.01)
            active += np.count_nonzero(positive) > 1
        if not np.isclose(active / len(spike_times), value, atol=1e-12):
            raise ValueError(f"Burst participation check disagrees at {center / GRID_HZ:.6f} s")


def detect(spike_times, duration):
    centers, fraction = participation(spike_times, duration)
    sigma = SIGMA_SECONDS * GRID_HZ / STEP_TICKS
    kernel = np.exp(-0.5 * (np.arange(-80, 81) / sigma) ** 2)
    smooth = convolve(fraction, kernel / kernel.sum(), mode="same", method="direct")
    candidates, _ = find_peaks(smooth)
    checks = candidates[np.argsort(smooth[candidates])[-3:]]
    check_centers(spike_times, centers[checks], fraction[checks])
    step = STEP_TICKS / GRID_HZ
    time = centers / GRID_HZ
    methods = {}
    for key, percentile in (("released", 97.5), ("paper99", 99)):
        threshold = float(np.percentile(smooth[candidates], percentile, method="hazen"))
        peaks = candidates[smooth[candidates] >= threshold]
        widths, _, left, right = peak_widths(smooth, peaks)
        onsets, offsets = time[0] + left * step, time[0] + right * step
        sizes = np.zeros(len(peaks), dtype=np.int64)
        recruited = np.zeros(len(peaks), dtype=np.int64)
        for spikes in spike_times:
            count = np.searchsorted(spikes, offsets) - np.searchsorted(spikes, onsets)
            sizes += count
            recruited += count > 0
        events = [{"peak_s": float(time[peak]), "onset_s": float(onsets[i]),
                   "offset_s": float(offsets[i]), "duration_s": float(widths[i] * step),
                   "event_spikes": int(sizes[i]), "participating_fraction": float(fraction[peak]),
                   "spiking_unit_fraction": float(recruited[i] / len(spike_times))}
                  for i, peak in enumerate(peaks)]
        methods[key] = {"events": events, "threshold": threshold, "peak_height_percentile": percentile}
    return methods, time, fraction, [float(time[i]) for i in checks]


def summarize(events, time, fraction, bounds):
    """Use only burst intervals fully inside the same continuously scored segment."""
    events = sorted(events, key=lambda pair: pair["peak_s"])
    peaks = np.array([e["peak_s"] for e in events])
    keep = []
    gaps = []
    synchrony = []
    for start, stop in bounds:
        left, right = np.searchsorted(peaks, [start, stop])
        here = [e for e in events[left:right] if e["onset_s"] >= start and e["offset_s"] <= stop]
        keep.extend(here)
        gaps.extend(np.diff([e["peak_s"] for e in here]).tolist())
        left, right = np.searchsorted(time, [start + 0.025, stop - 0.025])
        synchrony.extend(fraction[left:right])
    seconds = sum(stop - start for start, stop in bounds)
    return {"events": len(keep), "event_rate_per_min": len(keep) * 60 / seconds,
            "median_duration_ms": float(1000 * np.median([e["duration_s"] for e in keep])) if keep else None,
            "median_size_spikes": float(np.median([e["event_spikes"] for e in keep])) if keep else None,
            "mean_participation": float(np.mean([e["participating_fraction"] for e in keep])) if keep else None,
            "mean_synchrony": float(np.mean(synchrony)) if synchrony else None,
            "interburst_seconds": gaps}


if __name__ == "__main__":
    spikes = [np.array([0.4, 0.404, 0.408, 0.9]), np.array([0.39, 0.395, 0.4, 0.95])]
    centers, values = participation(spikes, 1.0)
    check_centers(spikes, centers, values)
    assert values.max() == 1 and values[0] == 0
