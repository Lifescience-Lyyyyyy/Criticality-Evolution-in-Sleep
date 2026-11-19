"""Summarize sensitivity and within-recording window associations from MATLAB."""

from collections import defaultdict
from pathlib import Path
import json
import sys

import numpy as np
from scipy.stats import rankdata


STATES = ("ACTIVE", "NON_ACTIVE", "SWS")
PAIRS = (("NON_ACTIVE", "ACTIVE"), ("SWS", "ACTIVE"), ("SWS", "NON_ACTIVE"))
METRICS = ("event_rate_per_min", "median_size_spikes", "mean_participation",
           "zero_fraction", "actual_rate_hz", "fano", "mean_synchrony", "paper99_event_rate_per_min")


def valid(row):
    return row["exit_status"] == 0 and row["solver_flag"] is not None and row["solver_flag"] > 0 and row["d2"] is not None


def distribution(values):
    values = np.asarray(values)
    if not len(values):
        return {"n": 0, "low": None, "median": None, "high": None}
    lo, mid, hi = np.quantile(values, [0.025, 0.5, 0.975])
    return {"n": len(values), "low": float(lo), "median": float(mid), "high": float(hi)}


def spearman(left, right):
    if len(left) < 4 or min(np.ptp(left), np.ptp(right)) == 0:
        return None
    return float(np.corrcoef(rankdata(left), rankdata(right))[0, 1])


def block_association(rows, name, seed, block_seconds=300):
    rows = [r for r in rows if valid(r) and r.get(name) is not None]
    if not rows:
        return {"status": "BLOCKED", "n": 0, "blocks": 0, "rho": None}
    xs = np.array([r[name] for r in rows])
    ys = np.array([r["d2"] for r in rows])
    blocks = np.array([int(r["start_s"]) // block_seconds for r in rows])
    groups = [np.flatnonzero(blocks == block) for block in np.unique(blocks)]
    rho = spearman(xs, ys)
    if rho is None or len(groups) < 5:
        return {"status": "BLOCKED", "n": len(rows), "blocks": len(groups), "rho": rho}
    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(2000):
        idx = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
        estimate = spearman(xs[idx], ys[idx])
        if estimate is not None:
            boot.append(estimate)
    return {"status": "DESCRIPTIVE", "n": len(rows), "blocks": len(groups), "rho": rho,
            "bootstrap": distribution(boot), "block_seconds": block_seconds, "seed": seed}


def match_burst_windows(rows, rules):
    rows = [r for r in rows if valid(r)]
    if len(rows) < 6:
        return {"status": "BLOCKED", "reason": "insufficient_windows"}
    low, high = np.quantile([r["event_rate_per_min"] for r in rows], [1/3, 2/3])
    if low >= high:
        return {"status": "BLOCKED", "reason": "tied_tertiles", "thresholds": [float(low), float(high)]}
    poor = [r for r in rows if r["event_rate_per_min"] <= low]
    rich = [r for r in rows if r["event_rate_per_min"] >= high]
    choices = []
    for i, more in enumerate(rich):
        for j, less in enumerate(poor):
            gap = abs(more["start_s"] - less["start_s"])
            a, b = sorted([more["actual_rate_hz"], less["actual_rate_hz"]])
            if a > 0 and b <= a * (1 + rules["rate_matching_caliper_fraction"]) and gap <= rules["time_matching_caliper_seconds"]:
                choices.append((b/a, gap, i, j))
    used_rich, used_poor, pairs = set(), set(), []
    for _, gap, i, j in sorted(choices):
        if i not in used_rich and j not in used_poor:
            more, less = rich[i], poor[j]
            pairs.append({"rich_start_s": more["start_s"], "poor_start_s": less["start_s"],
                          "rich_rate_hz": more["actual_rate_hz"], "poor_rate_hz": less["actual_rate_hz"],
                          "window_separation_s": gap, "rich_minus_poor_d2": more["d2"] - less["d2"]})
            used_rich.add(i)
            used_poor.add(j)
    return {"status": "DESCRIPTIVE" if len(pairs) >= 3 else "BLOCKED_INSUFFICIENT_COMMON_SUPPORT",
            "poor_n": len(poor), "rich_n": len(rich), "thresholds": [float(low), float(high)],
            "matched_pairs": pairs, "matched_difference": distribution([p["rich_minus_poor_d2"] for p in pairs]),
            "unmatched_difference_of_medians": float(np.median([r["d2"] for r in rich]) -
                                                     np.median([r["d2"] for r in poor]))}


def summarize(output, preparation, baseline):
    assert output["schema_version"] == preparation["schema_version"] == 2
    assert output["data_version"] == preparation["source_commit"]
    diagnostics = {(x["bird"], x["state"], x["variant"]): x for x in output["diagnostics"]}
    windows = defaultdict(list)
    for row in output["windows"]:
        state = preparation["states"][row["bird"]][row["state"]]
        feature, = [x for x in state["window_features"] if x["start_s"] == row["start_s"]]
        windows[row["bird"], row["state"]].append({**row, **feature})
    output["windows"] = [x for rows in windows.values() for x in rows]
    output["state_metrics"] = {}
    output["contrasts"] = {}
    output["sws_association"] = {}
    for bird, states in preparation["states"].items():
        status = output["strict_primary_status"]
        output["state_metrics"][bird] = {
            state: {k: info[k] for k in ("seconds", "units", "strict_active", "strict_status", "burst_summary", "paper99_burst_summary")}
            for state, info in states.items()}
        current = {}
        by_control = {(state, kind): {r["iteration"]: r for r in output["resampling"]
                                    if r["bird"] == bird and r["state"] == state and r["variant"] == kind}
                      for state in STATES for kind in ("unit_subsampling", "rate_thinning")}
        for (state, kind), values in by_control.items():
            assert list(values) == list(range(1, 101)), (bird, state, kind)
            assert len({r["seed"] for r in values.values()}) == 100
            rows = [r for r in values.values() if valid(r)]
            base = diagnostics[bird, state, "broad"]
            output["state_metrics"][bird][state][kind] = {
                "d2": distribution([r["d2"] for r in rows]),
                "rate_hz": distribution([r["actual_rate_hz"] for r in rows]),
                "count_variance_change": distribution([r["count_variance"] - base["count_variance"] for r in rows]),
                "lag1_autocorrelation_change": distribution([r["acf"][0] - base["observed_acf"][0] for r in rows])}
        for hi, lo in PAIRS:
            outcomes = {"strict": {"status": status}, "broad": diagnostics[bird, hi, "broad"]["d2"] -
                        diagnostics[bird, lo, "broad"]["d2"]}
            for kind in ("reach_0p40", "unmasked", "literal_mask"):
                h, l = diagnostics.get((bird, hi, kind)), diagnostics.get((bird, lo, kind))
                if kind in ("unmasked", "literal_mask") and hi == "SWS":
                    h = diagnostics[bird, hi, "broad"]
                outcomes[kind] = h["d2"] - l["d2"] if h and l and valid(h) and valid(l) else None
            for kind in ("unit_subsampling", "rate_thinning"):
                h, l = by_control[hi, kind], by_control[lo, kind]
                outcomes[kind] = distribution([h[i]["d2"] - l[i]["d2"] for i in h
                                               if valid(h[i]) and valid(l[i])])
            current[f"{hi}_minus_{lo}"] = outcomes
        output["contrasts"][bird] = current
        rows = windows[bird, "SWS"]
        seed = output["settings"]["seed"] + int(bird[-1])
        output["sws_association"][bird] = {
            "recording": next(x["session"] for x in output["settings"]["recordings"]
                              if x["bird"] == bird and x["state"] == "SWS"),
            "fixed_window_seconds": 60, "within_session": True,
            "correlations": {name: block_association(rows, name, seed) for name in METRICS},
            "matched_burst_rich_minus_poor": match_burst_windows(rows, output["settings"]["reanalysis"]),
            "cross_state_equal_window_matching": {"status": "BLOCKED_ACTIVE_HAS_NO_60S_WINDOWS",
                                                   "counts": {s: len(windows[bird, s]) for s in STATES}}}
    output["historical_reproduction"] = baseline
    output["burst_detector_validation"] = {k: {f: x[f] for f in ("burst_threshold", "paper99_threshold", "checked_centers_s",
        "archived_event_validation", "archived_event_validation_paper99") if f in x}
        for k, x in preparation["recording_qc"].items()}
    return output


if __name__ == "__main__":
    fit, preparation, baseline, output = map(Path, sys.argv[1:])
    revised = summarize(*(json.loads(p.read_text()) for p in (fit, preparation, baseline)))
    output.write_text(json.dumps(revised, indent=2, allow_nan=False) + "\n")
