"""Convert already audited CRCNS and NWB spike intervals to the shared MATLAB schema."""

from collections import defaultdict
from pathlib import Path
import argparse
import json

import numpy as np

from crcns_adapter import audited_nwb, clip, hc23_prepared, hc24_matlab, neurosuite, sha256
from prepare_sleep_dandi import clock_bins, matlab_segments, save_prepared, whole_windows


def paired_intervals(group):
    """Use the same scored bin support in all co-recorded regions."""
    states = set.intersection(*(set(r.sleep_intervals) for r in group))
    common = {}
    for state in states:
        bounds = group[0].sleep_intervals[state]
        for r in group[1:]:
            bounds = clip(bounds, r.sleep_intervals[state])
        common[state] = bounds
    return common


def quality_sets(record, cfg):
    keep = []
    for unit, times in zip(record.unit_ids, record.spike_times, strict=True):
        count = 0
        short = 0
        for x, y in record.recording_intervals:
            a, b = np.searchsorted(times, (x, y))
            diffs = np.diff(times[a:b])
            count += len(diffs)
            short += int((diffs < cfg["refractory_seconds"]).sum())
        if count and short / count <= cfg["max_refractory_fraction"]:
            keep.append(unit)
    return keep



EQUAL_DURATION_STATES = {
    "equal_duration_NW": ("Wake", "NREM"),
    "equal_duration_RN": ("NREM", "REM"),
    "equal_duration_RW": ("Wake", "REM"),
    "equal_duration_WNR": ("Wake", "NREM", "REM"),
}


def equal_duration_controls(windows):
    """Select first, lower median and last complete 60 s windows per eligible state."""
    for variant, states in EQUAL_DURATION_STATES.items():
        if any(len(windows.get(state, ())) < 3 for state in states):
            continue
        for state in states:
            candidates = windows[state]
            yield variant, state, [candidates[i] for i in (0, (len(candidates)-1)//2, -1)]


def prepare_group(group, cfg, mode):
    states = paired_intervals(group)
    if mode == "nrem":
        states = {"NREM": states["NREM"]} if "NREM" in states else {}
    width = cfg["bin_seconds"]
    if mode == "pairwise" and (width != .04 or cfg["window_seconds"] != 60):
        raise ValueError("Pairwise controls require 40 ms bins and 60 s windows")
    n = round(cfg["window_seconds"] / width)
    effective = {s: clock_bins(ranges, width, cfg["min_segment_seconds"]) for s, ranges in states.items()}
    states = {s: states[s] for s in ("Wake", "NREM", "REM") if s in effective
              and width * sum(b - a for a, b in effective[s]) >= cfg["min_state_seconds"]}
    if (mode != "pairwise" and "NREM" not in states) or (mode != "nrem" and len(states) < 2):
        return []
    windows = {s: whole_windows(effective[s], n) for s in states}
    count = min(map(len, windows.values()))
    exact_control_states = ("Wake", "NREM") if "Wake" in states else ("NREM", "REM")
    if mode == "controls" and min(len(windows[s]) for s in exact_control_states) < 3:
        return []
    minimum_units = min(len(r.unit_ids) for r in group)
    result = []
    for r in group:
        strict = quality_sets(r, cfg)
        meta = {"dandiset": r.source_dataset, "subject": r.subject_id,
                "session": r.session_id, "region": r.region, "context": r.context,
                "asset_id": r.asset_id, "path": r.provenance["path"], "digest": r.provenance["sha256"],
                "clock": r.timebase, "state_path": r.provenance.get("state_source", "validated sidecar"),
                "unit_ids": list(r.unit_ids), "strict_unit_ids": strict, "bin_seconds": width}
        pooled = np.sort(np.concatenate(r.spike_times))
        if mode == "pairwise":
            for variant, state, spans in equal_duration_controls(windows):
                result.append(({**meta, "state": state, "variant": variant,
                                "intervals_s": [[a*width,b*width] for a,b in spans]},
                               matlab_segments(pooled,spans,width)))
            continue
        strict_spikes = (np.sort(np.concatenate([t for u, t in zip(r.unit_ids, r.spike_times, strict=True)
                                                if u in strict]))
                         if cfg["min_units"] <= len(strict) < len(r.unit_ids) else None)
        for s, ranges in states.items():
            state = {**meta, "state": s}
            matched = windows[s]
            if mode == "controls":
                if s not in exact_control_states:
                    continue
                chosen = [matched[i] for i in np.linspace(0, len(matched)-1, 3, dtype=int)]
                cells = np.empty((len(chosen), 1), dtype=object)
                for j, interval in enumerate(chosen):
                    cells[j, 0] = np.stack([matlab_segments(times, [interval], width)[0, 0]
                                             for times in r.spike_times])
                target = minimum_units if len(group) > 1 else max(cfg["min_units"], len(r.unit_ids)//2)
                result.append(({**state, "variant": "matched_controls", "target_units": target,
                                "intervals_s": [[a*width, b*width] for a, b in chosen]}, cells))
                continue

            segments = effective[s]
            variants = [("reference", segments, pooled, r.unit_ids)]
            trimmed = clock_bins(ranges, width, cfg["min_segment_seconds"], cfg["edge_trim_seconds"])
            if width * sum(b-a for a, b in trimmed) >= cfg["min_state_seconds"]:
                variants.append(("edge_trim", trimmed, pooled, r.unit_ids))
            if strict_spikes is not None:
                variants.append(("strict_units", segments, strict_spikes, strict))
            if mode != "nrem" and count >= 3:
                chosen = [matched[i] for i in np.linspace(0, len(matched)-1, count, dtype=int)]
                variants.append(("equal_duration", chosen, pooled, r.unit_ids))
            for variant, spans, spikes, selected in variants:
                result.append(({**state, "variant": variant, "unit_ids": list(selected),
                                "intervals_s": [[a*width, b*width] for a,b in spans]},
                               matlab_segments(spikes, spans, width)))
    return result


def load_specs(path, cfg):
    for row in json.loads(path.read_text()):
        if row["adapter"] == "nwb":
            yield audited_nwb(row["audit_row"], row["region"], row["context"], cfg)
        else:
            directory = path.parent / row["directory"]
            if row["adapter"] == "hc24":
                yield hc24_matlab(directory, row)
            elif row["adapter"] == "neurosuite":
                yield neurosuite(directory, row)
            else:
                raise ValueError(f"Unknown adapter: {row['adapter']}")


def build(records, path, cfg, mode, source_hashes):
    grouped = defaultdict(list)
    for r in records:
        if len(r.unit_ids) >= cfg["min_units"]:
            grouped[(r.source_dataset, r.subject_id, r.session_id, r.context)].append(r)
    if mode == "controls":
        best = {}
        for (dataset, subject, session, context), group in sorted(grouped.items()):
            rows = prepare_group(group, cfg, mode)
            if rows:
                score = min(len(r.unit_ids) for r in group), session
                key = dataset, subject, context
                if key not in best or score > best[key][0]:
                    best[key] = score, rows
        data = [job for _, rows in best.values() for job in rows]
    else:
        data = [job for group in grouped.values() for job in prepare_group(group, cfg, mode)]
    return save_prepared(data, path, mode, sorted({m["asset_id"] for m, _ in data}),
                         source_sha256=source_hashes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument("--hc23-audit", type=Path)
    sources.add_argument("--verified-index", type=Path)
    parser.add_argument("--hc23-mat", type=Path)
    parser.add_argument("--hc23-overnight", type=Path)
    parser.add_argument("--asset")
    parser.add_argument("--mode", choices=("states","controls","pairwise"), default="states")
    parser.add_argument("--config", type=Path, default=Path("config/sleep_dandi.json"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    if args.hc23_audit:
        if not args.hc23_mat or not args.hc23_overnight:
            parser.error("Provide the verified MAT and sorted-spike cache")
        audit = json.loads(args.hc23_audit.read_text())
        hashes = {"author_mat": sha256(args.hc23_mat), "sorted_t_cache": sha256(args.hc23_overnight)}
        if hashes["sorted_t_cache"] != audit["source"]["cache_sha256"]:
            raise ValueError("Sorted-spike cache differs from audited input")
        eligible = {r["asset_id"]: r for r in audit["records"] if r["eligible_comparisons"]}
        records = []
        for r in hc23_prepared(args.hc23_mat, args.hc23_overnight):
            if r.asset_id in eligible and (not args.asset or args.asset == r.asset_id):
                if len(r.unit_ids) != eligible[r.asset_id]["units"] or r.quality_metrics != eligible[r.asset_id]["quality"]:
                    raise ValueError(f"Source QC no longer matches eligibility audit: {r.asset_id}")
                records.append(r)
    else:
        records = [r for r in load_specs(args.verified_index,cfg) if not args.asset or r.asset_id == args.asset]
        hashes = {"verified_index": sha256(args.verified_index)}
    print(build(records,args.out,cfg,args.mode,hashes))
