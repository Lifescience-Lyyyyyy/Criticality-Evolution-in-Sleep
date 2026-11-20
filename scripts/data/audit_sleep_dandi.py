"""Audit published DANDI assets by reading NWB annotations and sorted spike tables only."""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from extract_sleep_spikes import (DANDI, STATES, bounds, fetch_json, labeled_intervals, listing, open_asset,
                                  read_units, recording_clock, ripple_times, scalar_text, state_table, version_path)


def subtract(intervals, excluded):
    for x, y in sorted(excluded):
        intervals = [(s, t) for a, b in intervals
                     for s, t in ((a, min(b, x)), (max(a, y), b)) if t > s]
    return intervals


def intersect(a, b):
    return [(max(x, u), min(y, v)) for x, y in a for u, v in b if min(y, v) > max(x, u)]


def label_states(rows):
    states = {"Wake": [], "NREM": [], "REM": []}
    ordered = sorted(rows)
    for i, (a, b, raw) in enumerate(ordered):
        if i > 0 and a < ordered[i-1][1] - 0.001:
            raise ValueError(f"State overlap: {ordered[i-1]} and {(a,b,raw)}")
        canonical = STATES.get(raw, STATES.get(raw.lower()))
        if canonical:
            states[canonical].append((max(a, ordered[i-1][1]) if i else a, b))
    return states


def fields(file, row, cfg):
    ds = row["dandiset"]
    desc = scalar_text(file["session_description"])
    row.update(subject=scalar_text(file["general/subject/subject_id"]),
               session=scalar_text(file["general/session_id"]), description=desc,
               nwb_id=scalar_text(file["identifier"]),
               session_start_time=scalar_text(file["session_start_time"]))
    path, states_raw = state_table(file)
    row.update(state_path=path, state_labels=sorted({r[-1] for r in states_raw}),
               state_intervals_source=states_raw)
    groups, trains = read_units(file, ds, cfg, load_times=bool(states_raw))
    row["units_total"] = len(groups)
    reasons = []
    if not states_raw:
        reasons.append("No readable sleep-state table")
    if not groups:
        reasons.append("No sortable NWB spike times")
    if reasons:
        return {**row, "status": "BLOCKED", "reason": "; ".join(reasons)}, []
    row["units"] = groups
    clock = recording_clock(file)
    if clock is None:
        return {**row, "status": "PARTIAL", "reason": "Recording boundary unavailable separately from states/spikes"}, []
    row["clock"] = clock
    start, stop = clock["seconds"]
    if not stop > start:
        raise ValueError("Nonpositive recording length")
    if any(r["first"] is not None and r["first"] < start-1 or
           r["last"] is not None and r["last"] > stop+1 for r in groups):
        raise ValueError("Unit spike times outside independent recording clock by >1 s")
    ordered_spikes = np.sort(np.concatenate([x for x, u in zip(trains, groups, strict=True) if not u["nonmonotonic"]]))
    if not len(ordered_spikes):
        raise ValueError("No finite sorted spike trains")
    row["population_spike_range_s"] = [float(ordered_spikes[0]), float(ordered_spikes[-1])]
    spike_coverage = np.diff(row["population_spike_range_s"])[0]
    if spike_coverage < (stop-start) * 0.5:
        raise ValueError("Spike span < half the LFP clock; possible conversion/coverage mismatch")
    breaks = np.flatnonzero(np.diff(ordered_spikes) > 10)
    gaps = [(float(ordered_spikes[i]), float(ordered_spikes[i+1])) for i in breaks]
    row["global_spike_gaps_over_10s"] = gaps
    row["duplicate_spikes"] = sum(u["duplicates"] for u in groups)
    row["nonmonotonic_spikes"] = sum(u["nonmonotonic"] for u in groups)
    row["unit_location_labels"] = dict(Counter(u["region_raw"] for u in groups))
    row["spontaneous_contexts"] = {}
    row["eligible_regions"] = {key: c for key, c in Counter(u["region"] for u in groups if u["eligible"]).items()
                               if c >= cfg["min_units"]}
    row["strict_units"] = dict(Counter(u["region"] for u in groups if u["strict"]))
    state_ranges = label_states(states_raw)
    if max(r[1] for r in states_raw) > stop+2 or min(r[0] for r in states_raw) < start-2:
        raise ValueError("Stage intervals outside independent LFP clock by >2 s")
    # Clipping to population observation bounds avoids treating unobserved edges as silence.
    coverage = [(max(start, ordered_spikes[0]), min(stop, ordered_spikes[-1]))]
    independent = subtract(coverage, gaps)
    context = {"HOME": independent}
    if ds == "000061":
        if "intervals/epochs" not in file:
            return {**row,"status":"BLOCKED","reason":"PRE/POST home-cage epochs missing"}, []
        ep = labeled_intervals(file["intervals/epochs"], "tags", 1000)
        if not (max(r[1] for r in ep) <= stop + 2 and max(r[1] for r in ep) >= 0.75*stop):
            raise ValueError("Epochs milliseconds-to-seconds conversion fails recording-clock reconciliation")
        row["epoch_scale_raw_to_seconds"] = 0.001
        row["epoch_intervals"] = ep
        spontaneous = {"PRE": "presleep", "POST": "postsleep"}
        context = {k: intersect(independent, [(a,b) for a,b,raw in ep if raw.lower() == v])
                   for k,v in spontaneous.items()}
        task = [(a,b) for a,b,raw in ep if raw.lower() not in spontaneous.values()]
        context = {k: subtract(t,task) for k,t in context.items()}
    elif ds == "000552":
        if "home cage" not in desc.lower() or any(x in desc.lower() for x in ("maze", "linear track", "novel environment")):
            return {**row,"status":"PARTIAL","reason":"Wake environment/task boundaries cannot be determined as home cage"}, []
        p = "intervals/OptogeneticPulsePeriods"
        if p in file or "intervals/OptogeneticPulses" in file:
            p = p if p in file else "intervals/OptogeneticPulses"
            emitted = bounds(file[p])
            d = cfg["opto_guard_seconds"]
            row["opto_excluded_table"] = p
            row["opto_excluded_intervals"] = [[float(a)-d,float(b)+d] for a,b in emitted]
            context = {"HOME": subtract(independent, row["opto_excluded_intervals"])}
        elif "opto" in desc.lower() or "light pulses" in desc.lower():
            return {**row,"status":"BLOCKED","reason":"Optogenetic stimulation declared without an interval table"}, []
    eligible = False
    for name, usable in context.items():
        cells = {}
        for state, times in state_ranges.items():
            segments = [(float(a),float(b)) for a,b in intersect(usable,times) if b-a >= cfg["min_segment_seconds"]]
            cells[state] = {"intervals": segments, "seconds": sum(b-a for a,b in segments), "segments":len(segments)}
        row["spontaneous_contexts"][name] = cells
        eligible |= (bool(row["eligible_regions"]) and cells["NREM"]["seconds"] >= cfg["min_state_seconds"]
                     and cells["Wake"]["seconds"] >= cfg["min_state_seconds"])
    try:
        path, t = ripple_times(file)
        if len(t):
            if t[:,0].min() < start-1 or t[:,2].max() > stop+1:
                raise ValueError("Ripple clock outside recording")
            row["ripples"] = {"path":path, "peak_count":len(t), "duration_available":True}
    except (ValueError, KeyError) as exc:
        row["event_reason"] = str(exc)
    sufficient_sleep = any(c["NREM"]["seconds"] >= cfg["min_state_seconds"] for c in row["spontaneous_contexts"].values())
    if eligible:
        status, reason = "READY", "Paired author Wake/NREM and sorted regional spikes on a verified continuous recording clock"
    elif sufficient_sleep and row["eligible_regions"]:
        status, reason = "PARTIAL", "NREM qualifies; no same-context Wake duration for primary contrast"
    else:
        status, reason = "BLOCKED", "Insufficient spontaneous state duration or anatomical sorting QC"
    return {**row,"status":status,"reason":reason}, trains


def audit_asset(row, cfg):
    try:
        base = version_path(row["dandiset"], row["version"])
        meta = fetch_json(f"{base}/assets/{row['asset_id']}/")
        row["digest"] = meta["digest"]
        row["access"] = meta.get("access", [])
        row["metadata_url"] = f"{base}/assets/{row['asset_id']}/"
        row["download_url"] = meta["contentUrl"][0]
        with open_asset(meta) as file:
            row["root_groups"] = list(file)
            return fields(file, row, cfg)[0]
    except (ValueError, KeyError, OSError, IndexError) as error:
        public_access_failed = hasattr(error, "code") and error.code in (401, 403)
        return {**row,"status":"RESTRICTED" if public_access_failed else "BLOCKED",
                "reason":f"{type(error).__name__}: {error}"}


def audit(out, cfg, workers, filter_datasets=None):
    out.mkdir(parents=True, exist_ok=True)
    metadata = {}
    candidates = []
    for ds, version in cfg["versions"].items():
        if filter_datasets and ds not in filter_datasets:
            continue
        url = version_path(ds, version)
        version_meta = fetch_json(url+"/")
        public = fetch_json(f"{DANDI}/{ds}/")
        versions = [v["version"] for v in listing(f"{DANDI}/{ds}/versions/?page_size=100")]
        metadata[ds] = {"version": version,"published_versions": [v for v in versions if v != "draft"],
                        "access":public.get("embargo_status"), "title":version_meta.get("name"),
                        "license":version_meta.get("license", []), "dandiset_version_url":url+"/"}
        for a in listing(url+"/assets/?page_size=200"):
            if a["path"].endswith(".nwb"):
                candidates.append({"dandiset":ds, "version":version,"asset_id":a["asset_id"],"path":a["path"],
                                   "bytes":a["size"],"subject":a["path"].split("/")[0].removeprefix("sub-"),
                                   "session":a["path"].split("_ses-")[-1].split("_")[0]})
    receipts = out / "asset_receipts"
    receipts.mkdir(exist_ok=True)
    records = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        jobs = {}
        for row in candidates:
            path = receipts / (row["asset_id"]+".json")
            if path.exists():
                records.append(json.loads(path.read_text()))
            else:
                jobs[executor.submit(audit_asset, row, cfg)] = path
        for future in as_completed(jobs):
            record = future.result()
            record["audit_settings"] = {k:cfg[k] for k in ("min_segment_seconds","min_state_seconds","min_units",
              "opto_guard_seconds","refractory_seconds","max_refractory_fraction")}
            target = jobs[future]
            partial = target.with_suffix(".part")
            partial.write_text(json.dumps(record, allow_nan=False))
            partial.replace(target)
            records.append(record)
            print(f"{len(records)}/{len(candidates)} {record['dandiset']} {record['session']} "
                  f"{record['status']} {record['reason']}", flush=True)
    records.sort(key=lambda r:(r["dandiset"],r["path"]))
    return metadata,records


def export(out, metadata, records):
    (out / "sleep_sources.json").write_text(json.dumps(
        {"schema_version":1,"audited_utc":datetime.now(timezone.utc).isoformat(),
         "dandisets":metadata,"assets":[{k:v for k,v in r.items() if k not in
           ("units","state_intervals_source","opto_excluded_intervals","spontaneous_contexts")} for r in records]},
        indent=2,allow_nan=False)+"\n")
    keys = ("dandiset","version","asset_id","path","bytes","subject","session","status","reason","state_path",
            "units_total","state_labels","eligible_regions","strict_units","contexts","ripples","recording_seconds")
    with (out / "sleep_dataset_eligibility.csv").open("w",newline="") as handle:
        writer = csv.DictWriter(handle,fieldnames=keys,extrasaction="ignore")
        writer.writeheader()
        for row in records:
            slim = dict(row)
            for k in ("state_labels","eligible_regions","strict_units","ripples"):
                if k in row:
                    slim[k] = json.dumps(row[k], separators=(",",":"))
            slim["contexts"] = json.dumps({p:{s:round(e["seconds"],3) for s,e in q.items()}
                                     for p,q in row.get("spontaneous_contexts",{}).items()},separators=(",",":"))
            if "clock" in row:
                lo, hi = row["clock"]["seconds"]
                slim["recording_seconds"] = round(hi-lo, 3)
            writer.writerow(slim)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=Path("config/sleep_dandi.json"))
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--workers",type=int,default=6)
    parser.add_argument("--dandiset",action="append")
    o = parser.parse_args()
    export(o.out,*audit(o.out,json.loads(o.config.read_text()),o.workers,o.dandiset))
