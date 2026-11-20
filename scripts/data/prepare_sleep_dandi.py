"""Build versioned segment-aware MATLAB count fixtures from audited DANDI sessions."""

from concurrent.futures import ProcessPoolExecutor
from hashlib import file_digest
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.io import savemat

from extract_sleep_spikes import fetch_json, open_asset, read_units, ripple_times


def clock_bins(intervals, dt, min_seconds=1, trim=0):
    output = []
    for start, stop in intervals:
        a = int(np.ceil((start+trim)/dt - 1e-8))
        b = int(np.floor((stop-trim)/dt + 1e-8))
        if (b-a)*dt >= min_seconds-1e-10:
            output.append((a, b))
    return output


def matlab_segments(spikes, segments, width):
    cells = np.empty((len(segments), 1), dtype=object)
    for i, (a, b) in enumerate(segments):
        left = np.searchsorted(spikes, a*width-1e-10)
        right = np.searchsorted(spikes, b*width-1e-10)
        ticks = np.floor(spikes[left:right]/width+1e-9).astype(np.int64) - a
        if np.any(ticks < 0) or np.any(ticks >= b-a):
            raise ValueError("Spike-to-bin endpoint disagreement")
        counts = np.bincount(ticks, minlength=b-a)
        if np.max(counts, initial=0) > np.iinfo(np.uint16).max:
            raise ValueError("Population count overflow")
        cells[i, 0] = counts.astype(np.uint16)
    return cells


def whole_windows(segments, n_bins):
    return [(start, start+n_bins) for a,b in segments for start in range(a,b-n_bins+1,n_bins)]


def event_features(events,start,stop):
    starts, _, ends = events.T
    keep = (starts < stop) & (ends > start)
    # Intersections of already validated author events; merge any overlap before occupancy.
    edges = sorted((max(a,start),min(b,stop)) for a,b in zip(starts[keep],ends[keep],strict=True))
    seconds = 0.
    previous = start
    for a,b in edges:
        seconds += max(0,b-max(a,previous))
        previous = max(previous,b)
    n = int(np.sum((events[:,1]>=start) & (events[:,1]<stop)))
    return {"swr_count":n, "swr_rate_hz":n/(stop-start), "swr_occupancy":seconds/(stop-start)}


def prepare(job):
    row, cfg, events_mode = job
    with open_asset(fetch_json(row["metadata_url"])) as f:
        qc, trains = read_units(f,row["dandiset"],cfg)
        required = [k for k in ("id","region","eligible","strict","spikes","duplicates","nonmonotonic")]
        signature = lambda cells:[tuple(u[k] for k in required) for u in cells]
        if signature(qc) != signature(row["units"]):
            raise ValueError(f"Unit data changed since audit: {row['asset_id']}")
        events = np.empty((0,3))
        if events_mode and row.get("ripples"):
            _, events = ripple_times(f)
            if len(events) != row["ripples"]["peak_count"]:
                raise ValueError("Ripple count changed after audit")
    width = cfg["bin_seconds"]
    wn = round(cfg["window_seconds"]/width)
    n_min = cfg["min_state_seconds"]/width
    data = []
    for region in sorted(row["eligible_regions"]):
        unit_idx = [i for i,u in enumerate(qc) if u["eligible"] and u["region"]==region]
        strict_idx = [i for i in unit_idx if qc[i]["strict"]]
        pooled = np.sort(np.concatenate([trains[i] for i in unit_idx]))
        strict_spikes = np.sort(np.concatenate([trains[i] for i in strict_idx])) if len(strict_idx)>=cfg["min_units"] else None
        provenance = {k:row[k] for k in ("dandiset","version","asset_id","path","digest","subject","session",
                                        "state_path","state_labels","clock","nwb_id")}
        for context, states in row["spontaneous_contexts"].items():
            bins = {s:clock_bins(x["intervals"],width,cfg["min_segment_seconds"]) for s,x in states.items()}
            states_to_fit = [s for s,ranges in bins.items() if sum(b-a for a,b in ranges)>=n_min-1e-8]
            # All primary comparisons use eligible within-context Wake and NREM.
            if not {"Wake","NREM"}<=set(states_to_fit):
                continue
            windows = {s:whole_windows(bins[s],wn) for s in states_to_fit}
            least = min(len(windows[s]) for s in ("Wake","NREM"))
            for state in states_to_fit:
                if events_mode and (state!="NREM" or region!="CA1" or not len(events)):
                    continue
                base = bins[state]
                meta = {**provenance,"region":region,"context":context,"state":state,
                        "bin_seconds":width,"unit_ids":[qc[i]["id"] for i in unit_idx],
                        "strict_unit_ids":[qc[i]["id"] for i in strict_idx]}
                if events_mode:
                    for start, stop in windows[state]:
                        start_s,stop_s = start*width,stop*width
                        m = {**meta,"variant":"swr_window","start_s":start_s,**event_features(events,start_s,stop_s)}
                        data.append((m,matlab_segments(pooled,[(start,stop)],width)))
                    continue
                data.append(({**meta,"variant":"reference","intervals_s":[[a*width,b*width] for a,b in base]},
                             matlab_segments(pooled,base,width)))
                trimmed = clock_bins(states[state]["intervals"],width,cfg["min_segment_seconds"],cfg["edge_trim_seconds"])
                if sum(b-a for a,b in trimmed)>=n_min-1e-8:
                    data.append(({**meta,"variant":"edge_trim"},matlab_segments(pooled,trimmed,width)))
                if strict_spikes is not None and len(strict_idx)<len(unit_idx):
                    data.append(({**meta,"variant":"strict_units","unit_ids":meta["strict_unit_ids"]},
                                 matlab_segments(strict_spikes,base,width)))
                if state in ("Wake","NREM") and least>=3:
                    origins = windows[state]
                    ids = np.linspace(0,len(origins)-1,least,dtype=int)
                    matched = [origins[i] for i in ids]
                    data.append(({**meta,"variant":"equal_duration"},matlab_segments(pooled,matched,width)))
    print(row["dandiset"],row["session"],len(data),flush=True)
    return data


def write_fixture(rows,cfg,path,workers,events_mode):
    jobs = ((r,cfg,events_mode) for r in rows)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        data = [item for results in pool.map(prepare,jobs) for item in results]
    if not data:
        raise ValueError("No eligible fits survived audited preparation")
    path.parent.mkdir(parents=True,exist_ok=True)
    saved = {f"r{i:05d}":{"metadata":json.dumps(meta,allow_nan=False),"counts":counts}
             for i,(meta,counts) in enumerate(data,1)}
    temp = path.with_suffix(".part")
    savemat(temp,{"cohort":saved},appendmat=False,do_compression=True,long_field_names=True)
    with temp.open("rb") as stream:
        digest=file_digest(stream,"sha256").hexdigest()
    if temp.stat().st_size>=25*2**20:
        raise ValueError(f"Fixture {temp.stat().st_size} bytes requires a separate transfer batch")
    receipt = {"fixture":path.name,"sha256":digest,"bytes":temp.stat().st_size,"records":len(data),
               "asset_ids":[r["asset_id"] for r in rows],
               "reference_fits":sum(m["variant"]=="reference" for m,_ in data),"swr_windows":sum(m["variant"]=="swr_window" for m,_ in data)}
    temp.replace(path)
    path.with_suffix(".json").write_text(json.dumps(receipt,indent=2)+"\n")
    return receipt


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--dandiset",required=True)
    p.add_argument("--asset")
    p.add_argument("--config",type=Path,default=Path("config/sleep_dandi.json"))
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--paired-regions",action="store_true")
    p.add_argument("--events",action="store_true")
    a=p.parse_args()
    cfg=json.loads(a.config.read_text())
    rows=[json.loads(f.read_text()) for f in sorted((a.audit/"asset_receipts").glob("*.json"))]
    rows=[r for r in rows if r["dandiset"]==a.dandiset and r["status"]=="READY"
          and (not a.asset or r["asset_id"]==a.asset)
          and (not a.paired_regions or {"CA1","BLA"}<=set(r["eligible_regions"]))]
    if a.asset and len(rows)!=1:
        raise ValueError("Pilot input must match one audited READY asset")
    print(write_fixture(rows,cfg,a.out,a.workers,a.events))
