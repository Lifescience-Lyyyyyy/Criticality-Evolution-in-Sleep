#!/usr/bin/env python3
"""Independently inspect selected published NWB tables and prepare a small audit fixture."""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import argparse
import hashlib
import json
from pathlib import Path
import sys
from time import sleep

import h5py
import numpy as np
from scipy.io import loadmat

from build_atlas_data import ROOT, read, write
from catalog import FIT_FILES, V1_FILE

sys.path.insert(0, str(ROOT / "scripts/data"))
from extract_sleep_spikes import DANDI, fetch_json, labeled_intervals, open_asset
from prepare_phase10 import scores
from prepare_sleep_dandi import matlab_segments, save_prepared

METRICS = ("cell_type", "putative_cell_type", "quality", "lratio", "mahalanobis_distance",
           "isi_index", "in_lab_cell_list")
SELECTION = "results/validation/selected_sessions.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def refs_by_asset():
    out = defaultdict(list)
    for name in FIT_FILES:
        rows = read(name)["records"]
        for row in rows if isinstance(rows, list) else [rows]:
            m = row["record"]
            if m["variant"] == "reference" and m["context"] in {"HOME", "REST"}:
                out[m["asset_id"]].append(m)
    return out


def clock_groups(file):
    clocks = []
    for name in ("acquisition", "processing/ecephys/LFP"):
        obj = file.get(name)
        if not isinstance(obj, h5py.Group):
            continue

        def clock(_, node):
            if not isinstance(node, h5py.Group) or "data" not in node:
                return
            if "starting_time" in node:
                start = node["starting_time"]
                t0, rate = float(start[()]), float(start.attrs["rate"])
                clocks.append(dict(path=node.name, start=t0, end=t0 + node["data"].shape[0] / rate,
                                   samples=node["data"].shape[0], rate=rate))
        obj.visititems(clock)
    return clocks


def audit_nwb(case, refs, version):
    dataset, asset = case["dandiset"], case["asset_id"]
    url = f"{DANDI}/{dataset}/versions/{version}/assets/{asset}/"
    for attempt in range(3):
        try:
            meta = fetch_json(url)
            with open_asset(meta) as f:
                return extract(f, case, refs, meta, url)
        except (OSError, TimeoutError):
            if attempt == 2:
                raise
            sleep(2 * (attempt + 1))


def extract(f, case, refs, meta, url):
    print(case["case"], "opened NWB", flush=True)
    u = f["units"]
    ids = [str(i) for i in u["id"][:]]
    columns = {}
    for k in METRICS:
        if k not in u:
            continue
        d = u[k]
        values = (d.asstr()[:] if d.dtype.kind in "OSU" else d[:]).tolist()
        values = [v if not isinstance(v, float) or np.isfinite(v) else None for v in values]
        columns[k] = {"description": str(d.attrs.get("description", "")), "values": dict(zip(ids, values, strict=True))}
    clocks = clock_groups(f)
    indexes = np.r_[0, u["spike_times_index"][:]]
    all_spikes = u["spike_times"][:] if case["case"] in ("A", "B", "D") else None
    stop_spikes = [float(all_spikes[int(i)-1]) for i in set(indexes[1:]) if i] if case["case"] == "D" else []
    stages = {}
    for name, column in (("processing/behavior/sleep_states_buzcode", "label"),
                         ("processing/behavior/sleep_states", "label"),
                         ("intervals/sleep_states", "state"), ("processing/behavior/SleepStates", "label")):
        if name in f:
            stages[name] = scores(f[name], column)
    epochs = labeled_intervals(f["intervals/epochs"], "tags") if "intervals/epochs" in f else []
    subrecords = (labeled_intervals(f["intervals/sub_recordings"], "tags")
                  if "intervals/sub_recordings/tags_index" in f else
                  [[float(a), float(b)] for a, b in zip(f["intervals/sub_recordings/start_time"][:],
                      f["intervals/sub_recordings/stop_time"][:], strict=True)] if "intervals/sub_recordings" in f else [])
    audit = dict(case=case["case"], asset=case["asset_id"], source=url, digest=meta["digest"],
                 clocks=clocks, spikes_last_s=max(stop_spikes, default=None), metrics=columns,
                 unitCount=len(ids), stages=stages, epochs=epochs, subrecordings=subrecords)
    tasks = []
    if case["case"] not in ("A", "B"):
        return audit, tasks
    by_id = {u: i for i, u in enumerate(ids)}
    for region in case["regions"]:
        rows = [r for r in refs if r["region"] == region]
        assert len(rows) == 3
        original = list(map(str, rows[0]["unit_ids"]))
        assert all(r["unit_ids"] == rows[0]["unit_ids"] for r in rows)
        spikes = {s: all_spikes[indexes[by_id[s]]:indexes[by_id[s]+1]] for s in original}
        hasher = hashlib.sha256()
        for s in original:
            hasher.update(s.encode() + b"\0")
            hasher.update(np.asarray(spikes[s], dtype="<f8").tobytes())
        assert hasher.hexdigest() == rows[0]["extracted_unit_sha256"], (case["session"], region)
        lr = columns["lratio"]["values"]
        xs = [lr[s] for s in original if lr[s] is not None]
        assert len(xs) >= 5
        for label, quantile in (("lratio_q75", .75), ("lratio_q50", .5)):
            cutoff = float(np.quantile(xs, quantile))
            units = [s for s in original if lr[s] is not None and lr[s] <= cutoff]
            assert len(units) >= 5
            pooled = np.sort(np.concatenate([spikes[s] for s in units]))
            for r in rows:
                intervals = [(round(a/.04), round(b/.04)) for a, b in r["intervals_s"]]
                tasks.append(({**r, "variant": label, "unit_ids": units, "lratio_cutoff": cutoff,
                               "audit_case": case["case"]}, matlab_segments(pooled, intervals, .04)))
    return audit, tasks


def v1_tasks(case):
    session = case["session"]
    fixture = ROOT / "data/fixtures/mouse_population_40ms.mat"
    prep = read("manifests/mouse-preparation.json")
    assert digest(fixture) == prep["fixture_sha256"]
    population = loadmat(fixture, simplify_cells=True)["cohort"][session]
    source = next(s for s in prep["sessions"] if s["session"] == session)
    fits = [r for r in read(V1_FILE)["primary"] if r["session"] == session]
    for f in fits:
        state = f["state"]
        cells = np.empty((len(population[state]), 1), dtype=object)
        for i, segment in enumerate(population[state]):
            cells[i, 0] = np.asarray(segment)
        assert sum(map(len, cells[:, 0])) == f["samples"]
        meta = dict(dandiset="PNAS-V1", session=session, subject=None, region="V1", context="PUBLISHED",
                    asset_id=source["file"], path=source["file"], digest=source["md5"],
                    state="Wake" if state == "WAKE" else state, variant="reference", bin_seconds=.04,
                    audit_case=case["case"], input_count_fixture_sha256=digest(fixture))
        yield meta, cells


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--receipt", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    selection = read(SELECTION)
    cases = [r for r in selection["cases"] if r["dandiset"].isdecimal()]
    versions = {**read("manifests/sleep_sources.json")["dandisets"],
                **read("manifests/phase10_sources.json")["sources"]}
    references = refs_by_asset()
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(audit_nwb, r, references[r["asset_id"]],
                               versions[r["dandiset"]]["version"]) for r in cases]
        outputs = []
        for f in futures:
            audit, tasks = f.result()
            print(audit["case"], "checked", len(tasks), "fit inputs", flush=True)
            outputs.append((audit, tasks))
    write(a.receipt, [x[0] for x in outputs])
    v1 = next(r for r in selection["cases"] if r["case"] == "I")
    tasks = [task for _, group in outputs for task in group] + list(v1_tasks(v1))
    assert len(tasks) == 20
    path = a.out / "validation_input.mat"
    meta = save_prepared(tasks, path, "states", sorted({m["asset_id"] for m, _ in tasks}),
                         source_sha256={"selection": digest(ROOT / SELECTION), "remote_receipt": digest(a.receipt)})
    print(meta["records"], "fits; input SHA-256", meta["sha256"])


if __name__ == "__main__":
    main()
