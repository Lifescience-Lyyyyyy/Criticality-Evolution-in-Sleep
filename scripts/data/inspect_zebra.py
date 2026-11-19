"""Audit release-matched LMAN spike recordings and author state labels."""

import csv
from datetime import datetime, timezone
from hashlib import file_digest
from pathlib import Path
import json

import numpy as np
from scipy.io import loadmat


def fetch(root, relpath, assets):
    record = assets[relpath]
    target = root / record["drive_path"]
    staged = target.with_suffix(target.suffix + ".part")
    path = target if target.exists() else staged
    with path.open("rb") as source_file:
        digest = file_digest(source_file, "sha256").hexdigest()
    if digest != record["sha256"] or path.stat().st_size != record["bytes"]:
        raise ValueError(f"Drive transfer digest mismatch: {relpath}")
    if path == staged:
        staged.replace(target)
    return target


def lman_and_states(root, item, params, settings, assets):
    bird, session = item["bird"], item["session"]
    base = f"dataset/sub-{bird}/ses-{session}"
    derived = base + "/derived/"
    meta_path = fetch(root, base + "/ephys/" + item["ap_meta"], assets)
    meta = dict(line.split("=", 1) for line in meta_path.read_text().splitlines() if "=" in line)
    rate, duration = float(meta["imSampRate"]), float(meta["fileTimeSecs"])
    scores = fetch(root, derived + item["score"], assets)
    with scores.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if abs(len(rows) - duration) > 2:
        raise ValueError(f"State and electrophysiology durations disagree for {bird} {session}")
    if item["state"] == "WAKE":
        if [int(float(row["sec"])) for row in rows] != list(range(len(rows))):
            raise ValueError(f"Nonconsecutive awake scorer seconds for {bird} {session}")
        valid = [int(row["soundDetected"]) == settings["wake_sound_code"] for row in rows]
    else:
        valid = [int(float(row["final_states_new"])) == settings["sws_code"]
                 and int(float(row["new_mic"])) == 0 for row in rows]
    valid = np.asarray(valid[:int(duration)], dtype=bool)
    bounds = np.flatnonzero(np.diff(np.r_[False, valid, False])).reshape(-1, 2)

    cluster_path = fetch(root, derived + item["prefix"] + "cluster_info.tsv", assets)
    times_path = fetch(root, derived + item["prefix"] + "spike_times.npy", assets)
    labels_path = fetch(root, derived + item["prefix"] + "spike_clusters.npy", assets)
    times = np.load(times_path, allow_pickle=False).reshape(-1).astype(np.int64)
    labels = np.load(labels_path, allow_pickle=False).reshape(-1)
    if (len(times) != len(labels) or not len(times) or np.any(np.diff(times) < 0)
            or times[0] < 0 or times[-1] > (duration + 1) * rate):
        raise ValueError(f"AP spike timeline or sorting arrays failed QC for {bird} {session}")

    low, high = np.asarray(params[bird]["lman"]).astype(int)
    with cluster_path.open(newline="") as stream:
        table = list(csv.DictReader(stream, delimiter="\t"))
    units = [int(row["cluster_id"] if "cluster_id" in row else row["id"]) for row in table
             if low <= int(float(row["ch"])) <= high
             and row["group"].strip().lower() != "noise" and float(row["fr"]) > 0]
    select = np.isin(labels, units)
    lman_times, lman_labels = times[select], labels[select]
    order = np.lexsort((lman_labels, lman_times))
    lman_times, lman_labels = lman_times[order], lman_labels[order]
    duplicates = (lman_times[1:] == lman_times[:-1]) & (lman_labels[1:] == lman_labels[:-1])
    keep = np.r_[True, ~duplicates]
    lman_times, lman_labels = lman_times[keep], lman_labels[keep]
    duplicate_count = int(duplicates.sum())
    unit_count = int(len(np.unique(lman_labels)))
    good = unit_count >= 2 and int(valid.sum()) >= settings["minimum_state_seconds"]
    metrics = {**item, "clock_hz": rate, "ap_duration_s": duration,
               "ap_creation_time": meta["fileCreateTime"], "score_rows": len(rows),
               "probe_map": meta.get("~imroTbl", ""),
               "eligible_seconds": int(valid.sum()), "continuous_intervals": len(bounds),
               "lman_channel_range": [int(low), int(high)], "lman_units": unit_count,
               "lman_spikes": len(lman_times), "all_spikes": len(times),
               "last_spike_s": float(times[-1] / rate),
               "longest_all_unit_interspike_gap_s": float(np.diff(times).max() / rate),
               "duplicate_sample_unit_pairs_removed": duplicate_count,
               "eligible": good}
    return metrics, lman_times, lman_labels, bounds, rate


def audit():
    paths = json.loads(Path("config/paths.json").read_text())
    settings = json.loads(Path("config/zebra_finch.json").read_text())
    source = json.loads(Path(paths["source_manifest"]).read_text())["zebra_finch"]
    root = Path(paths["drive_root"])
    transfer = json.loads((root / "manifests/zebra-finch-transfer.json").read_text())
    actual, assets = transfer["version"], transfer["files"]
    assert actual == source["commit"], f"Transferred version {actual} is outside this audit."
    param_file = fetch(root, "code/UtilsMatlab/birdParams.mat", assets)
    params = loadmat(param_file, simplify_cells=True)["params"]
    recordings = []
    for item in settings["recordings"]:
        metrics, _, _, _, _ = lman_and_states(root, item, params, settings, assets)
        recordings.append(metrics)
        print(json.dumps({key: value for key, value in metrics.items() if key != "probe_map"}), flush=True)

    result = {"source_commit": actual, "verified_at_utc": datetime.now(timezone.utc).isoformat(),
              "parameters": settings, "recordings": recordings, "files": assets}
    path = root / "manifests" / "zebra-finch-audit.json"
    staged = path.with_suffix(".tmp")
    staged.write_text(json.dumps(result, indent=2) + "\n")
    staged.replace(path)
    print(f"Wrote {len(assets)} checksummed files to the zebra finch manifest.")


if __name__ == "__main__":
    audit()
