"""Prepare the validated three-state LMAN data and burst evidence on Colab."""

from collections import Counter
import csv
from datetime import datetime, timezone
from hashlib import file_digest
from pathlib import Path
import json
from tempfile import TemporaryDirectory

import numpy as np
from scipy.io import loadmat, savemat

from inspect_zebra import fetch, lman_and_states
from zebra_bursts import detect, summarize


STATES = ("ACTIVE", "NON_ACTIVE", "SWS")
BINS_PER_SECOND = 25


def read_table(path, delimiter=","):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def continuous(mask):
    return np.flatnonzero(np.diff(np.r_[False, mask, False])).reshape(-1, 2)


def segment_matrix(matrix, intervals):
    result = np.empty((len(intervals), 1), dtype=object)
    for i, (start, stop) in enumerate(intervals):
        result[i, 0] = matrix[:, start * BINS_PER_SECOND:stop * BINS_PER_SECOND]
    return result


def trace_summary(events, source):
    windows = loadmat(source, variable_names=["locs", "sstate", "t", "lfsb", "lfrnd"],
                      squeeze_me=True)
    ref_states = windows["sstate"] == 3
    ref_peaks = windows["locs"][ref_states]
    peaks = np.array([e["peak_s"] for e in events])
    idx = np.searchsorted(peaks, ref_peaks)
    errors = np.minimum(abs(peaks[np.clip(idx - 1, 0, len(peaks) - 1)] - ref_peaks),
                        abs(peaks[np.clip(idx, 0, len(peaks) - 1)] - ref_peaks))
    t = windows["t"]
    center = abs(t) <= 0.025
    base = (t >= -0.20) & (t <= -0.10)
    def contrast(key):
        x = windows[key][ref_states]
        return x[:, center].mean(axis=1) - x[:, base].mean(axis=1)
    matched = int(np.count_nonzero(errors <= 0.02000001))
    if matched < min(10, len(ref_peaks)):
        raise ValueError(f"Insufficient timing corroboration against {source.name}: {matched}")
    return {"archived_sws_events": int(len(ref_peaks)), "within_20ms": matched,
            "median_nearest_error_s": float(np.median(errors)),
            "lfp_event_center_minus_baseline_native": float(np.median(contrast("lfsb"))),
            "lfp_random_center_minus_baseline_native": float(np.median(contrast("lfrnd")))}


def prepare_record(item, root, params, config, assets, supplemental):
    metrics, ticks, groups, _, rate = lman_and_states(root, item, params, config, assets)
    bird, session = item["bird"], item["session"]
    record_dir = f"dataset/sub-{bird}/ses-{session}/derived/"
    rows = read_table(fetch(root, record_dir + item["score"], assets))
    complete = min(int(metrics["ap_duration_s"]), int(metrics["last_spike_s"]))
    rows = rows[:complete]
    all_units = np.unique(groups)
    seconds = ticks / rate
    trains = [seconds[groups == unit] for unit in all_units]
    catalog = []
    rules = config["reanalysis"]
    low, high = metrics["lman_channel_range"]
    for row in read_table(fetch(root, record_dir + item["prefix"] + "cluster_info.tsv", assets), "\t"):
        if not low <= int(float(row["ch"])) <= high:
            continue
        unit = int(row.get("cluster_id", row.get("id")))
        ind = np.flatnonzero(all_units == unit)
        contamination = float(row["ContamPct"]) if row.get("ContamPct") else None
        isi = np.diff(trains[ind[0]]) if len(ind) else np.array([])
        refractory = float(np.mean(isi < rules["refractory_seconds"])) if len(isi) else None
        reasons = []
        if row["group"].strip().lower() != "good":
            reasons.append("not_manually_labeled_good")
        if not len(ind):
            reasons.append("no_eligible_spikes_or_nonpositive_rate")
        if contamination is None or contamination > rules["maximum_contamination_percent"]:
            reasons.append("missing_or_high_contamination")
        if refractory is None or refractory > rules["maximum_refractory_fraction"]:
            reasons.append("missing_or_high_refractory_violation")
        catalog.append({"id": unit, "label": row["group"] or "unlabeled", "channel": int(row["ch"]),
                        "published_rate_hz": float(row["fr"]), "contamination_percent": contamination,
                        "amplitude": float(row["Amplitude"]) if row.get("Amplitude") else None,
                        "refractory_fraction": refractory, "broad": bool(len(ind)), "strict": not reasons,
                        "strict_exclusion_reasons": reasons})
    strict_ids = [row["id"] for row in catalog if row["strict"]]
    strict_mask = np.isin(all_units, strict_ids)

    width = complete * BINS_PER_SECOND
    bins = np.floor(seconds * BINS_PER_SECOND).astype(np.int64)
    valid = bins < width
    indices = np.searchsorted(all_units, groups[valid]) * width + bins[valid]
    counts = np.bincount(indices, minlength=len(trains) * width).reshape(len(trains), width)
    assert counts.max() <= np.iinfo(np.uint16).max
    counts = counts.astype(np.uint16)
    events, time, participation, threshold, checks = detect(trains, metrics["ap_duration_s"])
    context = {"recording_checks": metrics, "trimmed_at_s": complete, "unit_qc": catalog,
               "unit_label_counts": dict(Counter(row["label"] for row in catalog)),
               "burst_threshold": threshold, "checked_centers_s": checks}
    if item["state"] == "WAKE":
        artifact = read_table(fetch(root, f"{bird}_adResults.csv", supplemental))
        assert [int(r["sec"]) for r in artifact] == list(range(len(artifact)))
        assert len(artifact) - complete in range(0, 15)
        mask = np.array([int(r["artefactDetect"]) == 1 for r in artifact[:complete]])
        shifted = np.r_[False, mask[:-1]]
        score = np.array([int(r["soundDetected"]) for r in rows])
        assert set(np.unique(score)) <= {0, 1}
        masks = {"ACTIVE": score == 1, "NON_ACTIVE": score == 0}
        context["artifact_seconds_released"] = int(mask.sum())
        context["artifact_seconds_union"] = int((mask | shifted).sum())
    else:
        scores = np.asarray([int(float(r["final_states_new"])) for r in rows])
        mic = np.asarray([int(float(r["new_mic"])) for r in rows])
        masks = {"SWS": (scores == config["sws_code"]) & (mic == 0)}
        mask = shifted = np.zeros(complete, dtype=bool)
        context["archived_event_validation"] = trace_summary(
            events, fetch(root, f"LF_{bird}.mat", supplemental))

    state_data, state_info, event_tracks = {}, {}, {}
    for state, original in masks.items():
        bounds = continuous(original & ~(mask | shifted))
        exposure = int(sum(end - start for start, end in bounds))
        if exposure < config["minimum_state_seconds"]:
            state_info[state] = {"status": "BLOCKED", "seconds": exposure}
            continue
        segments = segment_matrix(counts, bounds)
        strict_segments = segment_matrix(counts[strict_mask], bounds)
        strict_active = int(np.count_nonzero(sum(x.sum(axis=1) for x in strict_segments[:, 0])))
        w = rules["window_seconds"]
        windows = np.array([(t, t + w) for start, end in bounds
                            for t in range(start, end - w + 1, w)], dtype=np.int64).reshape(-1, 2)
        window_features = [{"start_s": int(t), **summarize(events, time, participation, [(t, t + w)])}
                           for t, _ in windows]
        lag_guard = min(x.shape[1] for x in segments[:, 0])
        assert lag_guard > 10
        state_data[state] = {
            "segments": segments, "broad_unit_ids": all_units.astype(np.int32),
            "strict": strict_segments if strict_active >= rules["minimum_good_units"] else np.empty((0, 1), dtype=object),
            "unmasked": segment_matrix(counts, continuous(original)),
            "literal_mask": segment_matrix(counts, continuous(original & ~mask)),
            "windows": segment_matrix(counts, windows), "window_starts": windows[:, 0].astype(float)}
        state_info[state] = {"status": "PASS_BROAD", "strict_status": "ELIGIBLE_UNPAIRED" if len(state_data[state]["strict"]) else "BLOCKED",
                             "seconds": exposure, "intervals": bounds.tolist(),
                             "strict_active": strict_active, "strict_ids": strict_ids,
                             "units": int(len(all_units)), "unmasked_seconds": int(original.sum()),
                             "literal_mask_seconds": int((original & ~mask).sum()),
                             "window_seconds": w, "window_features": window_features,
                             "burst_summary": summarize(events, time, participation, bounds)}
        event_tracks[state] = [e for e in events if any(
            start <= e["onset_s"] and e["offset_s"] <= stop for start, stop in bounds)]
    return state_data, state_info, context, event_tracks


def atomic_text(path, content):
    stage = path.with_suffix(".tmp")
    stage.write_text(json.dumps(content, indent=2, allow_nan=False) + "\n")
    stage.replace(path)


def prepare():
    paths = json.loads(Path("config/paths.json").read_text())
    config = json.loads(Path("config/zebra_finch.json").read_text())
    root = Path(paths["drive_root"])
    transfer = json.loads((root / "manifests/zebra-finch-transfer.json").read_text())
    extra = json.loads(Path("manifests/zebra-reanalysis-sources.json").read_text())
    assert transfer["version"] == extra["commit"]
    assets = transfer["files"]
    params = loadmat(fetch(root, "code/UtilsMatlab/birdParams.mat", assets), simplify_cells=True)["params"]
    fixture, evidence, recording_qc, tracks = {}, {}, {}, {}
    for item in config["recordings"]:
        bird, session = item["bird"], item["session"]
        states, stats, context, events = prepare_record(item, root, params, config, assets, extra["files"])
        fixture.setdefault(bird, {}).update(states)
        evidence.setdefault(bird, {}).update(stats)
        tracks.setdefault(bird, {}).update(events)
        recording_qc[f"{bird}_{session}"] = context
        print(f"{bird} {session}: " + str({s: (v.get("seconds"), v.get("strict_active"),
              len(v.get("window_features", [])), v.get("burst_summary", {}).get("events")) for s, v in stats.items()}), flush=True)

    repetitions = config["reanalysis"]["repetitions"]
    seed = config["seed"]
    for bird, states in fixture.items():
        assert set(states) == set(STATES)
        ids_by_record = {}
        minimum = min(len(row["broad_unit_ids"]) for row in states.values())
        for state, row in states.items():
            key = "night" if state == "SWS" else "light"
            ids = row["broad_unit_ids"]
            if key not in ids_by_record:
                rng = np.random.default_rng(np.random.SeedSequence([seed, int(bird[-1]), int(key == "night")]))
                picks = np.stack([rng.choice(len(ids), minimum, replace=False) for _ in range(repetitions)])
                schedules = np.zeros((repetitions, len(ids)), dtype=bool)
                schedules[np.arange(repetitions)[:, None], picks] = True
                ids_by_record[key] = ids, schedules
            previous, schedules = ids_by_record[key]
            assert np.array_equal(ids, previous)
            row["unit_masks"] = schedules.astype(np.uint8)
    target = root / "matlab_inputs/zebra_lman_reanalysis_v2.mat"
    with TemporaryDirectory(prefix="zebra-v2-", dir="/content") as scratch:
        path = Path(scratch) / target.name
        savemat(path, {"cohort": fixture}, do_compression=True)
        size = path.stat().st_size
        if size >= 25 * 2**20:
            raise ValueError(f"Derived fixture exceeds validated Git transfer budget: {size}")
        with path.open("rb") as handle:
            sha = file_digest(handle, "sha256").hexdigest()
        stage = target.with_suffix(".part")
        stage.write_bytes(path.read_bytes())
        with stage.open("rb") as handle:
            assert file_digest(handle, "sha256").hexdigest() == sha
        stage.replace(target)

    for bird in fixture:
        maps = [row["recording_checks"]["probe_map"] for name, row in recording_qc.items()
                if name.startswith(bird + "_")]
        assert len(maps) == 2 and maps[0] and maps[0] == maps[1], f"Probe map mismatch: {bird}"
    for row in recording_qc.values():
        row["recording_checks"].pop("probe_map")
    result = {"schema_version": 2, "source_commit": transfer["version"],
              "created_at_utc": datetime.now(timezone.utc).isoformat(),
              "settings": config, "fixture_bytes": size, "fixture_sha256": sha,
              "states": evidence, "recording_qc": recording_qc,
              "strict_population_across_states": {"status": "BLOCKED", "reason": "insufficient good LMAN units"}}
    atomic_text(root / "manifests/zebra-finch-preparation-v2.json", result)
    atomic_text(root / "processed/zebra_finch/lman-burst-events-v2.json", tracks)
    print(f"Prepared {size} bytes, SHA-256 {sha}")


if __name__ == "__main__":
    prepare()
