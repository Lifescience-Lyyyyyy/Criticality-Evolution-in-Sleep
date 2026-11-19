"""Build a compact state-segmented count fixture from the 13 PNAS sessions."""

from datetime import datetime, timezone
from hashlib import file_digest
from pathlib import Path
import json
import shutil
from tempfile import TemporaryDirectory

import numpy as np
from scipy.io import loadmat, savemat

from inspect_mouse import source_file


def sample_ticks(seconds, sample_rate):
    scaled = np.asarray(seconds) * sample_rate
    ticks = np.rint(scaled)
    if not np.all(np.isfinite(scaled)) or np.max(np.abs(scaled - ticks)) > 1e-4:
        raise ValueError("Timestamps do not align with the published acquisition clock")
    return ticks.astype(np.int64)


def bin_epoch(spike_ticks, start, stop, width):
    edges = np.arange(start, stop + 1, width, dtype=np.int64)
    indices = np.searchsorted(spike_ticks, edges, side="left")
    last_edge_count = int(np.searchsorted(spike_ticks, edges[-1], side="right") - indices[-1])
    return np.append(np.diff(indices), last_edge_count)


def prepare_session(root, record, config):
    raw = source_file(root, record)
    states = config["states"]
    data = loadmat(raw, variable_names=["Session", "st", "clu", *states])
    rate = config["sample_rate_hz"]
    width = round(config["bin_seconds"] * rate)
    spikes = sample_ticks(data["st"].reshape(-1), rate)
    labels = data["clu"].reshape(-1)
    if len(spikes) != len(labels) or np.any(np.diff(spikes) < 0):
        raise ValueError(f"Cluster alignment or spike order failed in {record['session']}")

    times = {state: np.asarray(data[state]) for state in states}
    combined = np.concatenate(list(times.values()), axis=0)
    if np.any(combined[:, 1] < combined[:, 0]):
        raise ValueError(f"Reversed state intervals in {record['session']}")
    combined = combined[combined[:, 1] > combined[:, 0]]
    combined = combined[np.argsort(combined[:, 0])]
    if np.any(combined[1:, 0] < combined[:-1, 1]):
        raise ValueError(f"Overlapping state intervals in {record['session']}")

    segments = {}
    diagnostics = {}
    for state, intervals in times.items():
        ticks = sample_ticks(intervals, rate)
        duration = float(np.diff(intervals, axis=1).sum())
        if duration < config["minimum_state_seconds"]:
            raise ValueError(f"{record['session']} {state} has only {duration} scored seconds")
        zero = ticks[ticks[:, 0] == ticks[:, 1], 0]
        if zero.size:
            # histc would return a singleton count for an interval of zero length.
            segments[f"{state}_zero_duration_counts"] = [
                int(bin_epoch(spikes, boundary, boundary, width)[0]) for boundary in zero
            ]
        ticks = ticks[ticks[:, 1] > ticks[:, 0]]
        counts = np.empty((len(ticks), 1), dtype=object)
        for index, (start, stop) in enumerate(ticks):
            sequence = bin_epoch(spikes, start, stop, width)
            if sequence.max() > np.iinfo(np.uint16).max:
                raise ValueError(f"{record['session']} {state} exceeds uint16 bin capacity")
            counts[index, 0] = sequence.astype(np.uint16)
        pooled = np.concatenate(counts[:, 0])
        if pooled.std() == 0:
            raise ValueError(f"{record['session']} {state} has no population count variability")
        segments[state] = counts
        diagnostics[state] = {"seconds": duration, "segments": len(counts), "bins": len(pooled),
                              "spikes": int(pooled.sum()), "mean_count": float(pooled.mean()),
                              "variance": float(pooled.var()),
                              "zero_duration_start_s": (zero / rate).tolist()}

    summary = {**record, "status": "prepared", "spikes": len(spikes),
               "source_session_label": str(data["Session"].item()),
               "units": len(np.unique(labels)), "first_spike_s": float(spikes[0] / rate),
               "last_spike_s": float(spikes[-1] / rate), "states": diagnostics}
    return segments, summary


def main():
    paths = json.loads(Path("config/paths.json").read_text())
    config = json.loads(Path("config/pnas.json").read_text())
    manifest = json.loads(Path(paths["source_manifest"]).read_text())
    root = Path(paths["drive_root"])
    if not (root / "manifests").is_dir():
        raise RuntimeError("Mount the verified project Drive before preparing mouse data.")
    destination = root / "matlab_inputs" / "mouse_population_40ms.mat"
    inventory_path = root / "manifests" / "mouse-preparation.json"
    with Path(__file__).open("rb") as code:
        preparation_sha = file_digest(code, "sha256").hexdigest()
    if destination.exists() and inventory_path.exists():
        previous = json.loads(inventory_path.read_text())
        with destination.open("rb") as fixture:
            verified = file_digest(fixture, "sha256").hexdigest() == previous["fixture_sha256"]
        if (verified and previous["parameters"] == config
                and previous["source_version"] == manifest["mouse"]["doi"]
                and previous["preparation_code_sha256"] == preparation_sha):
            print(f"Reused verified fixture {destination.name}")
            return
    np.testing.assert_array_equal(bin_epoch(np.array([0, 799, 800, 1599, 1600]), 0, 1600, 800),
                                  [2, 2, 1])
    cohort = {}
    summaries = []
    for record in manifest["mouse"]["files"]:
        cohort[record["session"]], summary = prepare_session(root, record, config)
        summaries.append(summary)
        print(f"{record['session']} ready: {summary['units']} units", flush=True)

    with TemporaryDirectory(prefix="mouse-fixture-", dir="/content") as directory:
        scratch = Path(directory) / destination.name
        savemat(scratch, {"cohort": cohort}, do_compression=True)
        if scratch.stat().st_size >= 25 * 2**20:
            raise ValueError("The prepared fixture is too large for Git transfer")
        with scratch.open("rb") as source:
            digest = file_digest(source, "sha256").hexdigest()
        stage = destination.with_suffix(".mat.part")
        shutil.copyfile(scratch, stage)
        with stage.open("rb") as saved:
            assert file_digest(saved, "sha256").hexdigest() == digest
        stage.replace(destination)

    inventory = {"source_version": manifest["mouse"]["doi"], "parameters": config,
                 "preparation_code_sha256": preparation_sha,
                 "created_at_utc": datetime.now(timezone.utc).isoformat(),
                 "fixture_bytes": destination.stat().st_size, "fixture_sha256": digest,
                 "sessions": summaries}
    stage = inventory_path.with_suffix(".tmp")
    stage.write_text(json.dumps(inventory, indent=2) + "\n")
    stage.replace(inventory_path)
    print(f"Prepared {len(cohort)} sessions, {inventory['fixture_bytes']} bytes, SHA-256 {digest}")


if __name__ == "__main__":
    main()
