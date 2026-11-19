"""Export small LMAN population-count fixtures once the paired audit passes."""

from datetime import datetime, timezone
from hashlib import file_digest
from pathlib import Path
import json
import shutil
from tempfile import TemporaryDirectory

import numpy as np
from scipy.io import loadmat, savemat

from inspect_zebra import fetch, lman_and_states


def counts_by_interval(times, bounds, rate, bin_seconds):
    counts = np.empty((len(bounds), 1), dtype=object)
    bins_per_second = round(1 / bin_seconds)
    for i, (start, stop) in enumerate(bounds):
        # Author scores are one-second rows. Half-open bins keep only
        # full bins and never turn a segment boundary into a singleton bin.
        edges = (start + np.arange((stop - start) * bins_per_second + 1)
                 / bins_per_second) * rate
        indices = np.searchsorted(times, np.ceil(edges).astype(np.int64), side="left")
        sequence = np.diff(indices)
        if sequence.max() > np.iinfo(np.uint16).max:
            raise ValueError("Population count exceeds uint16 fixture capacity")
        counts[i, 0] = sequence.astype(np.uint16)
    return counts


def prepare():
    paths = json.loads(Path("config/paths.json").read_text())
    root = Path(paths["drive_root"])
    audit = json.loads((root / "manifests/zebra-finch-audit.json").read_text())
    source = json.loads(Path(paths["source_manifest"]).read_text())["zebra_finch"]
    assert audit["source_commit"] == source["commit"]
    settings = audit["parameters"]
    assets = audit["files"]
    params_file = fetch(root, "code/UtilsMatlab/birdParams.mat", assets)
    params = loadmat(params_file, simplify_cells=True)["params"]
    pairs = sorted({row["bird"] for row in audit["recordings"]
                    if row["eligible"] and row["state"] == "WAKE"}
                   & {row["bird"] for row in audit["recordings"]
                      if row["eligible"] and row["state"] == "SWS"})
    assert len(pairs) >= 2, "At least two eligible birds are required for paired exploration."
    rng = np.random.default_rng(settings["seed"])
    fixture = {}
    summaries = []
    for bird in pairs:
        per_state = {}
        for item in settings["recordings"]:
            if item["bird"] == bird:
                per_state[item["state"]] = lman_and_states(root, item, params, settings, assets)
        assert set(per_state) == {"WAKE", "SWS"}
        assert per_state["WAKE"][0]["probe_map"] == per_state["SWS"][0]["probe_map"], +            f"AP probe mapping mismatch in {bird}"
        unit_count = min(out[0]["lman_units"] for out in per_state.values())
        fixture[bird] = {}
        for state, (record, times, labels, bounds, rate) in per_state.items():
            primary = counts_by_interval(times, bounds, rate, settings["primary_bin_seconds"])
            narrow = counts_by_interval(times, bounds, rate, 0.02)
            units = rng.choice(np.unique(labels), size=unit_count, replace=False)
            matched = counts_by_interval(times[np.isin(labels, units)], bounds, rate, 0.04)
            fixture[bird][state] = {"primary": primary, "bin20": narrow, "matched_units": matched}
            pooled = np.concatenate(primary[:, 0])
            assert pooled.var() > 0
            record = {k: v for k, v in record.items() if k != "probe_map"}
            summaries.append({**record, "matched_unit_count": unit_count,
                              "matched_unit_ids": units.tolist(), "population_spikes_in_state": int(pooled.sum()),
                              "mean_count_per_bin": float(pooled.mean()), "count_variance": float(pooled.var()),
                              "silent_bin_fraction": float(np.mean(pooled == 0)),
                              "population_rate_hz": float(pooled.sum() / record["eligible_seconds"]),
                              "rate_per_unit_hz": float(pooled.sum() / record["eligible_seconds"] / record["lman_units"])})
        mean_count = {state: float(np.concatenate(data["primary"][:, 0]).mean())
                      for state, data in fixture[bird].items()}
        target_mean = min(mean_count.values())
        for state, data in fixture[bird].items():
            counts = data["primary"]
            downsampled = np.empty_like(counts)
            for index, sequence in enumerate(counts[:, 0]):
                downsampled[index, 0] = rng.binomial(sequence, target_mean / mean_count[state]).astype(np.uint16)
            data["matched_rate"] = downsampled
        print(f"{bird} prepared with {unit_count} matched units", flush=True)

    destination = root / "matlab_inputs/zebra_lman_population.mat"
    with TemporaryDirectory(prefix="finch-fixture-", dir="/content") as scratch_dir:
        scratch = Path(scratch_dir) / destination.name
        savemat(scratch, {"cohort": fixture}, do_compression=True)
        size = scratch.stat().st_size
        assert size < 25 * 2**20
        with scratch.open("rb") as source_file:
            digest = file_digest(source_file, "sha256").hexdigest()
        stage = destination.with_suffix(".mat.part")
        shutil.copyfile(scratch, stage)
        with stage.open("rb") as saved:
            assert file_digest(saved, "sha256").hexdigest() == digest
        stage.replace(destination)
    info = {"source_commit": audit["source_commit"], "parameters": settings,
            "created_at_utc": datetime.now(timezone.utc).isoformat(), "pairs": pairs,
            "fixture_bytes": size, "fixture_sha256": digest, "sessions": summaries}
    target = root / "manifests/zebra-finch-preparation.json"
    stage = target.with_suffix(".tmp")
    stage.write_text(json.dumps(info, indent=2) + "\n")
    stage.replace(target)
    print(f"Prepared {len(pairs)} birds, {size} bytes, SHA-256 {digest}")


if __name__ == "__main__":
    prepare()
