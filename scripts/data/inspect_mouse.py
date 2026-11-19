"""Stage a published mouse session in Drive and report source-level QC."""

from datetime import datetime, timezone
from hashlib import file_digest
from pathlib import Path
import json
import shutil
import sys
from tempfile import TemporaryDirectory
from urllib.request import urlopen

import numpy as np
from scipy.io import loadmat, whosmat


STATES = ("WAKE", "NREM", "REM")


def verify_file(path, record):
    with path.open("rb") as source:
        md5 = file_digest(source, "md5").hexdigest()
    if path.stat().st_size != record["bytes"] or md5 != record["md5"]:
        raise ValueError(f"Size or Figshare MD5 mismatch for {record['file']}")


def source_file(root, record):
    destination = root / record["drive_path"]
    if destination.exists():
        verify_file(destination, record)
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="figshare-", dir="/content") as directory:
        local = Path(directory) / record["file"]
        with urlopen(record["url"], timeout=60) as response, local.open("wb") as target:
            shutil.copyfileobj(response, target)
        verify_file(local, record)
        staged = destination.with_suffix(".mat.part")
        shutil.copyfile(local, staged)
        verify_file(staged, record)
        staged.replace(destination)
    return destination


def inspect(session, config, mouse):
    record = next(item for item in mouse["files"] if item["session"] == session)
    root = Path(config["drive_root"])
    if not (root / "manifests").is_dir():
        raise RuntimeError("Mount the verified project Drive with 00_drive_setup.ipynb first.")
    path = source_file(root, record)
    fields = whosmat(path)
    data = loadmat(path, variable_names=["st", "clu", *STATES])
    spikes = np.asarray(data["st"]).reshape(-1)
    labels = np.asarray(data["clu"]).reshape(-1)
    if not len(spikes) or len(labels) != len(spikes) or not np.all(np.isfinite(spikes)):
        raise ValueError(f"Invalid spike times or cluster labels in {session}")
    states = {}
    for state in STATES:
        epochs = np.asarray(data[state])
        if epochs.ndim != 2 or epochs.shape[1] != 2 or not np.all(np.isfinite(epochs)):
            raise ValueError(f"Invalid {state} intervals in {session}")
        lengths = epochs[:, 1] - epochs[:, 0]
        if np.any(lengths <= 0):
            raise ValueError(f"Invalid {state} durations in {session}")
        states[state] = {"intervals": len(epochs), "seconds": float(lengths.sum()),
                         "min_seconds": float(lengths.min()), "max_seconds": float(lengths.max())}

    result = {**record, "status": "downloaded_and_verified", "source_version": mouse["doi"],
              "verified_at_utc": datetime.now(timezone.utc).isoformat(),
              "scratch_path": f"/content/figshare-*/{record['file']}",
              "mat_fields": fields, "spikes": len(spikes), "units": len(np.unique(labels)),
              "spike_time_unit": "seconds", "spike_time_range_s": [float(spikes.min()), float(spikes.max())],
              "spike_times_nondecreasing": bool(np.all(spikes[1:] >= spikes[:-1])),
              "states": states}
    output = root / "manifests" / f"mouse-{session}-inventory.json"
    staged = output.with_suffix(".tmp")
    staged.write_text(json.dumps(result, indent=2) + "\n")
    staged.replace(output)
    return result


if __name__ == "__main__":
    config = json.loads(Path("config/paths.json").read_text())
    source = json.loads(Path(config["source_manifest"]).read_text())
    print(json.dumps(inspect(sys.argv[1], config, source["mouse"]), indent=2))
