"""Combine targeted L-ratio subsets with an exact-input shorter-window stability check."""
import argparse
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np
from scipy.io import loadmat

from build_atlas_data import ROOT, read
sys.path.insert(0, str(ROOT / "scripts/data"))
from prepare_sleep_dandi import save_prepared


def count_entries(data):
    for entry in loadmat(BytesIO(data), simplify_cells=True)["cohort"].values():
        m = json.loads(entry["metadata"])
        segments = entry["counts"]
        if segments.dtype != object:
            segments = [segments]
        yield m, list(segments)


def cells(segments):
    block = np.empty((len(segments), 1), dtype=object)
    block[:, 0] = segments
    return block


def paired_window_tasks(archive, cases):
    receipt_file = "cohort/phase10_input.json"
    with ZipFile(archive) as z:
        receipt = json.loads(z.read(receipt_file))
        raw = z.read("cohort/" + receipt["fixture"])
        assert sha256(raw).hexdigest() == receipt["sha256"]
    fit = read("results/phase10/fits/posub-adn-supplement.json")
    source = next(x for x in read("manifests/phase10_experiments.json")["runs"] if x["batch"] == "posub-adn-supplement")
    assert fit["project_commit"] == source["project_commit"]
    source_rows = {(row["record"]["asset_id"], row["record"]["state"]): row for row in fit["records"]
                   if row["record"]["dandiset"] == "000939" and row["record"]["variant"] == "reference"}
    selected = []
    assets = {x["asset_id"]: x["case"] for x in cases}
    for m, segments in count_entries(raw):
        if m["asset_id"] not in assets or m["variant"] != "reference":
            continue
        assert source_rows[(m["asset_id"], m["state"])]["record"] == m
        fit_row = source_rows[(m["asset_id"], m["state"])]["fit"]
        lengths = [len(s) for s in segments]
        assert lengths == fit_row["segment_bins"]
        assert sum(int(seg.sum()) for seg in segments) == round(fit_row["rate_hz"] * fit_row["duration_s"])
        for seconds, n in ((60, 3), (45, 4)):
            if assets[m["asset_id"]] == "E" and seconds == 60:
                continue
            bins = round(seconds / .04)
            windows = [
                (seg[i:i+bins], [round((a + .04 * i), 8), round(a + .04 * (i+bins), 8)])
                for (a, _), seg in zip(m["intervals_s"], segments, strict=True)
                for i in range(0, len(seg)-bins+1, bins)
            ]
            assert len(windows) >= n
            chosen = [windows[i] for i in np.linspace(0,len(windows)-1,n,dtype=int)]
            spans, intervals = zip(*chosen, strict=True)
            assert sum(map(len, spans)) == 4500
            selected.append(({**m, "variant": f"window_{seconds}s_180s",
                              "intervals_s": list(intervals), "audit_case": assets[m["asset_id"]]}, cells(spans)))
    assert len(selected) == (6 if "C" in assets.values() else 3)
    return selected, receipt["sha256"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("original_000939_archive", type=Path)
    p.add_argument("--sorting", type=Path)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--case", choices=("C", "E"), default="C")
    a = p.parse_args()
    cases = [x for x in read("results/validation/selected_sessions.json")["cases"] if x["case"] == a.case]
    tasks, sha = paired_window_tasks(a.original_000939_archive, cases)
    inputs = {"000939_cohort_counts": sha}
    if a.sorting:
        p = a.sorting.with_suffix(".json")
        original = json.loads(p.read_text())
        raw = a.sorting.read_bytes()
        assert sha256(raw).hexdigest() == original["sha256"]
        sort_rows = [(m, cells(rows)) for m, rows in count_entries(raw) if m["variant"].startswith("lratio_")]
        assert len(sort_rows) == 18
        tasks.extend(sort_rows)
        inputs.update(sorting_source_counts=original["sha256"], **original["source_sha256"])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    info = save_prepared(tasks, a.out, "states", sorted({m["asset_id"] for m, _ in tasks}), source_sha256=inputs)
    relative = a.out.resolve().relative_to(ROOT)
    a.out.with_suffix(".sha256").write_text(info["sha256"] + "  " + str(relative) + "\n", encoding="utf-8")
    print(info["records"], "fits", info["sha256"])


if __name__ == "__main__":
    main()
