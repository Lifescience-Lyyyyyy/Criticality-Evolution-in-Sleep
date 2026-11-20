"""Audit the exact CRCNS hc-24 release and prepare only author-scored NREM."""
import argparse
import json
from hashlib import sha256 as digest_bytes
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.io import loadmat, whosmat

from crcns_adapter import _verified_record, clean_clock, clip, intervals, sha256
from prepare_crcns_sleep import build
from prepare_sleep_dandi import clock_bins

DATA_SHA256 = "187e25461f9c2cb118eec825957ed2938e34809ce4211442d9229b249061d738"
SESSIONS = {
    "Session1_PETHphaseLocking": ("Rat-1", "REST"),
    "Session2_ThalPairCorrelations": ("Rat-2", "REST"),
    "Session3_SleepAfterRun": ("Rat-3", "POST_RUN"),
}
REGIONS = {"CA1": "CA1", "PFC": "mPFC", "THAL": "Midline thalamus"}


def audit_session(package, session, qualified_names, cfg):
    subject, context = SESSIONS[session]
    prefix = f"data/{session}/"
    names = sorted(name for name in qualified_names if name.startswith(prefix))
    timestamp_path, = (name for name in names if name.endswith("/Timestamp_LFP_600Hz.mat"))
    lfp = loadmat(BytesIO(package.read(timestamp_path)), variable_names=["timestamp"])["timestamp"].ravel()
    if not clean_clock(lfp) or not np.allclose(np.diff(lfp), 1 / 600, rtol=0, atol=1e-8):
        raise ValueError(f"{session}: irregular native 600 Hz clock")
    lfp_signals = [name for name in names if name.endswith(("_LFP_600Hz.mat", "/Timestamp_LFP_600Hz.mat"))
                   and name != timestamp_path]
    for name in lfp_signals:
        _, shape, _ = whosmat(BytesIO(package.read(name)))[0]
        if int(np.prod(shape)) != len(lfp):
            raise ValueError(f"{name}: signal length conflicts with timestamp vector")
    bounds = [(float(lfp[0]), float(lfp[-1] + 1 / 600))]
    score_path = prefix + "DerivedData/SleepTimes.mat"
    score_bytes = package.read(score_path)
    scores = loadmat(BytesIO(score_bytes))
    if {key for key in scores if not key.startswith("__")} != {"nremSleep"}:
        raise ValueError(f"{session}: unexpected author state variables")
    nrem = intervals(scores["nremSleep"])
    if clip(nrem, bounds) != nrem:
        raise ValueError(f"{session}: author score outside recorded support")
    width = cfg["bin_seconds"]
    support = clock_bins(nrem, width, cfg["min_segment_seconds"])
    seconds = round(width * sum(stop - start for start, stop in support), 8)
    trains = {region: {} for region in REGIONS.values()}
    units = {region: [] for region in REGIONS.values()}
    for name in names:
        path = Path(name)
        if not ("/PrimaryData/" in name and path.suffix == ".mat" and "cl-" in path.name):
            continue
        # Session 2's author demo identifies its two top-level clusters as thalamic.
        key = path.parent.name.split("_", 1)[0].upper()
        region = "Midline thalamus" if session == "Session2_ThalPairCorrelations" else REGIONS[key]
        data = loadmat(BytesIO(package.read(name)), variable_names=["unitdata"])["unitdata"]
        if data.ndim != 2 or not data.shape[0] or data.shape[1] != 5 or not np.isfinite(data[:, 4]).all():
            raise ValueError(f"{name}: invalid five-column cluster matrix")
        t = data[:, 4].astype(float)
        in_support = t[(t >= bounds[0][0]) & (t < bounds[0][1])]
        diffs = np.diff(t)
        nrem_spikes = sum(int(np.diff(np.searchsorted(in_support, (start, stop)))[0])
                          for start, stop in nrem) if clean_clock(t) else None
        uid = str(path.relative_to(prefix + "PrimaryData"))
        if nrem_spikes:
            trains[region][uid] = in_support
        units[region].append({"id": uid, "spikes": len(t), "outside_lfp_spikes": len(t) - len(in_support),
                              "ties": int(np.count_nonzero(diffs == 0)),
                              "reversals": int(np.count_nonzero(diffs < 0)),
                              "nrem_spikes": nrem_spikes, "accepted": bool(nrem_spikes)})
    eligible = [region for region, spikes in trains.items()
                if len(spikes) >= cfg["min_units"] and seconds >= cfg["min_state_seconds"]]
    source = {"asset_id": f"hc-24:{session}", "sha256": DATA_SHA256,
              "path": prefix.rstrip("/"), "state_source": score_path + ":nremSleep"}
    records = []
    for region in eligible:
        records.append(_verified_record(
            {"subject_id": subject, "session_id": session, "source_dataset": "hc-24", "region": region,
             "context": context, "timebase": {"units": "seconds", "lfp_hz": 600},
             "provenance": source, "sleep_intervals": {"NREM": nrem},
             "recording_intervals": bounds,
             "state_evidence": {"method": "author_scored", "reference": "https://doi.org/10.7554/eLife.48881"}},
            trains[region],
        ))
    result = {"subject": subject, "session": session, "context": context,
              "clock_source": timestamp_path, "clock_samples": len(lfp),
              "recording_seconds": round(float(lfp[-1] - lfp[0]), 8),
              "lfp_signals": lfp_signals, "state_source": source["state_source"],
              "state_sha256": digest_bytes(score_bytes).hexdigest(),
              "nrem_bouts": len(nrem), "nrem_seconds": seconds,
              "regions": [{"region": region, "released_units": len(units[region]),
                           "valid_units": len(trains[region]), "eligible_nrem": region in eligible,
                           "unit_qc": units[region]} for region in units]}
    return result, records


def audit(archive, cfg):
    if sha256(archive) != DATA_SHA256:
        raise ValueError("Archive does not match the verified CRCNS hc-24 release")
    with ZipFile(archive) as package:
        names = package.namelist()
        if len(set(names)) != len(names) or package.testzip() is not None:
            raise ValueError("CRCNS archive has a repeated or corrupt member")
        results = [audit_session(package, name, names, cfg) for name in SESSIONS]
    return {"archive_sha256": DATA_SHA256, "archive_bytes": archive.stat().st_size,
            "sessions": [audit_row for audit_row, _ in results]}, [
                record for _, records in results for record in records]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--config", type=Path, default=Path("config/sleep_dandi.json"))
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--session", choices=tuple(SESSIONS))
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    receipt, recordings = audit(args.archive, config)
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if args.out:
        selected = [r for r in recordings if args.session is None or args.session == r.session_id]
        if not selected:
            parser.error("No eligible NREM recordings for the selected session")
        print(build(selected, args.out, config, "nrem", {"archive": DATA_SHA256}))
