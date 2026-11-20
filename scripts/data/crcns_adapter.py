"""Adapt audited sources to unit timestamps and independent state intervals."""

from collections import Counter, defaultdict
from dataclasses import dataclass
from hashlib import file_digest
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile
import re

import numpy as np
from scipy.io import loadmat

from audit_crcns_sleep import HC23_MEMBER


STATES = ("Wake", "NREM", "REM")
HC23_MAT_SHA256 = "e6c335643b74ad60bc65c984d3046d52711ecc82ca10a2086f64c6ba581b87f7"


@dataclass(frozen=True)
class SleepRecording:
    subject_id: str
    session_id: str
    source_dataset: str
    region: str
    unit_ids: tuple
    spike_times: tuple
    sleep_intervals: dict
    recording_intervals: tuple
    timebase: dict
    provenance: dict
    quality_metrics: dict
    context: str

    @property
    def asset_id(self):
        return self.provenance["asset_id"]


def sha256(path):
    with Path(path).open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def mat_items(value):
    if isinstance(value, dict):
        return [value]
    return list(value) if np.size(value) else []


def intervals(value):
    array = np.asarray(value, dtype=float).reshape(-1, 2)
    if not np.isfinite(array).all() or np.any(array[:, 0] >= array[:, 1]):
        raise ValueError("Intervals must have finite, increasing endpoints")
    return [tuple(e) for e in array[np.argsort(array[:, 0], kind="stable")]]


def overlaps(ranges):
    return any(b > c + 1e-8 for (_, b), (c, _) in zip(ranges, ranges[1:]))


def merge_same_state(ranges):
    merged = []
    for left, right in ranges:
        if merged and left < merged[-1][1] - 1e-8:
            start, stop = merged[-1]
            merged[-1] = start, max(stop, right)
        else:
            merged.append((left, right))
    return merged


def clip(ranges, bounds):
    return [(max(a, x), min(b, y)) for a, b in ranges for x, y in bounds if max(a, x) < min(b, y)]


def clean_clock(spikes):
    return np.isfinite(spikes).all() and np.all(np.diff(spikes) > 0)


def parse_tfile(data):
    end = data.find(b"%%ENDHEADER")
    if end < 0 or not data.startswith(b"%%BEGINHEADER"):
        raise ValueError("Unexpected Neuralynx t-file header")
    offset = data.index(b"\n", end) + 1
    return np.frombuffer(data[offset:], dtype=">u4").astype(np.int64)


def _raw_bouts(path):
    groups = defaultdict(dict)
    with ZipFile(path) as z:
        for name in z.namelist():
            match = HC23_MEMBER.fullmatch(name)
            if not match:
                raise ValueError(f"Unexpected member: {name}")
            rat, session, bout, tetrode, cluster, nongrid = match.groups()
            key = rat, session, int(bout[5:])
            unit = (int(tetrode), int(cluster), 2 if nongrid else 1)
            if unit in groups[key]:
                raise ValueError(f"Duplicate sorted cluster and type in {name}")
            data = z.read(name)
            groups[key][unit] = None if data and not any(data) else parse_tfile(data)
    return groups


def hc23_prepared(mat, extracted, source_verified=False):
    """Iterate sessions using author scores and independently exported .t clocks."""
    digest = HC23_MAT_SHA256 if source_verified else sha256(mat)
    if digest != HC23_MAT_SHA256:
        raise ValueError("hc-23 processed MAT does not match Figshare file 14850551")
    raw = _raw_bouts(extracted)
    mat_rats = loadmat(mat, simplify_cells=True, variable_names=["rat"])["rat"]
    for animal in mat_items(mat_rats):
        slug = animal["name"]
        for j, session in enumerate(mat_items(animal["session"]), 1):
            day = mat_items(session["day"])[0]
            task, overnight = mat_items(day["task"])
            if task["name"] != "OpenField" or overnight["name"] != "Overnight":
                raise ValueError(f"Unexpected tasks for {slug} {j}")
            session_id = f"Session{j}"
            yield _hc23_session(slug, session_id, mat_items(overnight["bout"]), raw, digest)


def _hc23_session(subject, session, bouts, raw, digest):
    from audit_sleep_dandi import subtract

    qc = Counter()
    unit_table = defaultdict(list)
    selected_bouts = []
    state_source = defaultdict(float)
    for k, b in enumerate(bouts, 1):
        qc["published_bouts"] += 1
        scores = {"NREM": intervals(b["nRemTimes"]), "REM": intervals(b["remTimes"])}
        for state, ranges in scores.items():
            state_source[state] += sum(y - x for x, y in ranges)
        all_scores = sorted([*scores["NREM"], *scores["REM"]])
        if not all_scores:
            continue
        qc["scored_bouts"] += 1
        merged = {state: merge_same_state(ranges) for state, ranges in scores.items()}
        if any(merged[state] != ranges for state, ranges in scores.items()):
            qc["same_state_overlap_bouts"] += 1
        contested = clip(merged["NREM"], merged["REM"])
        if contested:
            qc["cross_state_overlap_bouts"] += 1
            qc["cross_state_excluded_seconds"] += float(sum(b-a for a,b in contested))
        scores = {state: subtract(ranges, contested) for state, ranges in merged.items()}
        raw_units = raw.get((subject, session, k))
        if not raw_units:
            qc["raw_missing_bouts"] += 1
            continue
        verified = {}
        ambiguous = set()
        for unit in mat_items(b["unit"]):
            qc["processed_unit_bout_rows"] += 1
            identity = np.atleast_1d(unit["ID"]).astype(float)
            cell_code = unit["type"]
            if (identity.shape != (2,) or not np.isfinite(identity).all()
                    or np.any(identity != np.floor(identity)) or cell_code not in (1, 2)):
                qc["invalid_mat_identity"] += 1
                continue
            idx = tuple(map(int, identity))
            t = np.atleast_1d(unit["spkTms"]).astype(float)
            cell_type = int(cell_code)
            raw_key = (*idx, cell_type)
            if raw_key not in raw_units or len(idx) != 2:
                qc["unmatched_unit_ids_or_type"] += 1
                continue
            raw_t = raw_units[raw_key]
            if raw_t is None:
                qc["zero_filled_tfile_rows"] += 1
                continue
            if not len(t) or not clean_clock(t) or not clean_clock(raw_t):
                qc["unit_clock_or_class_conflicts"] += 1
                continue
            ticks = np.rint(t * 1e4).astype(np.int64)
            if not np.allclose(ticks * 1e-4, t, rtol=0, atol=2e-8) or not np.isin(ticks, raw_t).all():
                qc["processed_raw_clock_conflicts"] += 1
                continue
            if idx in verified:
                if verified[idx] == cell_type:
                    qc["duplicate_mat_references_to_same_tfile"] += 1
                else:
                    ambiguous.add(idx)
                    qc["ambiguous_mat_unit_ids"] += 1
            verified[idx] = cell_type
        for uid in ambiguous:
            verified.pop(uid)
        qc["verified_unit_bout_rows"] += len(verified)
        if not verified:
            qc["unverified_bouts"] += 1
            continue
        pooled = np.unique(np.concatenate([raw_units[(*uid, kind)] for uid, kind in verified.items()])) * 1e-4
        if len(pooled) < 2:
            qc["empty_coverage_bouts"] += 1
            continue
        gaps = np.flatnonzero(np.diff(pooled) > 10)
        windows = [(float(c[0]), float(c[-1] + 1e-4)) for c in np.split(pooled, gaps+1) if len(c) > 1]
        scores = {state: clip(ranges, windows) for state, ranges in scores.items()}
        if not any(scores.values()):
            qc["scores_without_coverage"] += 1
            continue
        for uid, cell_type in verified.items():
            unit_table[uid].append(cell_type)
        selected_bouts.append((k, windows, scores, verified, raw_units))

    shared = set.intersection(*(set(b[3]) for b in selected_bouts)) if selected_bouts else set()
    ids = tuple(sorted(k for k in shared if len(set(unit_table[k])) == 1))
    cover = tuple(sorted(w for _, windows, _, _, _ in selected_bouts for w in windows))
    states = {state: tuple(sorted(t for _, _, scores, _, _ in selected_bouts for t in scores.get(state, ())))
              for state in ("NREM", "REM")}
    if overlaps(cover) or overlaps(sorted([*states["NREM"], *states["REM"]])):
        raise ValueError(f"Overlapping source bouts in {subject} {session}")
    trains = tuple(np.concatenate([raw_units[(*uid, kind[uid])] for _, _, _, kind, raw_units in selected_bouts]) * 1e-4 for uid in ids)
    if not all(clean_clock(train) for train in trains):
        raise ValueError(f"Sorted cluster clocks roll back across bouts in {subject} {session}")
    quality = {**qc, "source_state_seconds": dict(state_source), "retained_bouts": len(selected_bouts),
               "grid_units": sum(unit_table[uid][0] == 1 for uid in ids),
               "raw_member_status": "CRC-verified source archive; exact 100 microsecond timestamp concordance"}
    provenance = {"asset_id": f"hc-23:{subject}:{session}", "doi": "10.6084/m9.figshare.7973339.v1",
                  "source_file_id": 14850551, "sorted_archive_file_id": 14834216,
                  "path": f"rat/{subject}/{session}/Day1/Overnight", "sha256": digest,
                  "state_source": "Overnight remTimes/nRemTimes; two-rater video then author LFP thresholds"}
    clock = {"units": "seconds", "scale_from_tfile_ticks": .0001,
             "source": "author MAT against raw Neuralynx .t clocks; pooled support bounded per source sleep bout"}
    return SleepRecording(subject, session, "hc-23", "MEC II/III",
                          tuple(f"TT{t}_{c}" for t, c in ids), trains, states,
                          cover, clock, provenance, quality, "OVERNIGHT")


def _verified_record(spec, unit_spikes):
    evidence = spec["state_evidence"]
    if evidence["method"] not in ("author_scored", "independently_validated") or not evidence["reference"]:
        raise ValueError("Author scores or independently validated state labels are required")
    bounds = intervals(spec["recording_intervals"])
    state_ranges = {state: intervals(rows) for state, rows in spec["sleep_intervals"].items()}
    combined = sorted(row for ranges in state_ranges.values() for row in ranges)
    if (set(state_ranges) - set(STATES) or not combined or overlaps(combined) or overlaps(bounds)
            or clip(combined, bounds) != combined):
        raise ValueError("State evidence must be disjoint and inside checked recording intervals")
    ids, spikes = zip(*sorted(unit_spikes.items()), strict=True)
    if not all(clean_clock(t) and len(t) and t[0] >= bounds[0][0] and t[-1] <= bounds[-1][1]
               for t in spikes):
        raise ValueError("Unit clock conflicts with recording support")
    return SleepRecording(str(spec["subject_id"]), str(spec["session_id"]), spec["source_dataset"],
                          spec["region"], ids, spikes, state_ranges, tuple(bounds), spec["timebase"],
                          spec["provenance"], {"state_evidence": evidence}, spec["context"])


def _verified_input(directory, entry):
    source = (directory / entry["path"]).resolve(strict=True)
    if not source.is_relative_to(directory.resolve()) or sha256(source) != entry["sha256"]:
        raise ValueError("Source file or recorded digest does not match the audited input")
    return source


def hc24_matlab(directory, spec):
    """Load exact hc-24 cluster matrices identified in a verified source sidecar."""
    if spec["source_dataset"] != "hc-24":
        raise ValueError("Wrong adapter for the source")
    unit_spikes = {}
    for unit in spec["units"]:
        path = _verified_input(directory, unit)
        if "cl-" not in path.name:
            raise ValueError("Expected an author-sorted cl- unit, excluding tetrode multi-unit data")
        data = loadmat(path, variable_names=[unit["field"]])[unit["field"]]
        if data.ndim != 2 or data.shape[1] != 5 or unit["id"] in unit_spikes:
            raise ValueError("Unexpected hc-24 unit dimensions or identifier")
        unit_spikes[unit["id"]] = np.asarray(data[:, 4], dtype=float)
    return _verified_record(spec, unit_spikes)


def neurosuite(directory, spec):
    """Load paired .res/.clu after channel mapping and state scores are verified."""
    xml = ET.parse(_verified_input(directory, spec["xml"])).getroot()
    sample_hz = float(xml.findtext("./acquisitionSystem/samplingRate"))
    if sample_hz <= 0 or sample_hz != spec["timebase"]["spike_sample_hz"]:
        raise ValueError("NeuroSuite timebase differs from audited XML")
    unit_spikes = {}
    for group in spec["groups"]:
        ticks = np.loadtxt(_verified_input(directory, group["res"]), dtype=np.int64, ndmin=1)
        clu = np.loadtxt(_verified_input(directory, group["clu"]), dtype=np.int64, ndmin=1)[1:]
        if ticks.shape != clu.shape or ticks.ndim != 1 or np.any(ticks < 0) or np.any(np.diff(ticks) < 0):
            raise ValueError("NeuroSuite RES and CLU do not share one sorted clock")
        for cluster in np.unique(clu[clu >= 2]):
            uid = f"{group['shank']}:{cluster}"
            if uid in unit_spikes:
                raise ValueError("Duplicate NeuroSuite unit identifier")
            unit_spikes[uid] = ticks[clu == cluster] / sample_hz
    return _verified_record(spec, unit_spikes)


def audited_nwb(row, region, context, cfg):
    """Reuse DANDI's verified identifiers, clocks, and state interval extraction."""
    from audit_sleep_dandi import subtract
    from extract_sleep_spikes import fetch_json, open_asset, read_units

    if row["status"] != "READY" or region not in row["eligible_regions"]:
        raise ValueError("Expected an audited READY region")
    with open_asset(fetch_json(row["metadata_url"])) as source:
        units, times = read_units(source, row["dandiset"], cfg)
    unit_keys = ("id", "region", "eligible", "strict", "spikes", "duplicates", "nonmonotonic")
    keys = lambda values: [tuple(u[k] for k in unit_keys) for u in values]
    if keys(units) != keys(row["units"]):
        raise ValueError("NWB identity or timestamp QC changed after audit")
    state_ranges = {state: ranges["intervals"] for state, ranges in row["spontaneous_contexts"][context].items()
                    if state in STATES}
    unit_spikes = {str(u["id"]): t for u, t in zip(units, times, strict=True)
                   if u["eligible"] and u["region"] == region}
    bounds = subtract([row["clock"]["seconds"]], row["global_spike_gaps_over_10s"])
    prov = {k: row[k] for k in ("asset_id", "path", "version", "digest", "state_path", "nwb_id")}
    prov["sha256"] = row["digest"].get("dandi:sha2-256") or row["digest"]["dandi:dandi-etag"]
    prov["state_source"] = row["state_path"]
    spec = {"subject_id": row["subject"], "session_id": row["session"], "source_dataset": row["dandiset"],
            "region": region, "context": context, "timebase": row["clock"], "provenance": prov,
            "recording_intervals": bounds, "sleep_intervals": state_ranges,
            "state_evidence": {"method": "author_scored", "reference": row["state_path"]}}
    return _verified_record(spec, unit_spikes)
