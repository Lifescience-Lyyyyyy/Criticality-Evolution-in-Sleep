"""Read only the small spike, annotation and channel tables of published NWB assets."""

from contextlib import contextmanager
import json
import re
from urllib.request import urlopen

import h5py
import numpy as np
import remfile


DANDI = "https://api.dandiarchive.org/api/dandisets"
STATES = {"Awake": "Wake", "Non-REM": "NREM", "REM": "REM",
          "wake": "Wake", "awake": "Wake", "sws": "NREM", "nrem": "NREM", "non-rem": "NREM", "rem": "REM", "WSWake": "Wake",
          "WSNonREM": "NREM", "WSREM": "REM"}


def fetch_json(url):
    with urlopen(url, timeout=30) as response:
        return json.load(response)


def version_path(dandiset, version):
    if version == "draft":
        raise ValueError("Analysis requires a published immutable Dandiset version")
    return f"{DANDI}/{dandiset}/versions/{version}"


def listing(url):
    while url:
        page = fetch_json(url)
        yield from page["results"]
        url = page["next"]


@contextmanager
def open_asset(metadata):
    url, = [url for url in metadata["contentUrl"] if url.startswith("https://dandiarchive.s3.amazonaws.com/")]
    remote = remfile.File(url)
    try:
        with h5py.File(remote, "r") as file:
            yield file
    finally:
        remote.close()


def scalar_text(node):
    if isinstance(node, h5py.Dataset):
        return node.asstr()[()]
    return str(node)


def labels(table, name):
    return table[name].asstr()[:].tolist()


def bounds(table, scale=1):
    t = np.column_stack((table["start_time"][:], table["stop_time"][:])) / scale
    if not np.isfinite(t).all() or np.any(t[:, 1] <= t[:, 0]):
        raise ValueError(f"Nonfinite or nonpositive intervals at {table.name}")
    return t


def labeled_intervals(table, name="label", scale=1):
    t = bounds(table, scale)
    values = labels(table, name)
    if name + "_index" in table:
        splits = np.r_[0, table[name + "_index"][:]]
        values = ["+".join(values[a:b]) for a, b in zip(splits[:-1], splits[1:], strict=True)]
    return [[float(start), float(stop), label] for (start, stop), label in zip(t, values, strict=True)]


def state_table(file):
    for path in ("processing/behavior/sleep_states", "processing/behavior/SleepStates", "processing/behavior/states"):
        if path in file:
            return path, labeled_intervals(file[path])
    return "", []


def recording_clock(file):
    group = file.get("processing/ecephys/LFP")
    if group:
        for series in group.values():
            if "starting_time" in series:
                start = series["starting_time"]
                hz = float(start.attrs["rate"])
                nt = series["data"].shape[0]
                return {"source": start.name, "seconds": [float(start[()]), float(start[()]) + nt/hz], "rate_hz": hz}
            if "timestamps" in series:
                times = series["timestamps"]
                return {"source": times.name, "seconds": [float(times[0]), float(times[-1])], "rate_hz": None}
    return None


def electrode_locations(file, units):
    if "location" in units:
        return labels(units, "location")
    electrodes = file.get("general/extracellular_ephys/electrodes")
    if electrodes is None:
        return ["unknown"] * len(units["id"])
    by_channel = labels(electrodes, "location")
    if "electrodes" in units and "electrodes_index" in units:
        members = units["electrodes"][:]
        stop = units["electrodes_index"][:]
        start = np.r_[0, stop[:-1]]
        result = [
            "/".join(sorted({by_channel[i] for i in members[a:b]})) or "unknown"
            for a, b in zip(start, stop, strict=True)
        ]
        if any(x.lower() not in ("unknown", "") for x in result):
            return result
    if "electrode_group" in units:
        out = []
        for ref in units["electrode_group"][:]:
            if ref:
                group = file[ref]
                location = group.get("location", group.attrs.get("location", "unknown"))
                out.append(location.decode() if isinstance(location, bytes) else scalar_text(location))
            else:
                out.append("unknown")
        return out
    return ["unknown"] * len(units["id"])


def canonical_region(dandiset, label):
    if dandiset == "000061":
        key = label.upper()
        if key in {"CA1", "HPC", "HIP"}:
            return "CA1"
        if key in {"BLA", "BLV", "BLP"}:
            return "BLA"
        return None
    if dandiset == "000552":
        return "CA1" if label.upper() == "CA1" else None
    return label if label.upper() not in ("", "UNKNOWN", "CTX", "NAN", "UNDEFINED") and "/" not in label else None


def read_units(file, dandiset, rules, load_times=True):
    units = file.get("units")
    if units is None or not {"id", "spike_times", "spike_times_index"} <= set(units):
        return [], []
    ids = units["id"][:]
    if len(set(ids)) != len(ids):
        raise ValueError("Nonunique NWB unit IDs")
    index = units["spike_times_index"][:].astype(np.int64)
    stop = units["spike_times"].shape[0]
    if len(index) != len(ids) or not len(index) or index[-1] != stop or np.any(np.diff(index) < 0):
        raise ValueError("NWB unit ragged index does not match spike data")
    if not load_times:
        return [{"id": int(i), "spikes": int(n)} for i, n in zip(ids, np.diff(np.r_[0, index]), strict=True)], []
    all_spikes = units["spike_times"][:]
    if not np.isfinite(all_spikes).all():
        raise ValueError("Nonfinite spike timestamps")
    trains = np.split(all_spikes, index[:-1])
    regions = electrode_locations(file, units)
    class_col = next((n for n in ("cell_type", "putative_cell_type", "quality") if n in units), None)
    classes = labels(units, class_col) if class_col else ["unspecified"] * len(ids)
    rows = []
    for id_, region, group, spikes in zip(ids, regions, classes, trains, strict=True):
        isi = np.diff(spikes)
        duplicate = int(np.sum(isi == 0))
        reverse = int(np.sum(isi < 0))
        refractory = float(np.mean(isi < rules["refractory_seconds"])) if len(isi) and not reverse else None
        label = canonical_region(dandiset, region)
        valid = bool(len(spikes) > 1 and not reverse and not duplicate and spikes[0] >= 0 and label
                     and group.lower() not in {"mua","mu","noise","artifact","unsorted","unknown","bad"})
        rows.append({"id": int(id_), "region_raw": region, "region": label, "class": group,
                     "spikes": len(spikes), "first": float(spikes[0]) if len(spikes) else None,
                     "last": float(spikes[-1]) if len(spikes) else None, "duplicates": duplicate,
                     "nonmonotonic": reverse, "refractory_fraction": refractory,
                     "eligible": valid, "strict": bool(valid and refractory is not None and refractory <= rules["max_refractory_fraction"])})
    return rows, trains

def ripple_times(file):
    for path in ("processing/ecephys/Ripples", "processing/ecephys/Ca1Ripples", "stimulus/presentation/rip"):
        if path not in file:
            continue
        table = file[path]
        if path == "stimulus/presentation/rip":
            times = table["timestamps"][:]
            raw = labels(table, "data")
            parsed = [re.fullmatch(r"Ripple (start|peak|stop) (\d+)", name, re.I) for name in raw]
            if len(times) != len(parsed) or len(parsed) % 3 or any(m is None for m in parsed):
                raise ValueError("Unknown or incomplete ripple markers")
            for i in range(0,len(parsed),3):
                group = parsed[i:i+3]
                if ([m[1].lower() for m in group] != ["start","peak","stop"]
                    or len({m[2] for m in group}) != 1):
                    raise ValueError("Author ripple markers do not form chronological start/peak/stop triplets")
            # Event IDs can repeat in a file. The marked time triplet is its identity.
            t = times.reshape(-1,3)
        elif all(key in table for key in ("start_time", "peaks", "stop_time")):
            t = np.column_stack([table[key][:] for key in ("start_time", "peaks", "stop_time")])
        else:
            raise ValueError(f"Ripple source at {path} lacks annotated peaks and duration")
        if (not np.isfinite(t).all() or np.any(np.diff(t,axis=1) < 0) or np.any(t[:,0] >= t[:,2])
                or np.any(np.diff(t[:,1]) <= 0)):
            raise ValueError("Invalid or unordered ripple triplets")
        return path, t
    return "", np.empty((0,3))
