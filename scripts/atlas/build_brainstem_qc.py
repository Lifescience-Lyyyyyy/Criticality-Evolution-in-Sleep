"""Build the WT PAG/LDT contamination sensitivities at the original state intervals."""
import csv
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from recheck_sources import matlab_segments, refs_by_asset


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def brainstem_tasks(directory, chosen, original):
    refs = refs_by_asset()
    output, receipt = [], []
    for case in (c for c in chosen if c["case"] in ("F", "G")):
        asset, region = case["asset_id"], case["regions"][0]
        audit = original[asset]
        base = directory / audit["NP_PATH"]
        arrays, sorting, locs, coarse = (base / f for f in ("1k_train.npz", "cluster_info.TSV", "channel_locations.json", "traind.csv"))
        assert digest(arrays) == audit["digest"]["extracted_sha256"]
        assert digest(coarse) == audit["companion_sha256"]
        loc = json.loads(locs.read_text())
        with sorting.open(newline="") as f:
            clusters = list(csv.DictReader(f, delimiter="\t"))
        original_fits = [x for x in refs[asset] if x["region"] == region]
        assert len(original_fits) == 3
        cells = original_fits[0]["unit_ids"]
        original_regions = {c["cluster_id"] + "_good": (c.get("group") or c["KSLabel"],
                            loc.get("channel_" + c["ch"], {}).get("brain_region"), float(c["ContamPct"]))
                            for c in clusters}
        assert all(original_regions[x][0:2] == ("good", region) for x in cells)
        expected_spikes = audit["clock"]["seconds"][1] * 1000
        assert expected_spikes.is_integer()
        source_end = int(expected_spikes)
        unit_bins = {}
        same = hashlib.sha256()
        with np.load(arrays, allow_pickle=False) as z:
            for uid in cells:
                signal = z[uid][:source_end]
                assert len(signal) == source_end
                t = np.flatnonzero(signal) * .001
                unit_bins[uid] = t
                same.update(uid.encode()+b"\0")
                same.update(np.asarray(t, dtype="<f8").tobytes())
        assert same.hexdigest() == original_fits[0]["extracted_unit_sha256"]
        cuts = []
        for name, q in (("contamination_q75", .75), ("contamination_q50", .50)):
            cutoff = float(np.quantile([original_regions[x][2] for x in cells], q))
            accept = [x for x in cells if original_regions[x][2] <= cutoff]
            assert len(accept) >= 5
            cuts.append(dict(variant=name, cutoff_pct=cutoff, retained=len(accept)))
            t = np.sort(np.concatenate([unit_bins[x] for x in accept]))
            for r in original_fits:
                spans = [(round(x/.04), round(y/.04)) for x, y in r["intervals_s"]]
                output.append(({**r, "variant": name, "unit_ids": accept, "contamination_cutoff_pct": cutoff,
                                "audit_case": case["case"]}, matlab_segments(t, spans, .04)))
        anatomical = [c for c in clusters if loc.get("channel_" + c["ch"], {}).get("brain_region") == region]
        types = {}
        for c in anatomical:
            k = c.get("group") or c["KSLabel"]
            types[k] = 1 + types.get(k, 0)
        receipt.append(dict(case=case["case"], asset=asset, region=region, original_source=audit["path"],
                            sorting_sha256=digest(sorting), locations_sha256=digest(locs),
                            high_resolution_sha256=audit["digest"]["extracted_sha256"],
                            unit_count_by_sort_label=types, accepted_original=len(cells), subsets=cuts))
    assert len(output) == 12
    return output, receipt
