"""Repeat the 60 s, ten-cell comparisons at matched observed spike counts."""
import argparse
from collections import defaultdict
import gzip
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


SEEDS = (31400, 31401, 31402)


def unit(row):
    if row["model"] == "JercogEtAl2017":
        return str(row["connectivity_seed"]), ("figure7", "figure7_supp2"), row["context"]
    return row["animal"], ("Wake", "NREM"), row["state"]


def prepare(source, receipt, out):
    raw = source.read_bytes()
    assert sha256(raw).hexdigest() == json.loads(receipt.read_text())["paired_input_sha256"]
    rows = [r for r in json.loads(gzip.decompress(raw))["cases"]
            if r["category"] in ("empirical_match", "empirical")]
    assert len(rows) == 14
    assert all(r["n_units"] == 10 and r["duration_seconds"] == 60
               and r["bin_seconds"] == .04 and len(r["counts"]) == 1
               and len(r["counts"][0]) == 1500 for r in rows)
    spikes = {r["id"]: sum(r["counts"][0]) for r in rows}
    members = defaultdict(set)
    pair_targets = {}
    for row in rows:
        name, contrast, condition = unit(row)
        assert condition not in members[name] and condition in contrast
        members[name].add(condition)
        total = spikes[row["id"]]
        pair_targets[name] = min(pair_targets.get(name, total), total)
    assert len(members) == 7 and all(len(pair) == 2 for pair in members.values())
    common = min(spikes.values())
    assert common > 0
    cases = []
    for row in rows:
        name, _, _ = unit(row)
        x = np.asarray(row["counts"][0], dtype=np.int64)
        for mode, target in (("paired_rate", pair_targets[name]), ("common_rate", common)):
            for seed in SEEDS:
                source_seed = int.from_bytes(sha256(f"{row['id']}|{mode}|{seed}".encode()).digest()[:8],
                                             "little")
                y = np.random.default_rng(source_seed).multivariate_hypergeometric(x, target)
                np.testing.assert_array_less(y, x + 1)
                assert int(y.sum()) == target
                cases.append(dict(row, id=f"{row['id']}_{mode}_{seed}", counts=[y.tolist()],
                                  category="rate_sensitivity", source_id=row["id"],
                                  rate_control=mode, thinning_seed=seed, original_spikes=int(x.sum()),
                                  retained_spikes=target))
    doc = {"physiology_pass": {"verified_baseline": True}, "cases": cases}
    out.write_bytes(gzip.compress(json.dumps(doc, separators=(",", ":")).encode(), mtime=0))
    print(f"{len(cases)} inputs; common target {common} spikes over 60 s.")


def estimate_summary(values):
    good = [v for v in values if v is not None]
    return {"median": float(np.median(good)) if good else None,
            "range": [min(good), max(good)] if good else None,
            "passed": len(good), "attempted": len(values)}


def summarize(baseline, result, out):
    base = {r["id"]: r for r in json.loads(baseline.read_text())["records"]}
    rows = json.loads(result.read_text())["records"]
    assert len(rows) == 84
    observations = defaultdict(list)
    pairs = defaultdict(dict)
    for r in rows:
        name, contrast, condition = unit(r)
        assert r["source_id"] in base
        observations[r["source_id"], r["rate_control"]].append(r)
        k = r["model"], name, contrast, r["rate_control"], r["thinning_seed"]
        assert condition not in pairs[k]
        pairs[k][condition] = r
    obs = []
    for (key, mode), group in observations.items():
        assert len(group) == len(SEEDS)
        orig = base[key]
        item = group[0]
        obs.append({"source_id": key, "model": item["model"], "condition": unit(item)[2],
                    "rate_control": mode, "original_hz": item["original_spikes"] / 60,
                    "retained_hz": item["retained_spikes"] / 60,
                    "observed": {o: orig[o]["d2"] for o in ("fit", "ar10")},
                    "thinned": {o: estimate_summary([a[o]["d2"] for a in group])
                                for o in ("fit", "ar10")}})
    changes = defaultdict(lambda: defaultdict(list))
    for (model, name, (lo, hi), mode, _), group in pairs.items():
        assert group[lo]["retained_spikes"] == group[hi]["retained_spikes"]
        for order in ("fit", "ar10"):
            a, b = group[lo][order]["d2"], group[hi][order]["d2"]
            changes[model, name, mode][order].append(b-a if a is not None and b is not None else None)
    contrasts = [{"model": model, "recording": name, "rate_control": mode,
                  "delta": {o: estimate_summary(vals) for o, vals in orders.items()}}
                 for (model, name, mode), orders in changes.items()]
    out.write_text(json.dumps({"observations": obs, "contrasts": contrasts}, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "summarize"))
    parser.add_argument("source", type=Path)
    parser.add_argument("receipt_or_fits", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    action = prepare if args.mode=="prepare" else summarize
    action(args.source, args.receipt_or_fits, args.output)
