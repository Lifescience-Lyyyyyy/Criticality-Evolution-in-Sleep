#!/usr/bin/env python3
"""Verify the distributed atlas files against independent saved summaries and fit receipts."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory
from statistics import mean

from build_atlas_data import ROOT, DEFAULT, aggregate, export, read
from catalog import BASE, PAIRS, PIN, SOURCE_COMMIT, STATES, V1_FILE

PAIR_NAMES = {"NREM_minus_Wake": "NW", "REM_minus_NREM": "RN", "REM_minus_Wake": "RW"}
CATALOGS = (
    ("manifests/sleep_experiments.json", "file_sha256"),
    ("manifests/phase9_experiments.json", "result_sha256"),
    ("manifests/phase10_experiments.json", "result_sha256"),
)


def checked(condition, *args):
    if not condition:
        raise ValueError(args)


def close(a, b):
    return math.isclose(a, b, rel_tol=0, abs_tol=1e-12)


def load(p):
    return json.loads(p.read_text(encoding="utf-8"))


def validate(out):
    with TemporaryDirectory() as t:
        tmp = Path(t)
        export(tmp)
        actual = {f.name for f in out.glob("*.json")}
        checked(actual == {f.name for f in tmp.glob("*.json")} | {"anatomy.json"}, actual)
        for f in tmp.glob("*.json"):
            checked(f.read_bytes() == (out / f.name).read_bytes(), "regenerate", f.name)

    published = {}
    for path, key in CATALOGS:
        published.update(read(path)[key])
    initial = read("manifests/experiments.json")["experiments"]
    mouse_fit = next(r for r in initial if r["experiment_id"] == "mouse-44557104540")
    published.update(mouse_fit["results_sha256"] if isinstance(mouse_fit["results_sha256"], dict) else
                     {V1_FILE: mouse_fit["results_sha256"]})

    atlas = load(out / "atlas.json")
    anatomy = load(out / "anatomy.json")
    checked(atlas["baseline"] == SOURCE_COMMIT and atlas["schema"] == 1, atlas["baseline"])
    checked(anatomy["species"] == "mouse" and anatomy["voxelMm"] == .05, "anatomy")
    for key, p in atlas["provenance"].items():
        expected = hashlib.sha256((ROOT / p["path"]).read_bytes()).hexdigest()
        checked(expected == p["sha256"] == published.get(p["path"]), "fit source digest", p["path"], expected)
        checked(p["run"] == f"{BASE}/actions/runs/{key}" and len(p["commit"]) == 40, p["path"])
        checked((p["toolkit"], p["beta"], p["binSeconds"], p["units"]) == (PIN, 2, .04, "bits/s"), p["run"])
    contains = {}
    for species in ("mouse", "rat"):
        for region in load(out / f"{species}.json"):
            checked(region["species"] == species, region["id"])
            if species == "mouse":
                checked(anatomy["allenIds"][region["id"]]["acronym"] == region["allen"], region["id"])
                for plane in anatomy["views"].values():
                    pos = plane["positions"][region["id"]]
                    checked(region["id"] in plane["slices"][str(pos)]["labels"], region["id"], pos)
            else:
                checked(region["allen"] is None, "rat on mouse anatomy", region["id"])
            if not region["path"]:
                checked(not region["cohorts"], region["id"])
                continue
            data = load(out / region["path"])["cohorts"]
            for meta, cohort in zip(region["cohorts"], data, strict=True):
                id_ = cohort["id"]
                ds = atlas["datasets"][cohort["dataset"]]
                checked(id_ not in contains and ds["species"] == species, id_)
                checked(all(v == cohort[k] for k, v in meta.items()), id_)
                rec = cohort["records"]
                checked(cohort["sessions"] == len(rec) == len({r["id"] for r in rec}), id_)
                subjects = {r["subject"] for r in rec}
                checked(cohort["animals"] == (None if None in subjects else len(subjects)), id_)
                checked(cohort["states"] == [s for s in STATES if any(s in r["fits"]["reference"] for r in rec)], id_)
                checked(ds["version"] and ds["report"] and ds["paper"].startswith("https://doi.org/"), ds["id"])
                for r in rec:
                    checked(r["asset"] and r["digest"], id_, r["id"])
                    for variant, fits in r["fits"].items():
                        for state, f in fits.items():
                            checked(state in STATES and f["order"] == (10 if variant == "AR10" else 5), id_, variant, state)
                            checked(math.isfinite(f["value"]) and f["value"] >= 0 and f["provenance"] in atlas["provenance"],
                                    id_, r["id"], variant, state)
                            checked(f["root"] > 1 and f["exit"] == 0 and f["solver"] > 0, id_, r["id"])
                contains[id_] = cohort
    validate_preliminary(atlas, contains)
    return contains



def validate_preliminary(atlas, contains):
    """Cross-check displayed supplementary controls against the verified MATLAB fits."""
    extra=atlas.get("preliminary")
    if extra is None:
        return
    from preliminary_atlas_data import P111
    from collections import defaultdict
    from statistics import median
    checked(set(contains)<=set(extra["evidence"]), "cohort evidence coverage")
    fit_path=P111+"pairwise_equal_duration/matlab-fits.json"
    raw=read(fit_path)
    receipt=read(P111+"pairwise_equal_duration/matlab-fits.receipt.json")
    checked(receipt["fit_sha256"]==hashlib.sha256((ROOT/fit_path).read_bytes()).hexdigest(), "MATLAB receipt")
    checked(receipt["action"]==extra["pairwise"]["provenance"]["action"], "MATLAB link")
    checked(raw["toolkit_commit"]==PIN and raw["settings"]["beta"]==2, "estimator pin")
    qualified=read(P111+"pairwise_equal_duration/eligibility.json")
    by_case=defaultdict(dict)
    for row in raw["records"]:
        m,f=row["record"],row["fit"]
        checked(m["dandiset"]=="000939" and m["region"]=="PoSub" and m["context"]=="REST", m["asset_id"])
        checked(m["state"] not in by_case[m["variant"], m["asset_id"]],m["asset_id"],m["variant"])
        checked(f["segment_bins"]==[1500,1500,1500] and f["ar_order"]==5 and f["min_ar_root"]>1
                and f["solver_flag"]>0, m["asset_id"])
        by_case[m["variant"],m["asset_id"]][m["state"]]=(f["d2"],tuple(m["unit_ids"]),m["subject"])
    for name,(a,b) in PAIRS.items():
        comparison=extra["pairwise"]["contrasts"][name]
        cohort=contains[extra["pairwise"]["cohort"]]
        checked(comparison["legacy"]["n"]==cohort["effects"]["equal_duration"][name]["n"],name,"legacy")
        for variant,view,states in (("equal_duration_"+name,"corrected",{a,b}),
                                    ("equal_duration_WNR","three_state_fixed",set(STATES))):
            stats=comparison[view]
            animals=defaultdict(list)
            for case in qualified:
                computed=all(case["candidates"].get(s,0)>=3 for s in states)
                available=by_case.get((variant,case["asset_id"]),{})
                checked(set(available)==(states if computed else set()),name,variant,case["asset_id"])
                checked(computed==(not case["exclusions"][variant]),name,case["asset_id"])
                if available:
                    checked(available[a][1:]==available[b][1:],name,case["asset_id"])
                    animals[available[a][2]].append(available[b][0]-available[a][0])
            effects=[mean(v) for v in animals.values()]
            checked(stats["n"]==len(effects) and stats["sessions"]==sum(map(len,animals.values())),name,view)
            checked(close(stats["mean"],mean(effects)) and close(stats["median"],median(effects)),name,view)
            checked(stats["positive"]==sum(x>0 for x in effects) and stats["negative"]==sum(x<0 for x in effects),name,view)
    additional=P111+"sorting_checks/matlab-fits.json"
    ext=read(additional)
    ext_receipt=read(P111+"sorting_checks/matlab-fits.receipt.json")
    checked(ext_receipt["fit_sha256"]==hashlib.sha256((ROOT/additional).read_bytes()).hexdigest(),
            "sorting/temporal MATLAB receipt")
    checked(ext["toolkit_commit"]==PIN and len(ext["records"])==ext_receipt["count"],"sorting/temporal provenance")
    extensions={(r["record"]["dandiset"],r["record"]["region"],r["record"]["session"],
                 r["record"]["variant"],r["record"]["state"]): r
                for r in ext["records"]}
    examples={(s["cohort"],s["session"]):s for s in extra["sorting"]}
    checked(len(examples)==len(extra["sorting"]), "duplicate sorting cases")
    for example in extra["sorting"]:
        if example["run"]!=ext_receipt["action"]:
            continue
        _,region,dataset,_=example["cohort"].split("|")
        for variant in example["variants"]:
            for name,(a,b) in PAIRS.items():
                left,right=[extensions[dataset,region,example["session"],variant["id"],state] for state in (a,b)]
                checked(left["record"]["unit_ids"]==right["record"]["unit_ids"],"subset units",example["session"])
                checked(variant["units"]==len(left["record"]["unit_ids"]), "subset count",example["session"])
                checked(close(variant["effects"][name],right["fit"]["d2"]-left["fit"]["d2"]),
                        "subset effect",example["session"],name)
    for cohort,session in (("mouse|PoSub|000056|REST","Mouse28-140313"),
                           ("mouse|LDT|zenodo:19462601|REST","DL108_122320n1")):
        case=examples[cohort,session]
        checked(len(case["variants"])>=2 and case["run"].startswith(BASE+"/actions/runs/"),cohort,session)
        source=next(s for s in contains[cohort]["records"] if s["session"]==session)
        for k,(a,b) in PAIRS.items():
            ref=source["fits"]["reference"]
            checked(close(case["reference"]["effects"][k],ref[b]["value"]-ref[a]["value"]),session,k)

def verify_summaries(cohorts, out):
    by_session, by_animal = {}, {}
    for id_, cohort in cohorts.items():
        _, region, ds, context = id_.split("|")
        for variant in cohort["specs"]:
            for contrast, (a, b) in PAIRS.items():
                head = ds, region, context, variant, contrast
                rows = aggregate(cohort["records"], variant, [a, b])
                summary = cohort["effects"][variant][contrast]
                checked(len(rows) == summary["n"], head)
                if rows:
                    checked(close(mean([r["effect"] for r in rows]), summary["mean"]), head)
                for r in rows:
                    by_animal[(*head, r["id"])] = (r["effect"], r["sessions"])
                for s in cohort["records"]:
                    fits = s["fits"].get(variant, {})
                    if a in fits and b in fits:
                        by_session[(*head, s["subject"], s["id"])] = fits[b]["value"] - fits[a]["value"]
    for level, lookup in (("animal", by_animal), ("session", by_session)):
        with (ROOT / f"results/phase10/{level}_effects.csv").open(encoding="utf-8", newline="") as f:
            saved = {}
            for r in csv.DictReader(f):
                key = (*(r[k] for k in ("dataset", "region", "context", "variant", "contrast", "subject")),
                       *((r["asset_id"],) if level == "session" else ()))
                saved[key] = (float(r["effect"]), int(r["sessions"])) if level == "animal" else float(r["effect"])
        actual = {k: v for k, v in lookup.items() if k[0] in {"000056", "000939", "zenodo:19462601"}}
        checked(saved.keys() == actual.keys(), level, len(saved), len(actual))
        for key, value in saved.items():
            result = actual[key]
            checked(close(value, result) if level == "session" else close(value[0], result[0]) and value[1] == result[1], key)

    old = ("sleep_dandi/000041/cortex", "sleep_dandi/000061/cross", "sleep_dandi/000552/mouse-ca1",
           "sleep_crcns/hc-23/cohort")
    for file in old:
        src = read(f"results/{file}.summary.json")
        for r in src["per_animal_region"]:
            if r["contrast"] not in PAIR_NAMES:
                continue
            key = (r["dandiset"], r["region"], r["context"], r["variant"], PAIR_NAMES[r["contrast"]], r["subject"])
            result, n = by_animal[key]
            checked(n == r["sessions"] and close(result, r["effect"]), file, key)
        for r in src["session_contrasts"]:
            if r["contrast"] not in PAIR_NAMES or r["variant"] != "reference" or r["long_effect"] is None:
                continue
            species = "mouse" if r["dandiset"] == "000552" else "rat"
            cohort = cohorts["|".join([species, r["region"], r["dandiset"], r["context"]])]
            session = next(x for x in cohort["records"] if x["subject"] == r["subject"] and x["session"] == r["session"])
            key = (r["dandiset"], r["region"], r["context"], "AR10", PAIR_NAMES[r["contrast"]], r["subject"], session["id"])
            checked(close(by_session[key], r["long_effect"]), file, key)
    v1 = cohorts["mouse|V1|PNAS-V1|PUBLISHED"]
    vals = [by_session[("PNAS-V1", "V1", "PUBLISHED", "reference", "NW", None, s["id"])]
            for s in v1["records"]]
    original_v1 = read(V1_FILE)
    checked(all(close(a, b) for a, b in zip(vals, original_v1["paired_nrem_minus_wake"], strict=True))
            and v1["animals"] is None and v1["states"] == ["Wake", "NREM"], "V1")
    h24 = read("results/sleep_crcns/hc-24/summary.json")["reference"]
    for row in h24:
        contexts = [c for id_, c in cohorts.items() if id_.split("|")[1:3] == [row["region"], "hc-24"]]
        s = next(s for c in contexts for s in c["records"] if s["session"] == row["session"])
        checked(close(s["fits"]["reference"]["NREM"]["value"], row["ar5"]), row["session"])
        checked(close(s["fits"]["AR10"]["NREM"]["value"], row["ar10"]), row["session"])
    checked(all(p.stat().st_size < 1_000_000 for p in out.glob("*.json")), "individual data file size")
    print(f"PASS: {len(cohorts)} region contexts, {len(by_session)} session contrasts, {len(by_animal)} replicate contrasts.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=DEFAULT)
    a = p.parse_args()
    verify_summaries(validate(a.out), a.out)
