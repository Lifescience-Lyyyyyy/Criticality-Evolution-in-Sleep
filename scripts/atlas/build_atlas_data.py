#!/usr/bin/env python3
"""Export compact, species-specific atlas data from the committed MATLAB outputs."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

from catalog import BASE, CONTEXTS, DANDI_PAPERS, DATASETS, FIT_FILES, MOUSE, PAIRS, PIN, RAT, SOURCE_COMMIT, STATES, V1_FILE

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "website/public/data"
MAIN_CONTEXTS = {"HOME", "REST", "PRE", "POST", "POST_RUN", "OVERNIGHT", "PUBLISHED"}


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n", encoding="utf-8")


def blob(path):
    return f"{BASE}/blob/{SOURCE_COMMIT}/{path}"


def session_groups(sessions):
    groups = defaultdict(list)
    for s in sessions:
        groups[s["subject"] or s["session"]].append(s)
    return dict(sorted(groups.items()))


def aggregate(sessions, variant, states):
    """When linked, each replicate's state means use exactly its complete common sessions."""
    rows = []
    for name, group in session_groups(sessions).items():
        paired = [s for s in group if set(states) <= s["fits"].get(variant, {}).keys()]
        if paired:
            vals = {state: mean(s["fits"][variant][state]["value"] for s in paired) for state in states}
            rows.append({"id": name, "sessions": len(paired), "values": vals,
                         "effect": vals[states[1]] - vals[states[0]] if len(states) == 2 else None})
    return rows


def summarize(rows):
    vals = [r["effect"] for r in rows]
    return {"n": len(vals), "sessions": sum(r["sessions"] for r in rows),
            "mean": mean(vals) if vals else None, "positive": sum(x > 0 for x in vals),
            "negative": sum(x < 0 for x in vals)}


def finite(value):
    return isinstance(value, (float, int)) and math.isfinite(value)


def fit_row(fit, units, prov, order, dated):
    value = fit.get("d2", fit.get("d2_bits_per_second"))
    root = fit.get("min_ar_root", fit.get("fitted_min_root"))
    flag = fit["solver_flag"]
    assert finite(value) and value >= 0 and finite(root) and root > 1 and fit["exit_status"] == 0 and flag > 0
    residual = fit.get("residual_acf")
    max_res = max((abs(x) for x in residual if finite(x)), default=None) if residual else None
    kurt = fit.get("residual_excess_kurtosis")
    rmse = fit.get("acf_rmse")
    warnings = []
    if finite(kurt) and abs(kurt) > 1:
        warnings.append("kurtosis")
    if finite(max_res) and max_res > .1:
        warnings.append("residualACF")
    if dated and finite(rmse) and rmse > .1:
        warnings.append("modelACF")
    return {"value": value, "order": order, "units": units, "seconds": fit["duration_s"],
            "root": root, "exit": fit["exit_status"], "solver": flag, "kurtosis": kurt,
            "residual": max_res, "rmse": rmse, "silence": fit.get("max_silence_s"),
            "rate": fit.get("rate_hz"), "warnings": warnings,
            "provenance": prov}


def dataset_meta():
    src = read("manifests/sources.json")
    p8 = read("manifests/sleep_sources.json")["dandisets"]
    p9 = read("manifests/phase9_sources.json")
    p10 = read("manifests/phase10_sources.json")["sources"]
    sources = {}
    for ds, (species, title, report, recording, scoring, wake, limits) in DATASETS.items():
        if ds == "PNAS-V1":
            version = "v1"
            url = f"https://doi.org/{src['mouse']['doi']}"
            paper = f"https://doi.org/{src['estimator']['paper_doi']}"
            manifest, access = "manifests/sources.json", "公开 Figshare"
        elif ds in p8:
            version = p8[ds]["version"]
            url, paper = f"https://doi.org/10.48324/dandi.{ds}/{version}", DANDI_PAPERS[ds]
            manifest, access = "manifests/sleep_sources.json", "公开 NWB"
        elif ds in p9:
            version = "Figshare v1" if ds == "hc-23" else "CRCNS 压缩包（发布页无版本号）"
            url, paper = f"https://doi.org/{p9[ds]['doi']}", p9[ds]["paper"]
            manifest, access = "manifests/phase9_sources.json", "公开衍生档案" if ds == "hc-23" else "CRCNS 原档需账户"
        else:
            version, url, paper = (p10[ds][key] for key in ("version", "release", "paper"))
            manifest, access = "manifests/phase10_sources.json", "公开 Zenodo" if ds.startswith("zenodo:") else "公开 NWB"
        sources[ds] = dict(id=ds, species=species, title=title, version=version, url=url, paper=paper,
                           report=blob(f"reports/{report}.md"), manifest=blob(manifest),
                           recording=recording, scoring=scoring, wake=wake, limitations=limits, access=access)
    return sources


def build():
    sources = dataset_meta()
    provenances = {}
    cohorts = defaultdict(dict)
    unit_sets = {}
    for path in (*FIT_FILES, V1_FILE):
        data = read(path)
        config = data["settings"]
        assert (config["bin_seconds"], config["ar_order"], config["beta"], config["fit_method"]) == (.04, 5, 2, "YuleWalker")
        assert data["toolkit_commit"] == PIN and "R2024b" in data["matlab"]
        run = str(data["run_id"])
        assert run not in provenances and run.isdecimal()
        provenances[run] = dict(
            run=f"{BASE}/actions/runs/{run}", file=blob(path), path=path,
            sha256=hashlib.sha256((ROOT / path).read_bytes()).hexdigest(),
            commit=data["project_commit"], toolkit=PIN, matlab=data["matlab"],
            binSeconds=.04, beta=2, method="Yule–Walker", units="bits/s",
        )
        if path == V1_FILE:
            prep = {s["session"]: s for s in read("manifests/mouse-preparation.json")["sessions"]}
            replay = {(s["session"], s["state"]): s for s in data["optimizer_replay"]}
            for f in data["primary"]:
                s = f["session"]
                verified = replay[(s, f["state"])]
                assert abs(verified["d2_bits_per_second"] - f["d2_bits_per_second"]) < 1e-8
                cohort = cohorts["mouse", "V1", "PNAS-V1", "PUBLISHED"]
                row = cohort.setdefault(s, dict(id=s, session=s, asset=prep[s]["file"], digest=prep[s]["md5"],
                                                subject=None, fits={"reference": {}}))
                state = f["state"].title() if f["state"] == "WAKE" else f["state"]
                assert state not in row["fits"]["reference"]
                row["fits"]["reference"][state] = fit_row(
                    {**f, "solver_flag": verified["solver_flag"], "duration_s": f["samples"] * .04},
                    prep[s]["units"], run, 5, False,
                )
            continue
        entries = data["records"] if isinstance(data["records"], list) else [data["records"]]
        for row in entries:
            r, fit = row["record"], row["fit"]
            ds, area, ctx, variant, state = (r[k] for k in ("dandiset", "region", "context", "variant", "state"))
            assert r["bin_seconds"] == .04 and state in STATES and ctx in CONTEXTS
            assert variant in ("reference", "equal_duration", "edge_trim", "strict_units")
            dataset = sources[ds]
            if "version" in r:
                assert r["version"] == dataset["version"], (ds, r["version"], dataset["version"])
            assert r["subject"] and r["asset_id"] and r["digest"] and r["state_path"]
            ids = r["unit_ids"] if isinstance(r["unit_ids"], list) else [r["unit_ids"]]
            assert len(set(ids)) == len(ids) and ids
            # The frontal roll-up reuses each recording exactly once.
            for region in (area, "Frontal") if ds == "000041" else (area,):
                cohort_key = (dataset["species"], region, ds, ctx)
                group = cohorts[cohort_key]
                a = r["asset_id"]
                session = group.setdefault(a, dict(id=a, subject=r["subject"], session=r["session"],
                                                   asset=r["path"], digest=r["digest"], fits={}))
                assert (session["subject"], session["session"], session["asset"]) == (r["subject"], r["session"], r["path"])
                sig = (*cohort_key, a, variant)
                assert unit_sets.setdefault(sig, ids) == ids, (sig, state)
                models = [(variant, fit, 5)]
                if variant == "reference" and row.get("long_reach"):
                    models.append(("AR10", row["long_reach"], 10))
                for model, result, order in models:
                    assert result["ar_order"] == order
                    target = session["fits"].setdefault(model, {})
                    assert state not in target, (sig, state)
                    entry = fit_row(result, len(ids), run, order, ds in {"000041", "000061", "000552"})
                    entry["score"] = r["state_path"]
                    target[state] = entry
    controls = read_controls(provenances)
    regions = all_regions()
    index = defaultdict(list)
    for (species, region, ds, ctx), rows in sorted(cohorts.items()):
        assert (species, region) in regions, (species, region, ds, ctx)
        group = [s for _, s in sorted(rows.items(), key=lambda p: (p[1]["subject"] or "", p[1]["session"]))]
        specs = sorted(set().union(*(s["fits"].keys() for s in group)))
        assert "reference" in specs and all("reference" in s["fits"] for s in group)
        states = [s for s in STATES if any(s in row["fits"]["reference"] for row in group)]
        id_ = "|".join((species, region, ds, ctx))
        effects = {variant: {key: summarize(aggregate(group, variant, pair)) for key, pair in PAIRS.items()}
                   for variant in specs}
        subjects = {s["subject"] for s in group}
        units = [x["units"] for s in group for x in s["fits"]["reference"].values()]
        index[species, region].append(dict(
            id=id_, dataset=ds, context=ctx, contextLabel=CONTEXTS[ctx], animals=None if None in subjects else len(subjects),
            sessions=len(group), states=states, specs=specs, effects=effects, unitRange=[min(units), max(units)],
            controls=controls.get((ds, region, ctx), []), records=group,
        ))
    for (species, region), cohorts_for_region in index.items():
        # Default is the main context with the most verified animals, then most sessions, then stable ID.
        cohorts_for_region.sort(key=lambda c: (c["context"] not in MAIN_CONTEXTS, -(c["animals"] or 0), -c["sessions"], c["id"]))
        # Keep PRE first if its coverage ties with POST.
        if species == "rat" and region in {"CA1", "BLA"}:
            cohorts_for_region.sort(key=lambda c: {"PRE": 0, "POST": 1}.get(c["context"], 2))
    return regions, sources, provenances, index


def read_controls(provenances):
    phase10 = read("manifests/phase10_experiments.json")
    repeated = defaultdict(list)
    csv_path = "results/phase10/control_effects.csv"
    run = str(next(r["run_id"] for r in phase10["runs"] if r["batch"] == "matched-controls"))
    with (ROOT / csv_path).open(newline="", encoding="utf-8") as f:
        by_anchor = defaultdict(dict)
        for row in csv.DictReader(f):
            key = tuple(row[k] for k in ("dataset", "region", "subject", "asset_id", "contrast"))
            by_anchor[key][row["variant"]] = row
        for (ds, region, subject, asset, contrast), values in sorted(by_anchor.items()):
            ref = values["matched_controls_reference"]
            for model in ("unit_count", "rate_thinning"):
                r = values[model]
                repeated[ds, region, "REST"].append(dict(
                    model=model, contrast=contrast, subject=subject, asset=asset, reference=float(ref["median"]),
                    value=float(r["median"]), repeats=int(r["valid_repeats"]), range=json.loads(r["quantiles"]),
                    run=f"{BASE}/actions/runs/{run}", source=blob(csv_path),
                ))
    control_sources = (
        ("results/sleep_dandi/controls.summary.json", "rate_thinning", "NW"),
        ("results/sleep_dandi/half-count.summary.json", "unit_count", "NW"),
        ("results/sleep_crcns/hc-23/controls.summary.json", "rate_thinning", "RN"),
        ("results/sleep_crcns/hc-23/controls.summary.json", "unit_count", "RN"),
    )
    for path, model, contrast in control_sources:
        doc = read(path)
        for row in doc["sensitivity"]:
            ds, region, ctx = (row[k] for k in ("dandiset", "region", "context"))
            if model not in row or not row[model]["effect"]["n"]:
                continue
            v = row[model]["effect"]
            for area in (region, "Frontal") if ds == "000041" else (region,):
                repeated[ds, area, ctx].append(dict(
                    model=model, contrast=contrast, subject=row["subject"], asset=row["session"],
                    reference=row["reference_effect"], value=v["median"], repeats=v["n"], range=[v["low"], v["high"]],
                    run=f"{BASE}/actions/runs/{doc['run_id']}", source=blob(path),
                ))
    return repeated


def all_regions():
    regions = {}
    for code, name, english, system, allen, about in MOUSE:
        regions["mouse", code] = dict(id=code, name=name, english=english, system=system,
                                       allen=allen, about=about, species="mouse")
    for code, name, english, system in RAT:
        regions["rat", code] = dict(id=code, name=name, english=english, system=system,
                                     allen=None, about=None, species="rat")
    return regions


def export(out):
    regions, datasets, provenance, index = build()
    for ds, entry in datasets.items():
        cohorts = [c for cs in index.values() for c in cs if c["dataset"] == ds and c["context"] in MAIN_CONTEXTS]
        subjects = {r["subject"] for c in cohorts for r in c["records"]}
        entry["animals"] = None if None in subjects else len(subjects)
        entry["sessions"] = len({r["id"] for c in cohorts for r in c["records"]})
    for species in ("mouse", "rat"):
        contents = []
        for (sp, region), meta in regions.items():
            if sp != species:
                continue
            cohorts = index[sp, region]
            if cohorts:
                slug = region.replace(" ", "-").replace("/", "-")
                path = f"{species}-{slug}.json"
                write(out / path, dict(cohorts=cohorts))
                summary = [{k: v for k, v in x.items() if k not in ("controls", "records")} for x in cohorts]
            else:
                path, summary = None, []
            contents.append({**meta, "path": path, "cohorts": summary})
        write(out / f"{species}.json", contents)
    audit = ROOT / "results/validation/impact_assessment.json"
    write(out / "atlas.json", dict(
        schema=1, baseline=SOURCE_COMMIT, datasets=datasets, provenance=provenance, states=STATES, pairs=PAIRS,
        audit=read(audit) if audit.exists() else [],
    ))
    return regions, datasets, provenance, index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT)
    args = parser.parse_args()
    _, datasets, _, index = export(args.out)
    for ds, meta in datasets.items():
        print(ds, meta["species"], meta["animals"], meta["sessions"])
    print(f"Exported {sum(map(len, index.values()))} separate region/recording contexts to {args.out}")


if __name__ == "__main__":
    main()
