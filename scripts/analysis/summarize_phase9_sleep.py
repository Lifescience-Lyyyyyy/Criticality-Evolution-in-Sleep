"""Summarize CRCNS controls and compare verified sleep contrasts without pooling datasets."""

from collections import defaultdict
from pathlib import Path
from statistics import mean
import argparse
import json

from summarize_sleep_dandi import distribution, state_summary, valid


def control_summary(rows):
    grouped = defaultdict(dict)
    keys = ("dandiset", "subject", "session", "region", "context")
    for r in rows:
        m = r["record"]
        key = tuple(m[k] for k in keys)
        job = (m["state"], m["variant"], m["iteration"])
        if job in grouped[key]:
            raise ValueError(f"Duplicate control: {key} {job}")
        grouped[key][job] = r
    output = []
    for key, jobs in sorted(grouped.items()):
        pair = ("NREM", "Wake") if ("Wake", "matched_controls_reference", 0) in jobs else ("REM", "NREM")
        high, low = pair
        originals = {s: jobs[s, "matched_controls_reference", 0] for s in pair}
        baseline = [originals[s]["fit"] for s in pair]
        result = {**dict(zip(keys,key,strict=True)), "contrast": f"{high}_minus_{low}",
                  "reference_effect": baseline[0]["d2"]-baseline[1]["d2"] if all(map(valid,baseline)) else None,
                  "reference_units": len(originals[high]["record"]["unit_ids"])}
        for variant in ("unit_count", "rate_thinning"):
            by_rep = defaultdict(dict)
            for (state, v, n), row in jobs.items():
                if v == variant:
                    by_rep[n][state] = row
            if not by_rep:
                continue
            paired = []
            for n, pair_results in sorted(by_rep.items()):
                fits = [pair_results[s]["fit"] for s in pair]
                if all(map(valid, fits)):
                    paired.append(fits[0]["d2"]-fits[1]["d2"])
            result[variant] = {"effect": distribution(paired), "state_diagnostics": {}}
            for state in pair:
                fit_rows = [v[state] for v in by_rep.values() if valid(v[state]["fit"])]
                original = originals[state]["fit"]
                result[variant]["state_diagnostics"][state] = {
                    "assessed": len(fit_rows), "seeds": sorted({r["record"]["seed"] for r in fit_rows}),
                    "unit_count": fit_rows[0]["record"]["unit_count"] if fit_rows else None,
                    "actual_rate_hz": distribution([r["fit"]["rate_hz"] for r in fit_rows]),
                    "variance_change": distribution([r["fit"]["count_variance"]-original["count_variance"] for r in fit_rows]),
                    "lag1_change": distribution([r["fit"]["observed_acf"][0]-original["observed_acf"][0] for r in fit_rows]),
                }
        output.append(result)
    return output


def comparison(pnas_path, cortex_path, mouse_path, cross_path, hc23_path):
    groups = []
    pnas = json.loads(pnas_path.read_text())
    sessions = defaultdict(dict)
    for r in pnas["primary"]:
        if r["exit_status"] != 0:
            raise ValueError("Invalid original V1 reference")
        sessions[r["session"]][r["state"]] = r["d2_bits_per_second"]
    effects = {s: r["NREM"]-r["WAKE"] for s,r in sorted(sessions.items())}
    if len(effects) != 13:
        raise ValueError("V1 session inventory differs from the validated replication")
    groups.append({"source": str(pnas_path), "dataset": "PNAS mouse V1", "region": "V1", "species":"Mouse", "context":"HOME",
                   "contrast":"NREM_minus_Wake", "replicate":"session", "animals":None,
                   "sessions":len(effects), "effect_by_replicate":effects,
                   "effect_mean":mean(effects.values()), "positive":sum(x>0 for x in effects.values())})
    for source, species, collapse in ((cortex_path,"Rat","Frontal cortex"),(mouse_path,"Mouse",None),
                                      (cross_path,"Rat",None),(hc23_path,"Rat",None)):
        d = json.loads(source.read_text())
        field = "per_animal_dataset" if collapse else "per_animal_region"
        cells = defaultdict(list)
        for r in d[field]:
            if r["variant"] == "reference":
                cells[r["dandiset"], collapse or r["region"], r["context"], r["contrast"]].append(r)
        for (dataset, region, context, contrast), rows in sorted(cells.items()):
            values = {r["subject"]: r["effect"] for r in rows}
            groups.append({"source":str(source),"dataset":dataset,"region":region,"species":species,
                           "context":context,"contrast":contrast,"replicate":"animal","animals":len(values),
                           "sessions":sum(r["sessions"] for r in rows),"effect_by_replicate":values,
                           "effect_mean":mean(values.values()),"positive":sum(x>0 for x in values.values())})
    return groups


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("--comparison-root",type=Path)
    args=p.parse_args()
    data=json.loads(args.input.read_text())
    summary={"schema_version":1, "run_id":data["run_id"]}
    if data["input"]["mode"] == "controls":
        summary["sensitivity"] = control_summary(data["records"])
    else:
        summary.update(state_summary(data["records"]))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    if args.comparison_root:
        root=args.comparison_root
        base=root/"results/sleep_dandi"
        paths=(root/"results/pnas/mouse.json",base/"000041/cortex.summary.json",
               base/"000552/mouse-ca1.summary.json",base/"000061/cross.summary.json")
        groups=comparison(*paths,args.output)
        (args.output.parent/"phase9_comparison.json").write_text(json.dumps(groups,indent=2)+"\n")
