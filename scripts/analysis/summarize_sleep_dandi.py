"""Summarize verified fits by session and animal, retaining unmatched outcomes."""

from collections import defaultdict
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.stats import rankdata


COMPARISONS = (("NREM", "Wake"), ("REM", "NREM"), ("REM", "Wake"))
CONTEXTS = ("HOME", "PRE", "POST")


def valid(fit):
    return bool(fit and fit["d2"] is not None and fit["d2"] >= -1e-10 and fit["exit_status"] == 0
                and fit["solver_flag"] is not None and fit["solver_flag"]>0
                and fit["min_ar_root"]>1 and fit["innovation_variance"]>0)


def distribution(values):
    if not len(values):
        return {"n":0,"low":None,"median":None,"high":None}
    lo,mid,hi = np.quantile(values,[.025,.5,.975])
    return {"n":len(values),"low":float(lo),"median":float(mid),"high":float(hi)}


def assess(fit):
    acf = [abs(a) for a in fit["residual_acf"] if a is not None]
    r_max = max(acf, default=None)
    return {"max_abs_residual_acf_20":r_max,
            "model_deviation":bool((r_max is not None and r_max>.1)
              or (fit["residual_excess_kurtosis"] is not None and abs(fit["residual_excess_kurtosis"])>1)
              or (fit["acf_rmse"] is not None and fit["acf_rmse"]>.1))}


def state_summary(rows):
    grouped = defaultdict(dict)
    for r in rows:
        record = r["record"]
        key = tuple(record[k] for k in ("dandiset","subject","session","context","region","variant"))
        if record["state"] in grouped[key]:
            raise ValueError(f"Duplicate job for {key} {record['state']}")
        grouped[key][record["state"]] = r
    contrasts = []
    for (ds,bird,session,context,region,variant), states in grouped.items():
        for high,low in COMPARISONS:
            if not (high in states and low in states):
                continue
            pair = [states[s] for s in (high,low)]
            if not all(valid(q["fit"]) for q in pair):
                continue
            first,second = (q["fit"] for q in pair)
            alts = [q["long_reach"] for q in pair]
            contrasts.append({"dandiset":ds,"subject":bird,"session":session,"context":context,"region":region,
                              "variant":variant,"contrast":f"{high}_minus_{low}","effect":first["d2"]-second["d2"],
                              "rate_difference_hz":first["rate_hz"]-second["rate_hz"],"seconds":[q["duration_s"] for q in (first,second)],
                              "long_effect":alts[0]["d2"]-alts[1]["d2"] if all(valid(q) for q in alts) else None,
                              "model_deviation":any(assess(q)["model_deviation"] for q in (first,second))})
    per_animal = defaultdict(list)
    overall = defaultdict(list)
    joint = defaultdict(dict)
    for item in contrasts:
        per_animal[tuple(item[k] for k in ("dandiset","context","subject","region","variant","contrast"))].append(item)
        overall[tuple(item[k] for k in ("dandiset","context","subject","variant","contrast"))].append(item)
        if item["dandiset"]=="000061" and item["region"] in ("CA1","BLA"):
            joint[tuple(item[k] for k in ("subject","session","context","variant","contrast"))][item["region"]]=item
    akeys = ("dandiset","context","subject","region","variant","contrast")
    okeys = ("dandiset","context","subject","variant","contrast")
    flatten = lambda keys, groups:[{**dict(zip(keys,key,strict=True)),"effect":float(np.mean([r["effect"] for r in items])),
                                    "sessions":len(items)} for key,items in sorted(groups.items())]
    aligned = []
    for (subject,session,context,variant,contrast), pair in joint.items():
        if set(pair)!={"CA1","BLA"}:
            continue
        a,b = pair["CA1"],pair["BLA"]
        if a["seconds"]!=b["seconds"]:
            raise ValueError(f"Nonmatching {subject} {session} {context} cross-region exposure")
        aligned.append({"subject":subject,"session":session,"context":context,"variant":variant,
                        "contrast":contrast,"ca1":a["effect"],"bla":b["effect"],"ca1_minus_bla":a["effect"]-b["effect"]})
    return {"session_contrasts":contrasts,"per_animal_region":flatten(akeys,per_animal),
            "per_animal_dataset":flatten(okeys,overall),"paired_cross_region_sessions":aligned}


def spearman(x,y):
    if len(x)<5 or np.ptp(x)==0 or np.ptp(y)==0:
        return None
    return float(np.corrcoef(rankdata(x),rankdata(y))[0,1])


def event_summary(rows,seed):
    groups = defaultdict(list)
    full_pairs = {}
    for r in rows:
        m=r["record"]
        if m["variant"]=="swr_stability" and valid(r["fit"]):
            key=tuple(m[k] for k in ("dandiset","subject","session","context","pair_start_s"))
            full_pairs[key]=r["fit"]["d2"]
        if m["variant"]=="swr_window" and valid(r["fit"]):
            groups[tuple(m[k] for k in ("dandiset","subject","session","context"))].append(r)
    output=[]
    for key,subset in sorted(groups.items()):
        times=np.array([x["record"]["start_s"] for x in subset])
        d2=np.array([x["fit"]["d2"] for x in subset])
        rates=np.array([x["fit"]["rate_hz"] for x in subset])
        swr=np.array([x["record"]["swr_count"] for x in subset])
        halves = defaultdict(list)
        for j,row in enumerate(subset):
            parent = row["record"].get("pair_start_s")
            if parent is not None:
                halves[parent].append(d2[j])
        half_diffs=[abs(y-x) for pair in halves.values() if len(pair)==2 for x,y in (pair,)]
        m=dict(zip(("dandiset","subject","session","context"),key,strict=True))
        m["windows"]=len(subset)
        m["adjacent_window_abs_delta"]=distribution(half_diffs)
        m["full_vs_mean_halves_abs_delta"]=distribution(
            [abs(full_pairs[(*key,at)]-np.mean(pair)) for at,pair in halves.items()
             if len(pair)==2 and (*key,at) in full_pairs])
        if len(subset)<12 or np.ptp(swr)==0:
            output.append({**m,"status":"INSUFFICIENT_ELIGIBLE_WINDOWS"})
            continue
        lo,hi=np.quantile(swr,[1/3,2/3])
        if lo>=hi:
            output.append({**m,"status":"TIED_THRESHOLDS"})
            continue
        rlow=np.flatnonzero(swr<=lo)
        rhigh=np.flatnonzero(swr>=hi)
        candidates=[]
        for i in rhigh:
            for j in rlow:
                r=min(rates[i],rates[j])
                gap=abs(times[i]-times[j])
                if r>0 and max(rates[i],rates[j])/r<=1.10 and gap<=1200:
                    candidates.append((max(rates[i],rates[j])/r,gap,i,j))
        used=set()
        matched=[]
        for _,gap,i,j in sorted(candidates):
            if i not in used and j not in used:
                used.update([i,j])
                matched.append({"rich_start_s":float(times[i]),"poor_start_s":float(times[j]),"time_gap_s":float(gap),
                                "rich_rate_hz":float(rates[i]),"poor_rate_hz":float(rates[j]),"rich_minus_poor":float(d2[i]-d2[j]),
                                "high_count_fraction_gap":abs(subset[i]["record"]["high_count_bin_fraction"]-
                                                              subset[j]["record"]["high_count_bin_fraction"])})
        rng=np.random.default_rng(seed)
        blocks=np.floor(times/300).astype(int)
        indices=[np.flatnonzero(blocks==b) for b in np.unique(blocks)]
        metrics={"swr_count":swr,"swr_occupancy":np.array([x["record"]["swr_occupancy"] for x in subset]),
                 "rate_hz":rates,"fano":np.array([x["fit"]["fano"] for x in subset]),
                 "burst_fraction":np.array([x["record"].get("high_count_bin_fraction",0) for x in subset])}
        associations={}
        for name,values in metrics.items():
            p=spearman(values,d2)
            boot=[]
            if len(indices)>=5 and p is not None:
                for _ in range(500):
                    selection=np.concatenate([indices[i] for i in rng.integers(0,len(indices),len(indices))])
                    value=spearman(values[selection],d2[selection])
                    if value is not None:
                        boot.append(value)
            associations[name]={"rho":p,"blocks":len(indices),"time_block_interval":distribution(boot)}
        output.append({**m,"status":"DESCRIPTIVE" if len(matched)>=3 else "INSUFFICIENT_MATCHED_PAIRS",
                       "thresholds":[float(lo),float(hi)],"matched_pairs":matched,
                       "unmatched_median_difference":float(np.median(d2[rhigh])-np.median(d2[rlow])),
                       "matched_difference":distribution([r["rich_minus_poor"] for r in matched]),
                       "matched_with_high_count_fraction_gap_le_0p05":distribution(
                           [r["rich_minus_poor"] for r in matched if r["high_count_fraction_gap"]<=.05]),
                       "associations":associations})
    return output


def controls_summary(rows):
    groups=defaultdict(dict)
    for r in rows:
        m=r["record"]
        key=tuple(m[k] for k in ("dandiset","subject","session","context","region"))
        entry=(m["state"],m["variant"],m["iteration"])
        if entry in groups[key]:
            raise ValueError(f"Duplicate control {key} {entry}")
        groups[key][entry]=r
    result=[]
    for key,group in sorted(groups.items()):
        row=dict(zip(("dandiset","subject","session","context","region"),key,strict=True))
        ref="matched_controls_reference"
        baseline={state:group[state,ref,0]["fit"] for state in ("Wake","NREM")}
        row["reference_effect"]=baseline["NREM"]["d2"]-baseline["Wake"]["d2"] if all(map(valid,baseline.values())) else None
        variants = ("unit_count","rate_thinning") if ("Wake","rate_thinning",1) in group else ("unit_count",)
        for variant in variants:
            states={state:{i:r for (s,v,i),r in group.items() if s==state and v==variant}
                    for state in ("Wake","NREM")}
            fits=[]
            for i in sorted(states["Wake"].keys() & states["NREM"].keys()):
                n,w=(states[s][i]["fit"] for s in ("NREM","Wake"))
                if valid(n) and valid(w):
                    fits.append(n["d2"]-w["d2"])
            status={"effect":distribution(fits),"states":{}}
            for state, iterations in states.items():
                successful=[r["fit"] for r in iterations.values() if valid(r["fit"])]
                base=baseline[state]
                status["states"][state]={"returned_seeds":len({r["record"]["seed"] for r in iterations.values()}),
                    "actual_rate_hz":distribution([f["rate_hz"] for f in successful]),
                    "variance_change":distribution([f["count_variance"]-base["count_variance"] for f in successful]),
                    "lag1_change":distribution([f["observed_acf"][0]-base["observed_acf"][0] for f in successful])}
            row[variant]=status
        result.append(row)
    return result


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",type=Path)
    p.add_argument("output",type=Path)
    o=p.parse_args()
    dataset=json.loads(o.input.read_text())
    rows=dataset["records"]
    if isinstance(rows,dict):
        rows=[rows]
    if dataset["input"].get("mode")=="controls":
        summary={"sensitivity":controls_summary(rows)}
    elif dataset["input"]["swr_windows"]:
        summary={"windows":event_summary(rows,dataset["settings"]["seed"])}
    else:
        summary=state_summary(rows)
    o.output.parent.mkdir(parents=True,exist_ok=True)
    o.output.write_text(json.dumps({"schema_version":1,"run_id":dataset["run_id"],**summary},
                                   allow_nan=False,indent=2)+"\n")
