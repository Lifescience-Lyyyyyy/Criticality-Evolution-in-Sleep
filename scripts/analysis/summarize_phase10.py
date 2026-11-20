"""Summarize aligned sleep-state fits using the animal as the replicate."""
from collections import Counter, defaultdict
import argparse
import csv
import json
from pathlib import Path

import numpy as np

CONTRASTS = {"NW": ("NREM","Wake"), "RN": ("REM","NREM"), "RW": ("REM","Wake")}
STATES = ("Wake","NREM","REM")


def valid(fit):
    return bool(fit and fit.get("d2") is not None and np.isfinite(fit["d2"]) and fit["d2"]>=-1e-10
                and fit["exit_status"]==0 and (fit["solver_flag"] or 0)>0
                and (fit["min_ar_root"] or 0)>1 and (fit["innovation_variance"] or 0)>0)


def load_rows(files):
    for path in files:
        data=json.loads(path.read_text())
        rows=data["records"]
        yield from rows if isinstance(rows,list) else [rows]


def summarize(files):
    sessions=defaultdict(dict)
    state_animals=defaultdict(list)
    diagnostics=[]
    invalid=[]
    for row in load_rows(files):
        meta=row["record"]
        context=meta["context"]
        base=(meta["dandiset"],meta["region"],context,meta["subject"],meta["asset_id"],meta["session"])
        groups=[(meta["variant"],row["fit"])]
        if meta["variant"]=="reference":
            groups.append(("AR10",row["long_reach"]))
        for variant,fit in groups:
            if not valid(fit):
                invalid.append({"dataset":base[0],"region":base[1],"subject":base[3],"asset_id":base[4],
                                "context":base[2],"state":meta["state"],"variant":variant,
                                "exit_status":fit.get("exit_status") if fit else None})
                continue
            sessions[(*base,variant)][meta["state"]]=(fit["d2"],meta["unit_ids"],meta["intervals_s"])
            if variant=="reference":
                state_animals[(*base[:4],meta["state"])].append(fit["d2"])
                residual=fit["residual_acf"]
                max_res=max(abs(i) for i in residual if i is not None)
                diagnostics.append({"dataset":base[0],"region":base[1],"context":context,"subject":base[3],
                    "asset_id":base[4],"state":meta["state"],"duration_s":fit["duration_s"],"units":len(meta["unit_ids"]),
                    "rate_hz":fit["rate_hz"],"fano":fit["fano"],"zero_fraction":fit["zero_fraction"],
                    "residual_excess_kurtosis":fit["residual_excess_kurtosis"],"max_residual_acf":max_res,
                    "acf_rmse":fit["acf_rmse"],"warning":bool(abs(fit["residual_excess_kurtosis"])>1 or max_res>0.1),
                    "observed_acf":fit["observed_acf"],"model_acf":fit["model_acf"],"residual_acf":residual})
    effects=defaultdict(list)
    by_session=[]
    for (dataset,region,context,subject,asset,sess,variant), states in sessions.items():
        for contrast,(a,b) in CONTRASTS.items():
            if not {a,b}<=states.keys():
                continue
            left,right=states[a],states[b]
            if left[1]!=right[1] and variant in ("reference","AR10","edge_trim","equal_duration","strict_units"):
                raise ValueError(f"{asset} {variant}: state unit identities differ")
            value=left[0]-right[0]
            by_session.append({"dataset":dataset,"region":region,"context":context,"subject":subject,
                               "session":sess,"asset_id":asset,"variant":variant,"contrast":contrast,"effect":value})
            effects[dataset,region,context,variant,contrast,subject].append(value)
    rows=[{"dataset":d,"region":r,"context":c,"variant":v,"contrast":k,"subject":a,
           "sessions":len(nums),"effect":float(np.mean(nums))}
          for (d,r,c,v,k,a),nums in sorted(effects.items())]
    groups=defaultdict(list)
    for item in rows:
        groups[tuple(item[k] for k in ("dataset","region","context","variant","contrast"))].append(item)
    estimates=[]
    rng=np.random.default_rng(2026)
    for (d,r,c,v,k),rep in sorted(groups.items()):
        numbers=np.array([u["effect"] for u in rep])
        sample=rng.choice(numbers,(10000,len(numbers)),replace=True).mean(axis=1)
        estimates.append({"dataset":d,"region":r,"context":c,"variant":v,"contrast":k,
                          "animals":len(numbers),"sessions":sum(u["sessions"] for u in rep),
                          "mean":float(numbers.mean()),"median":float(np.median(numbers)),
                          "positive":int(np.sum(numbers>0)),"negative":int(np.sum(numbers<0)),
                          "bootstrap95":np.quantile(sample,[.025,.975]).tolist() if len(numbers)>=3 else None})
    state_values=[{"dataset":d,"region":r,"context":c,"subject":a,"state":s,
                   "d2":float(np.mean(x)),"sessions":len(x)}
                  for (d,r,c,a,s),x in sorted(state_animals.items())]
    paired=[]
    primary=defaultdict(dict)
    for (d,r,c,subject,asset,sess,v),states in sessions.items():
        if d=="000056" and c=="REST" and r in ("ADN","PoSub"):
            primary[(subject,asset,sess,v)][r]=states
    for (subject,asset,session,v),areas in primary.items():
        if set(areas)!={"ADN","PoSub"}:
            continue
        for s in areas["ADN"].keys()&areas["PoSub"].keys():
            adn,pos=areas["ADN"][s],areas["PoSub"][s]
            if adn[2]!=pos[2]:
                raise ValueError(f"Mismatched same-time regional exposure: {asset} {s} {v}")
            paired.append({"subject":subject,"asset_id":asset,"session":session,"variant":v,"state":s,
                           "adn":adn[0],"posub":pos[0],"difference_adn_posub":adn[0]-pos[0]})
    return {"cohort":estimates,"animal_effects":rows,"session_effects":by_session,
            "animal_state_d2":state_values,"simultaneous_adn_posub":paired,"diagnostics":diagnostics,"invalid":invalid}


def controls(files):
    cells=defaultdict(dict)
    invalid=Counter()
    rate_target_error=defaultdict(list)
    for row in load_rows(files):
        meta=row["record"]
        fit=row["fit"]
        key=tuple(meta[x] for x in ("dandiset","region","subject","asset_id","variant","iteration"))
        if valid(fit):
            cells[key][meta["state"]]=fit["d2"]
            if meta["variant"]=="rate_thinning":
                target=meta["target_rate_hz"]
                rate_target_error[key[:2]].append(abs(fit["rate_hz"]-target)/target)
        else:
            invalid[key[:2]]+=1
    effects=defaultdict(list)
    for (d,r,subject,asset,variant,_),states in cells.items():
        for contrast,(a,b) in CONTRASTS.items():
            if {a,b}<=states.keys():
                effects[d,r,subject,asset,variant,contrast].append(states[a]-states[b])
    return {"animal_anchors":[{"dataset":d,"region":r,"subject":s,"asset_id":a,"variant":v,"contrast":c,
                                "valid_repeats":len(x),"median":float(np.median(x)),
                                "quantiles":np.quantile(x,[.025,.975]).tolist()}
                               for (d,r,s,a,v,c),x in sorted(effects.items())],
            "invalid":[{"dataset":d,"region":r,"fits":n} for (d,r),n in sorted(invalid.items())],
            "rate_target_error":[{"dataset":d,"region":r,"fits":len(x),
                                  "median_relative_error":float(np.median(x)),
                                  "p95_relative_error":float(np.quantile(x,.95))}
                                 for (d,r),x in sorted(rate_target_error.items())]}


def write_csv(path,rows):
    if rows:
        with path.open("w",newline="") as fd:
            writer=csv.DictWriter(fd,fieldnames=list(rows[0]),lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fits",type=Path,nargs="+",required=True)
    p.add_argument("--controls",type=Path,nargs="*",default=[])
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    data=summarize(args.fits)
    data["controls"]=controls(args.controls)
    (args.out/"summary.json").write_text(json.dumps(data,allow_nan=False,indent=2)+"\n")
    for name in ("animal_effects","session_effects","animal_state_d2","simultaneous_adn_posub"):
        write_csv(args.out/(name+".csv"),data[name])
    write_csv(args.out/"control_effects.csv",data["controls"]["animal_anchors"])


if __name__=="__main__":
    main()
