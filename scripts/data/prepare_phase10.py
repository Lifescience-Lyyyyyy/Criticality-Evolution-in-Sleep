"""Audit published sleep scores and reuse the segment-aware CRCNS preparation."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
from urllib.error import HTTPError

import numpy as np

from audit_sleep_dandi import intersect, subtract
from crcns_adapter import SleepRecording, overlaps
from extract_sleep_spikes import (DANDI, bounds, electrode_locations, fetch_json, labeled_intervals,
                                  labels, open_asset, recording_clock, scalar_text)
from prepare_crcns_sleep import prepare_group
from prepare_sleep_dandi import clock_bins, save_prepared

VERSIONS = {"000056": "0.261106.1949", "000939": "0.261016.2310",
            "000978": "0.240511.0307", "000987": "0.261003.2342"}
REGIONS = {"anterodorsal nucleus": "ADN", "anterodorsal nucleus of thalamus": "ADN", "adn": "ADN",
           "postsubiculum": "PoSub", "dorsal presubiculum": "PoSub", "posub": "PoSub",
           "presubiculum": "PoSub", "retrosplenial area": "RSC", "rsc": "RSC",
           "field ca1": "CA1", "ca1": "CA1", "prelimbic area": "PL", "pfc": "PFC"}
STATES = {"wake": "Wake", "wakestate": "Wake", "awake": "Wake", "nrem": "NREM", "nremstate": "NREM",
          "sws": "NREM", "non-rem": "NREM", "rem": "REM", "remstate": "REM"}
ORDER = ("Wake", "NREM", "REM")


def canon(value):
    return REGIONS.get(value.strip().lower())


def span(ranges):
    return sum(b-a for a,b in ranges)


def scores(table, column):
    a,b=table["start_time"][:],table["stop_time"][:]
    raw=labels(table,column)
    if column+"_index" in table:
        edges=np.r_[0,table[column+"_index"][:]]
        raw=["+".join(raw[x:y]) for x,y in zip(edges[:-1],edges[1:],strict=True)]
    return [(float(x),float(y),c) for x,y,c in zip(a,b,raw,strict=True)]


def cleaned(rows, context, clock, audit):
    lower,upper=clock["seconds"]
    tolerance=1/clock["rate_hz"] if clock.get("rate_hz") else 1e-6
    good,bad=[],[]
    for a,b,label in rows:
        s=STATES.get(label.lower())
        if not s:
            continue
        if not np.isfinite([a,b]).all() or a>=b:
            bad.append([a if np.isfinite(a) else None,b if np.isfinite(b) else None,label,"invalid endpoints"])
        elif a<lower-tolerance or b>upper+tolerance:
            bad.append([a,b,label,"outside independent recording clock"])
        else:
            good.append((a,b,s))
    good.sort()
    conflicts=[]
    for i,(a,b,label) in enumerate(good):
        for x,y,other in good[i+1:]:
            if x>=b:
                break
            conflicts.append([max(a,x),min(b,y),label,other])
    audit.update(score_counts=dict(Counter(s for a,b,s in rows)),invalid_score_rows=bad,overlapping_score_rows=conflicts)
    blocked=[(a,b) for a,b,_,_ in conflicts]
    return {s: sorted(set(intersect(subtract([(a,b) for a,b,l in good if l==s],blocked),context))) for s in ORDER}


def rest_epochs(rows, labels_on):
    return sorted(subtract([(a,b) for a,b,l in rows if l.lower() in labels_on],
                           [(a,b) for a,b,l in rows if l.lower() not in labels_on]))


def contexts(file, ds, clock, row):
    if ds != "000056":
        return {}
    hypnogram="processing/behavior/sleep_states_buzcode"
    paper="processing/behavior/sleep_states"
    if "intervals/epochs" not in file:
        return {}
    epochs=labeled_intervals(file["intervals/epochs"],"tags")
    parts=bounds(file["intervals/sub_recordings"]).tolist()
    if overlaps(sorted(parts)) or overlaps(sorted((a,b) for a,b,_ in epochs)):
        raise ValueError("Overlapping behavioral epochs or source sub-recordings")
    rest=rest_epochs(epochs,{"sleep_pre","sleep_post"})
    home=sorted(intersect(intersect(rest,parts),[clock["seconds"]]))
    row.update(epochs=epochs,source_intervals=parts,context_support=home,state_path=hypnogram)
    if hypnogram not in file or not home:
        return {}
    result={"REST":cleaned(scores(file[hypnogram],"label"),home,clock,row)}
    row["state_description"]=file[hypnogram].attrs["description"]
    if paper in file:
        evidence={"state_path":paper}
        original=cleaned(scores(file[paper],"label"),home,clock,evidence)
        evidence["overlap_seconds"]={s:span(intersect(original[s],result["REST"][s])) for s in ORDER}
        evidence["valid_seconds"]={s:span(original[s]) for s in ORDER}
        row["paper_sleep_score"]=evidence
        result["REST_PAPER"]=original
    return result


def unit_trains(file, ds, clock, cfg, row):
    u=file["units"]
    ids=u["id"][:]
    edges=np.r_[0,u["spike_times_index"][:]].astype(np.int64)
    data=u["spike_times"]
    if len(set(ids))!=len(ids) or len(edges)!=len(ids)+1 or edges[-1]!=len(data) or np.any(np.diff(edges)<0):
        raise ValueError("Invalid NWB sorted-spike ragged index")
    electrodes=electrode_locations(file,u)
    anatomy=labels(u,"region") if "region" in u else electrodes
    groups=[""]*len(ids)
    if "electrode_group" in u:
        for i,ref in enumerate(u["electrode_group"][:]):
            if ref and "location" in file[ref]:
                groups[i]=scalar_text(file[ref]["location"])
    quality=labels(u,"cell_type") if "cell_type" in u else ["unspecified"]*len(ids)
    lab_list=u["in_lab_cell_list"][:].tolist() if "in_lab_cell_list" in u else [None]*len(ids)
    lratio=u["lratio"][:].tolist() if "lratio" in u else [None]*len(ids)
    shank=u["shank"][:].tolist() if "shank" in u else [None]*len(ids)
    cluster=u["cluster_id"][:].tolist() if "cluster_id" in u else [None]*len(ids)
    targets={"000056":{"ADN","PoSub"},"000939":{"PoSub"},"000987":{"RSC"}}[ds]
    resolution=data.attrs.get("resolution")
    if ds=="000056" and resolution!=.00005:
        raise ValueError("Unexpected unit time resolution for 20 kHz source")
    a,b=clock["seconds"]
    output,mapping=defaultdict(dict),[]
    for i,uid in enumerate(ids):
        region=canon(anatomy[i])
        checks={k for k in (canon(electrodes[i]),canon(groups[i])) if k}
        match=checks=={region} and bool(canon(electrodes[i]))
        if ds=="000056" and region=="PoSub" and row["subject"] not in {"Mouse24","Mouse25","Mouse28"}:
            match=False
        n=int(edges[i+1]-edges[i])
        item={"id":str(uid),"source_region":anatomy[i],"electrode_region":electrodes[i],"group_region":groups[i],
              "region":region,"anatomy_agrees":match,"cell_type":quality[i],"spikes":n,
              "shank":shank[i],"cluster_id":cluster[i],"in_lab_cell_list":lab_list[i],
              "lratio":lratio[i] if lratio[i] is not None and np.isfinite(lratio[i]) else None,"eligible":False}
        mapping.append(item)
        if not match or region not in targets or n<2 or quality[i].lower() in {"mua","mu","noise","bad","unsorted"}:
            continue
        t=data[edges[i]:edges[i+1]].astype(np.float64)
        delta=np.diff(t)
        finite=np.isfinite(t).all()
        grid=bool(resolution is None or (finite and np.allclose(t/resolution,np.rint(t/resolution),rtol=0,atol=1e-6)))
        outside=int(np.count_nonzero((t<a)|(t>=b)))
        duplicates=int(np.count_nonzero(delta==0))
        reversals=int(np.count_nonzero(delta<0))
        ok=bool(finite and grid and not outside and not duplicates and not reversals)
        item.update(first=float(t[0]) if finite else None,last=float(t[-1]) if finite else None,
                    on_source_grid=grid,outside_spikes=outside,duplicates=duplicates,reversals=reversals,
                    refractory=float(np.mean(delta<cfg["refractory_seconds"])) if ok else None,eligible=ok)
        if ok:
            output[region][str(uid)]=t
    row.update(units=mapping,units_total=len(ids),units_by_region={r:len(t) for r,t in output.items()},
               spike_resolution_s=resolution)
    return output


def finish(row, state_contexts, regional, cfg):
    row.update(contexts={},eligible={},gaps=[])
    if not state_contexts:
        return {**row,"gate":"BLOCKED_METADATA","reason":"Missing verified state or spontaneous context tables"},[]
    good={r:t for r,t in regional.items() if len(t)>=cfg["min_units"]}
    if not good:
        return {**row,"gate":"BLOCKED_QUALITY","reason":"No adequately mapped stable sorted population"},[]
    pooled=np.sort(np.concatenate([t for units in good.values() for t in units.values()]))
    gap_points=np.flatnonzero(np.diff(pooled)>10)
    gaps=[(float(pooled[i]),float(pooled[i+1])) for i in gap_points]
    row["gaps"]=gaps
    row["pooled_bounds"]=[float(pooled[0]),float(pooled[-1])]
    support=subtract([tuple(row["pooled_bounds"])],gaps)
    width=cfg["bin_seconds"]
    extracted={}
    for region,units in good.items():
        h=sha256()
        for uid,t in units.items():
            h.update(uid.encode()+b"\0")
            h.update(np.asarray(t,dtype="<f8").tobytes())
        extracted[region]=h.hexdigest()
    row["extracted_unit_sha256"]=extracted
    result=[]
    for context,states in state_contexts.items():
        times={s:sorted(intersect(ranges,support)) for s,ranges in states.items()}
        row["contexts"][context]={}
        for state,ranges in times.items():
            bins=clock_bins(ranges,width,cfg["min_segment_seconds"])
            row["contexts"][context][state]={"intervals":ranges,"segments":len(bins),
                                             "seconds":round(width*sum(y-x for x,y in bins),8)}
        states={s:ranges for s,ranges in times.items() if row["contexts"][context][s]["seconds"]>=cfg["min_state_seconds"]}
        if "NREM" not in states or len(states)<2:
            continue
        row["eligible"][context]=list(states)
        for r,units in good.items():
            path=row["paper_sleep_score"]["state_path"] if context=="REST_PAPER" else row["state_path"]
            digest=row["digest"].get("dandi:sha2-256") or row["digest"].get("dandi:dandi-etag") or row["digest"]["extracted_sha256"]
            result.append(SleepRecording(row["subject"],row["session"],row["dandiset"],r,tuple(units),
                tuple(units.values()),states,tuple(support),row["clock"],
                {"asset_id":row["asset_id"],"path":row["path"],"sha256":digest,"state_source":path,
                 "extracted_unit_sha256":extracted[r]}, {},context))
    primary=len(row["eligible"].get("REST",[]))
    secondary=max(map(len,row["eligible"].values()),default=0)
    if primary==3:
        row.update(gate="READY_3STATE",reason="Three aligned spontaneous states")
    elif max(primary,secondary)>=2:
        row.update(gate="READY_2STATE",reason="Two aligned spontaneous states")
    else:
        row.update(gate="PARTIAL",reason="Fewer than two spontaneous states with sufficient exposure")
    return row,result


def inspect_nwb(asset,cfg):
    row=dict(asset)
    ds=row["dandiset"]
    try:
        metadata=fetch_json(row["metadata_url"])
        row.update(digest=metadata["digest"],access=metadata.get("access",[]))
        with open_asset(metadata) as f:
            row.update(nwb_id=scalar_text(f["identifier"]),description=scalar_text(f["session_description"]),
                       subject=scalar_text(f["general/subject/subject_id"]),
                       session=scalar_text(f["general/session_id"]) if "general/session_id" in f else Path(row["path"]).stem)
            row["clock"]=recording_clock(f)
            if not row["clock"] or row["clock"]["seconds"][0]>=row["clock"]["seconds"][1]:
                return {**row,"gate":"BLOCKED_METADATA","reason":"No valid independent LFP or acquisition clock"},[]
            states=contexts(f,ds,row["clock"],row)
            trains=unit_trains(f,ds,row["clock"],cfg,row)
        return finish(row,states,trains,cfg)
    except HTTPError as error:
        if error.code in (401,403):
            return {**row,"gate":"ACCESS_REQUIRED","reason":f"Public HTTP {error.code}"},[]
        raise
    except (KeyError,ValueError) as error:
        return {**row,"gate":"BLOCKED_METADATA","reason":str(error)},[]


def catalog(ds,meta):
    info=json.loads((meta/(ds+".json")).read_text())
    version=VERSIONS[ds]
    if info["version"]!=version:
        raise ValueError(f"{ds}: immutable version mismatch")
    return sorted([{"dandiset":ds,"version":version,"path":a["path"],"asset_id":a["asset_id"],"bytes":a["size"],
                    "subject":a["path"].split("/")[0].removeprefix("sub-"),
                    "metadata_url":f"{DANDI}/{ds}/versions/{version}/assets/{a['asset_id']}/"}
                   for a in info["assets"] if a["path"].endswith(".nwb")],key=lambda a:a["path"])


def prepared(regions,cfg,mode):
    by_context=defaultdict(list)
    for r in regions:
        if mode!="controls" or r.context=="REST":
            by_context[r.context].append(r)
    result=[]
    for group in by_context.values():
        evidence={r.region:r.provenance for r in group}
        for meta,counts in prepare_group(group,cfg,mode):
            meta["extracted_unit_sha256"]=evidence[meta["region"]]["extracted_unit_sha256"]
            result.append((meta,counts))
    return result


def save_json(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(".part")
    temp.write_text(json.dumps(data,allow_nan=False,indent=2)+"\n")
    temp.replace(path)


def run_stage(candidates,worker,cfg,mode,output,workers):
    def run(asset):
        row,regions=worker(asset,cfg)
        save_json(output/"audit"/(asset["asset_id"]+".json"),row)
        print(row["dandiset"],row.get("session",row["path"]),row["gate"],row.get("units_by_region",{}),flush=True)
        return (row, prepared(regions,cfg,"states") if mode!="audit" else [],
                prepared(regions,cfg,"controls") if mode=="cohort" else [])
    if mode=="pilot":
        for asset in candidates:
            row,fits,_=run(asset)
            if any(m["context"]=="REST" for m,_ in fits):
                save_prepared(fits,output/"phase10_input.mat","states",[asset["asset_id"]])
                return
        print("No eligible pilot")
        return
    audit,fits,controls=[],[],[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for row,primary,sensitivity in pool.map(run,candidates):
            audit.append(row)
            fits.extend(primary)
            controls.extend(sensitivity)
    save_json(output/"audit.json",audit)
    if fits:
        save_prepared(fits,output/"cohort"/"phase10_input.mat","states",sorted({m["asset_id"] for m,_ in fits}))
    if controls:
        strongest={}
        for m,_ in controls:
            key=m["subject"],m["region"]
            score=(len(m["unit_ids"]),m["session"])
            strongest[key]=max(strongest.get(key,score),score)
        chosen=[(m,c) for m,c in controls if m["session"]==strongest[(m["subject"],m["region"])][1]]
        save_prepared(chosen,output/"controls"/"phase10_input.mat","controls",sorted({m["asset_id"] for m,_ in chosen}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset",choices=tuple(VERSIONS),required=True)
    parser.add_argument("--meta",type=Path,required=True)
    parser.add_argument("--config",type=Path,default=Path("config/sleep_dandi.json"))
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--mode",choices=("pilot","audit","cohort"),required=True)
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    cfg=json.loads(args.config.read_text())
    run_stage(catalog(args.dataset,args.meta),inspect_nwb,cfg,args.mode,args.output,args.workers)


if __name__=="__main__":
    main()
