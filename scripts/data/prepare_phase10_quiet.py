"""Prepare a separate, video-confirmed low-motion subset of author-scored Wake."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path
import re

import numpy as np

from audit_sleep_dandi import intersect
from extract_sleep_spikes import fetch_json, open_asset
from prepare_phase10 import inspect_nwb, prepared, save_json
from prepare_sleep_dandi import clock_bins, save_prepared

SPEED_CM_S = 2.0
VELOCITY_SECONDS = 0.5


def low_motion(asset, states, cfg):
    with open_asset(fetch_json(asset["metadata_url"])) as f:
        blue=f["acquisition/blue_led_position_pixels"]
        red=f["acquisition/red_led_position_pixels"]
        descriptions=[g.attrs["description"] for g in (blue,red)]
        scales=[float(re.search(r"(\d+(?:\.\d+)?)\s*cm per pixel",s).group(1)) for s in descriptions]
        starts=[float(g["starting_time"][()]) for g in (blue,red)]
        rates=[float(g["starting_time"].attrs["rate"]) for g in (blue,red)]
        if len(set(scales))!=1 or len(set(starts))!=1 or len(set(rates))!=1 or blue["data"].shape!=red["data"].shape:
            raise ValueError("Camera timebase or pixel calibration conflict")
        hz, start, scale = rates[0], starts[0], scales[0]
        b,r=blue["data"][:],red["data"][:]
    good=np.isfinite(b).all(axis=1)&np.isfinite(r).all(axis=1)&(b>=0).all(axis=1)&(r>=0).all(axis=1)
    midpoint=(b+r)*0.5
    span=max(2,round(VELOCITY_SECONDS*hz))
    if len(midpoint)<=span:
        raise ValueError("Camera recording shorter than the speed window")
    prefix=np.r_[0,np.cumsum(~good)]
    valid=prefix[span+1:]==prefix[:-(span+1)]
    speed=np.linalg.norm(midpoint[span:]-midpoint[:-span],axis=1)*scale*hz/span
    low=valid&(speed<=SPEED_CM_S)
    centers=start+(np.arange(len(low))+span/2)/hz
    changes=np.diff(np.r_[False,low,False].astype(int))
    camera_intervals=[(centers[a],centers[b-1]+1/hz) for a,b in
                      zip(np.flatnonzero(changes==1),np.flatnonzero(changes==-1),strict=True)]
    eligible=sorted(intersect(states,camera_intervals))
    bins=clock_bins(eligible,cfg["bin_seconds"],cfg["min_segment_seconds"])
    checksum=sha256(b.astype("<f8").tobytes()+r.astype("<f8").tobytes()).hexdigest()
    return eligible, {"video_pair_sha256":checksum,"camera_rate_hz":hz,"camera_frames":len(b),
                      "cm_per_pixel":scale,"actual_speed_window_seconds":span/hz,
                      "quiet_speed_cm_s":SPEED_CM_S,"quiet_seconds":cfg["bin_seconds"]*sum(b-a for a,b in bins),
                      "quiet_bouts":len(bins)}


def process(assets,cfg):
    attempts=[]
    for asset in sorted(assets,key=lambda r:(-r["units_by_region"]["ADN"],r["path"])):
        r,records=inspect_nwb(asset,cfg)
        adn=next((x for x in records if x.context=="REST" and x.region=="ADN"),None)
        if not adn or "Wake" not in adn.sleep_intervals:
            continue
        intervals, evidence=low_motion(asset,adn.sleep_intervals["Wake"],cfg)
        evidence.update(subject=asset["subject"],session=asset["session"],asset_id=asset["asset_id"],
                        reference_wake_seconds=r["contexts"]["REST"]["Wake"]["seconds"])
        attempts.append(evidence)
        if evidence["quiet_seconds"]>=cfg["min_state_seconds"]:
            source={**adn.provenance,
                    "state_source":adn.provenance["state_source"]+"; Wake intersected with two-LED low-motion frames"}
            match=replace(adn,sleep_intervals={"Wake":intervals,"NREM":adn.sleep_intervals["NREM"]},
                          provenance=source,context="REST_QUIET_EXPLORATORY")
            jobs=prepared([match],cfg,"states")
            for record,_ in jobs:
                record["video_quality"]=evidence
            return attempts,jobs
    return attempts,[]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--config",type=Path,default=Path("config/sleep_dandi.json"))
    args=p.parse_args()
    cfg=json.loads(args.config.read_text())
    rows=json.loads(args.audit.read_text())
    cohort={}
    for r in rows:
        if r["gate"]=="READY_3STATE" and r["units_by_region"].get("ADN",0)>=cfg["min_units"]:
            cohort.setdefault(r["subject"],[]).append(r)
    jobs, attempts=[],[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for audit, rows in pool.map(lambda grp:process(grp,cfg),cohort.values()):
            attempts.extend(audit)
            jobs.extend(rows)
    save_json(args.out/"video_audit.json",attempts)
    if jobs:
        save_prepared(jobs,args.out/"phase10_input.mat","states",sorted({m["asset_id"] for m,_ in jobs}))


if __name__=="__main__":
    main()
