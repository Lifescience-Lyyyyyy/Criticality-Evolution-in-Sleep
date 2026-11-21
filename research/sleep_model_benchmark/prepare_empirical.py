"""Duration-matched 000041 mPFC reference from the already published Phase 8 count fixture."""
import argparse
import json
from hashlib import sha256
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.io import loadmat

from prepare_fits import ROOT,SOURCE,save_prepared,write,cells,stationarity

TAG="matlab-sleep-cortex-20261120"
SOURCE_COMMIT="5fe4ab5efbb4888dd63f7b2972d4195a628c83b0"
SOURCE_SHA="5276b132e0f2920012db817ba457d78718f745bffbcbcc52082fa25fd463a6e8"
BIN=.04
N=6000
MIN_UNIT,MAX_UNIT=20,36


def run(path,destination):
    assert sha256(path.read_bytes()).hexdigest()==SOURCE_SHA, "Archived count fixture SHA mismatch"
    records=loadmat(path,simplify_cells=True)["cohort"]
    reference={(r["record"]["asset_id"],r["record"]["state"]):r for r in
               json.loads((ROOT/SOURCE).read_text())["records"] if r["record"]["region"]=="mPFC"
               and r["record"]["variant"]=="reference"}
    animals=defaultdict(set)
    windows=[]
    exclusions=[]
    for item in records.values():
        m=json.loads(item["metadata"])
        if m["region"]!="mPFC" or m["variant"]!="reference":
            continue
        old=reference[m["asset_id"],m["state"]]
        digest={k.replace(":","_").replace("-","_"):v for k,v in m["digest"].items()}
        assert old["record"]==m|{"digest":digest}
        units=len(m["unit_ids"])
        if not MIN_UNIT<=units<=MAX_UNIT:
            exclusions.append(dict(subject=m["subject"],session=m["session"],state=m["state"],
                                   units=units,reason="outside frozen 20–36 unit band"))
            continue
        original=np.atleast_1d(item["counts"])
        segments=[np.atleast_1d(v) for v in original] if original.dtype==object else [original]
        intervals=np.asarray(m["intervals_s"],dtype=float).reshape(-1,2)
        assert [len(s) for s in segments]==np.atleast_1d(old["fit"]["segment_bins"]).tolist()
        fitted=old["fit"]
        assert sum(int(np.sum(s)) for s in segments)==round(fitted["rate_hz"]*fitted["duration_s"])
        assert len(segments)==len(intervals)
        candidates=[(round(a/BIN)+i, s[i:i+N])
                    for s,(a,b) in zip(segments,intervals,strict=True) for i in range(0,len(s)-N+1,N)]
        if not candidates:
            exclusions.append(dict(subject=m["subject"],session=m["session"],state=m["state"],
                                   units=units,reason="no complete contiguous 240 s window"))
            continue
        # The set removes repeats for recordings with one or two complete windows.
        selected=sorted({0,(len(candidates)-1)//2,len(candidates)-1})
        animals[m["state"]].add(m["subject"])
        for idx in selected:
            start,x=candidates[idx]
            fields={k:m[k] for k in ("dandiset","version","state","subject","session","asset_id","region","context","unit_ids")}
            fields.update(variant="reference",unit_count=units,window_index=idx,
                          candidate_count=len(candidates),limits_s=[start*BIN,(start+N)*BIN])
            windows.append((fields,cells(x)))
    assert len(animals["NREM"])==3
    assert len(animals["Wake"])==3
    assert len(animals["REM"])<=3
    assets=sorted({m["asset_id"] for m,_ in windows})
    destination.mkdir(parents=True,exist_ok=True)
    folder=destination/"empirical"
    saved=save_prepared(windows,folder/"counts.mat","states",assets,source_sha256=SOURCE_SHA)
    write(destination/"fit_plan.json",[dict(name="empirical",receipt=str(folder/"counts.json"),
          sha256=saved["sha256"],count=len(windows),analysis_spec=[BIN,5,10])])
    write(destination/"empirical_selection.json",dict(archive_commit=SOURCE_COMMIT,archive_sha256=SOURCE_SHA,
          sha256_reference_fits=sha256((ROOT/SOURCE).read_bytes()).hexdigest(),dandiset="000041",
          selection="20–36 units, first/lower-middle/last whole non-overlapping 240 s within author-scored bout",
          sample_seconds=N*BIN, animals={k:sorted(v) for k,v in animals.items()}, exclusions=exclusions,
          stationarity=[{"session":m["session"],"state":m["state"],"window_index":m["window_index"],
                         **stationarity(c[0],BIN,2100+i)} for i,(m,c) in enumerate(windows)]))
    print(f"{len(windows)} continuous windows in {len(assets)} sessions",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.archive,a.out)
