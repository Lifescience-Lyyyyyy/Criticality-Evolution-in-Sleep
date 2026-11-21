"""Read selected immutable NWB metadata and sorted timestamps for Phase 11.1."""
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import argparse
from pathlib import Path
from pickle import dump, load
from time import sleep

import numpy as np

from build_atlas_data import ROOT, read, write
from build_pairwise_fixture import OUT
from catalog import FIT_FILES, STATES
from recheck_sources import DANDI, extract, fetch_json, open_asset, save_prepared, matlab_segments
from build_audit_results import covered, spans, sorted_intervals, STAGING
from prepare_phase10 import adn_clock, contexts
from extract_sleep_spikes import recording_clock, bounds


def selected_reference():
    results = defaultdict(list)
    for file in FIT_FILES:
        records = read(file)["records"]
        for item in records if isinstance(records,list) else [records]:
            m=item["record"]
            if m["variant"]=="reference" and m["context"] in ("REST","HOME"):
                results[m["asset_id"],m["region"]].append(item)
    return results


def cut_spans(r):
    return [(round(a/.04),round(b/.04)) for a,b in spans(r["intervals_s"])]


def source_checks(case,f,audit,records,cfg):
    ds=case["dandiset"]
    clock=adn_clock(f) if ds=="000056" else recording_clock(f)
    original=records[0]["record"]["clock"]
    assert np.allclose(clock["seconds"],original["seconds"],rtol=0,atol=.04), (case["session"],clock,original)
    opto=[]
    state_audit={}
    if ds=="000552":
        for key in ("intervals/OptogeneticPulsePeriods","intervals/OptogeneticPulses"):
            if key in f:
                guard=cfg["opto_guard_seconds"]
                opto.extend([[float(a-guard),float(b+guard)] for a,b in bounds(f[key])])
    else:
        context=contexts(f,ds,clock,state_audit,cfg)["REST"]
        opto=state_audit.get("opto_exclusion",[])
    seen=[]
    for row in records:
        m=row["record"]
        sc=m["state"]
        published=audit["stages"][m["state_path"].lstrip("/")]
        timeline=sorted_intervals((a,b) for a,b,name in published if b>a and STAGING.get(name.lower())==sc)
        for a,b in spans(m["intervals_s"]):
            assert clock["seconds"][0]<=a<b<=clock["seconds"][1]+1e-7
            # Legacy 000552 used nearest 40 ms indices. Allow only one half-bin for its rounded endpoints.
            offset=.020000001 if ds=="000552" else 1e-8
            assert covered(a+offset,b-offset,timeline),(case["session"],sc,a,b)
            if ds!="000552":
                assert covered(a+offset,b-offset,context[sc]),(case["session"],sc,a,b)
            assert not any(min(b,v)-max(a,u)>offset for u,v in opto)
            seen.append((a,b,sc))
    seen.sort()
    assert all(x[1]<=y[0]+1e-8 for x,y in zip(seen,seen[1:]))
    wake_context="家笼说明与原分析背景核对" if ds=="000552" else "原作者休息场景并排除任务和干预"
    return {"independent_clock":clock,"published_stages":audit["stages"],"epochs":audit["epochs"],
            "within_stage_segments":len(seen),"excluded_opto_with_margin":opto,
            "outside_or_crossed_in_used_segments":0,
            "wake_context":wake_context,"source_score_checks":state_audit,
            "published_clock_observations":audit["clocks"]}


def inspect(case,records,version,cfg):
    ds=case["dandiset"]
    url=f"{DANDI}/{ds}/versions/{version}/assets/{case['asset_id']}/"
    print("READ",ds,case["session"],flush=True)
    for attempt in range(3):
        try:
            meta=fetch_json(url)
            with open_asset(meta) as f:
                rows=verify(case,records,f,meta,url,cfg)
            print("VERIFIED",ds,case["session"],len(rows[0]),flush=True)
            return rows
        except (OSError,TimeoutError):
            if attempt==2:
                raise
            sleep(3)


def verify(case,records,f,meta,url,cfg):
    ds=case["dandiset"]
    assert len(records)==3
    first=records[0]["record"]
    original=first["unit_ids"]
    assert all(x["record"]["unit_ids"]==original for x in records)
    digest=first["digest"]
    expected=set(digest.values()) if isinstance(digest,dict) else {digest}
    assert expected <= set(meta["digest"].values()), (case["session"],expected)
    # The existing Phase 11 reader only reads tables for unknown case codes.
    audit,_=extract(f,{**case,"case":"phase11_1"},[],meta,url)
    timing=source_checks(case,f,audit,records,cfg)
    u=f["units"]
    ids={str(v):i for i,v in enumerate(u["id"][:])}
    edge=np.r_[0,u["spike_times_index"][:]]
    raw=u["spike_times"][:]
    times={uid:raw[edge[ids[str(uid)]]:edge[ids[str(uid)]+1]] for uid in original}
    clock=timing["independent_clock"]["seconds"]
    refractory=cfg["refractory_seconds"]
    lratio=audit["metrics"].get("lratio",{}).get("values",{})
    enough=sum(lratio.get(str(uid)) is not None for uid in original)>=5
    method=case["method"] if case["method"]!="lratio_q50" or enough else "isi_q50"
    metrics={}
    used=sorted([s for item in records for s in cut_spans(item["record"])])
    h=sha256()
    for uid,t in times.items():
        assert np.isfinite(t).all() and len(t)>1 and np.all(np.diff(t)>0)
        # Curation at entry has checked the physical bounds. Legacy CA1 allows a 1 s tolerance.
        tolerance=1 if ds=="000552" else 1e-8
        assert clock[0]-tolerance<=t[0]<=t[-1]<clock[1]+tolerance
        if ds!="000552":
            assert np.max(np.abs(t/.00005-np.rint(t/.00005)))<1e-5
        h.update(str(uid).encode()+b"\0")
        h.update(np.asarray(t,dtype="<f8").tobytes())
        if method=="lratio_q50":
            v=lratio[str(uid)]
            if v is not None:
                metrics[uid]=v
        else:
            short,total=0,0
            for a,b in used:
                i,j=np.searchsorted(t,(a*.04,b*.04))
                gaps=np.diff(t[i:j])
                total+=len(gaps)
                short+=int(np.sum(gaps<refractory))
            if total:
                metrics[uid]=short/total
    if first.get("extracted_unit_sha256"):
        assert h.hexdigest()==first["extracted_unit_sha256"]
    pooled=np.sort(np.concatenate(list(times.values())))
    timing.update(units_checked=len(original), spike_range=[float(pooled[0]),float(pooled[-1])],
                  spike_time_sha256=h.hexdigest(), published_metadata=url, digest=meta["digest"])
    for item in records:
        f,m=item["fit"],item["record"]
        counts=matlab_segments(pooled,cut_spans(m),.04)
        assert sum(int(seg[0].sum()) for seg in counts)==round(f["duration_s"]*f["rate_hz"]), m["session"]
    assert len(metrics)>=5
    threshold=float(np.median(list(metrics.values())))
    subset=[u for u in original if u in metrics and metrics[u]<=threshold]
    assert 5<=len(subset)<=len(original)
    chosen=np.sort(np.concatenate([times[u] for u in subset]))
    result=[]
    if subset!=original:
        for item in records:
            m=item["record"]
            result.append(({**m,"unit_ids":subset,"variant":method,"threshold":threshold},
                           matlab_segments(chosen,cut_spans(m),.04)))
    quality=dict(dandiset=ds,session=case["session"],region=case["region"],subject=case["subject"],
                 asset_id=case["asset_id"],quality_method=method,requested_method=case["method"],cutoff=threshold,
                 fallback_reason="源文件缺少足量有限 L-ratio，改用预定的 1.5 ms ISI 相对排序" if method!=case["method"] else None,
                 original_units=len(original),n_with_metric=len(metrics),subset_units=len(subset),
                 units=subset,metrics={str(k):v for k,v in metrics.items()},
                 note="原始已分离单元内相对排序。并列保留。类型字段不参与。",
                 exact_original_spike_totals=True)
    return result,quality,timing


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    chosen=read(OUT+"sorting_checks/selection.json")["cases"]
    cases=[x for x in chosen if x["dandiset"].isdecimal()]
    refs=selected_reference()
    versions={**read("manifests/sleep_sources.json")["dandisets"],**read("manifests/phase10_sources.json")["sources"]}
    cfg=read("config/sleep_dandi.json")
    cache=a.out/"jobs"
    cache.mkdir(exist_ok=True)
    outputs={}
    with ProcessPoolExecutor(max_workers=3) as pool:
        futures={}
        for c in cases:
            asset=c["asset_id"]
            path=cache/(asset+".pickle")
            if path.exists():
                with path.open("rb") as stream:
                    outputs[asset]=load(stream)
            else:
                futures[pool.submit(inspect,c,refs[asset,c["region"]],versions[c["dandiset"]]["version"],cfg)]=asset
        for f in as_completed(futures):
            asset=futures[f]
            outputs[asset]=f.result()
            partial=cache/(asset+".part")
            with partial.open("wb") as stream:
                dump(outputs[asset], stream)
            partial.replace(cache/(asset+".pickle"))
    fits,quality,timing=[],[],[]
    for case in cases:
        r,q,t=outputs[case["asset_id"]]
        fits.extend(r)
        quality.append(q)
        timing.append({**{k:case[k] for k in ("dandiset","subject","session","region","asset_id")},**t})
    write(a.out/"nwb_sorting.json",quality)
    write(a.out/"nwb_timing.json",timing)
    save_prepared(fits,a.out/"nwb_counts.mat","states",sorted({m["asset_id"] for m,_ in fits}),
                  source_sha256={"sample":sha256((ROOT/(OUT+"sorting_checks/selection.json")).read_bytes()).hexdigest()})


if __name__=="__main__":
    main()
