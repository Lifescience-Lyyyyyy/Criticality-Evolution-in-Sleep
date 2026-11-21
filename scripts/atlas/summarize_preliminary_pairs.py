"""Summarize verified 000939 pairwise fits with animal-level reweighting."""
from collections import defaultdict
from statistics import mean, median

import numpy as np

from build_atlas_data import ROOT, read, write
from build_pairwise_fixture import OUT, REFERENCE
from catalog import PAIRS

DIR=OUT+"pairwise_equal_duration/"


def paired(rows,variant,states,strict=False):
    grouped=defaultdict(dict)
    for row in rows:
        m=row["record"]
        if m["dandiset"]=="000939" and m["region"]=="PoSub" and m["context"]=="REST" and m["variant"]==variant:
            key=(*[m[k] for k in ("subject","session","asset_id")],)
            assert m["state"] not in grouped[key]
            grouped[key][m["state"]]=row
    a,b=states
    sessions=[]
    for (subject,session,asset), rs in sorted(grouped.items()):
        if not {a,b} <= rs.keys():
            assert not strict,(variant,session,tuple(rs))
            continue
        left,right=rs[a],rs[b]
        assert left["record"]["unit_ids"]==right["record"]["unit_ids"]
        if strict:
            for row in rs.values():
                f=row["fit"]
                assert f["segment_bins"]==[1500,1500,1500] and abs(f["duration_s"]-180)<1e-8
                assert f["ar_order"]==5 and f["solver_flag"]>0 and f["min_ar_root"]>1
        effect=right["fit"]["d2"]-left["fit"]["d2"]
        sessions.append(dict(subject=subject,session=session,asset_id=asset,effect=effect,
                             state_d2={a:left["fit"]["d2"],b:right["fit"]["d2"]},
                             state_seconds={s:rs[s]["fit"]["duration_s"] for s in (a,b)}))
    return sessions


def by_animal(sessions):
    groups=defaultdict(list)
    for row in sessions:
        groups[row["subject"]].append(row["effect"])
    return {s:float(mean(v)) for s,v in sorted(groups.items())}


def metrics(rows):
    a=by_animal(rows)
    xs=list(a.values())
    n=len(xs)
    return dict(n=n,sessions=len(rows),mean=mean(xs) if n else None,median=median(xs) if n else None,
                positive=sum(v>0 for v in xs),negative=sum(v<0 for v in xs),zero=xs.count(0),
                bootstrap95=np.quantile(np.random.default_rng(11101).choice(xs,size=(10000,n),replace=True).mean(axis=1),
                                        [.025,.975]).tolist() if n>=3 else None,animal_effects=a)


def main():
    fit=read(DIR+"matlab-fits.json")
    assert fit["settings"]["ar_order"]==5 and fit["settings"]["bin_seconds"]==.04 and fit["settings"]["beta"]==2
    catalog=read(DIR+"eligibility.json")
    old=read(REFERENCE)["records"]
    new=fit["records"]
    actual=defaultdict(set)
    for row in new:
        m=row["record"]
        actual[m["variant"],m["asset_id"]].add(m["state"])
    from prepare_crcns_sleep import EQUAL_DURATION_STATES
    for case in catalog:
        for v,states in EQUAL_DURATION_STATES.items():
            assert actual[v,case["asset_id"]]==(set(states) if not case["exclusions"][v] else set())
    summary={}
    for name,pair in PAIRS.items():
        legacy=paired(old,"equal_duration",pair)
        variant="equal_duration_"+name
        result=paired(new,variant,pair,strict=True)
        three=paired(new,"equal_duration_WNR",pair,strict=True)
        ref=paired(old,"reference",pair)
        before,after,primary=map(metrics,(legacy,result,ref))
        shared=set(before["animal_effects"]) & set(after["animal_effects"])
        changes={s:after["animal_effects"][s]-before["animal_effects"][s] for s in sorted(shared)}
        old_assets={s["asset_id"] for s in legacy}
        excluded=[]
        for c in catalog:
            if c["exclusions"][variant] or c["asset_id"] not in old_assets:
                missing=[s for s in pair if s not in c["candidates"]]
                legacy_reason={s:f"只有 {n} 个完整 60 秒窗" for s,n in c["candidates"].items() if n<3}
                if missing:
                    legacy_reason.update({s:"缺少有效参考状态" for s in missing})
                excluded.append(dict(subject=c["subject"],session=c["session"],asset_id=c["asset_id"],
                    legacy_exclusion=None if c["asset_id"] in old_assets else legacy_reason,
                    pair_exclusion=c["exclusions"][variant] or None))
        summary[name]=dict(states=pair,legacy=before,corrected=after,three_state_fixed=metrics(three),
                           primary_reference=primary,matched_animals=len(shared),
                           corrected_minus_legacy_same_animals=mean(changes.values()), per_animal_paired_delta=changes,
                           corrected_minus_legacy_all_animals=after["mean"]-before["mean"],
                           exclusions=excluded,paired_sessions=result)
    write(ROOT/(DIR+"summary.json"),dict(dandiset="000939",region="PoSub",context="REST",
         configuration={"bin_seconds":.04,"window_seconds":60,"windows_per_state":3,"min_seconds":180,
                        "variant_legacy":"equal_duration","bootstrap_seed":11101,"bootstrap_draws":10000,"unit":"animal",
                        "window_choice":"ascending verified complete windows; first, lower median, last"},
         provenance=read(DIR+"matlab-fits.receipt.json"), contrasts=summary))
    for k,r in summary.items():
        print(k,r["legacy"]["n"],"->",r["corrected"]["n"],r["corrected"]["mean"],r["corrected"]["bootstrap95"])


if __name__=="__main__":
    main()
