"""Summarize targeted new fits, legacy model adequacy and editorial evidence classes."""
from collections import defaultdict
from math import copysign, isfinite, sqrt
from statistics import mean,median

from build_atlas_data import ROOT, build, read, write
from build_pairwise_fixture import OUT
from catalog import FIT_FILES, PAIRS, STATES, V1_FILE
from recheck_preliminary_sources import selected_reference

MATERIAL_DELTA=.25  # descriptive audit flag in bits/s, not a biological or statistical threshold


def diagnostic(row):
    f=row["fit"]
    observed=f.get("observed_acf") or []
    predicted=f.get("model_acf") or []
    exact=[(a,b) for a,b in zip(observed,predicted) if a is not None and b is not None]
    residual=[abs(v) for v in f.get("residual_acf") or [] if v is not None]
    return dict(d2=f["d2"], seconds=f["duration_s"], convergence=f["solver_flag"],ar_root=f["min_ar_root"],
                observed_acf=observed,fitted_acf=predicted,max_abs_residual_acf=max(residual,default=None),
                excess_kurtosis=f.get("residual_excess_kurtosis"),
                acf_rmse=sqrt(mean((a-b)**2 for a,b in exact)) if exact else None)


def effects(states):
    return {k:states[b]["d2"]-states[a]["d2"] for k,(a,b) in PAIRS.items() if a in states and b in states}


def impact(original,new):
    diff={k:new[k]-original[k] for k in new}
    flip=[k for k in diff if ((original[k]>0)-(original[k]<0))!=((new[k]>0)-(new[k]<0))]
    large=[k for k in diff if abs(diff[k])>=MATERIAL_DELTA]
    return dict(delta=diff,direction_changes=flip,material_delta=large,
                assessment="IMPACT_IN_SAMPLE" if flip or large else "NO_OBSERVED_IMPACT_IN_SAMPLE")


def spot_checks():
    cases=read(OUT+"sorting_checks/selection.json")["cases"]
    completed=read(OUT+"sorting_checks/matlab-fits.json")
    receipt=read(OUT+"sorting_checks/matlab-fits.receipt.json")
    source=read(OUT+"sorting_checks/source_metrics.json")
    index=selected_reference()
    group=defaultdict(dict)
    for row in completed["records"]:
        m=row["record"]
        k=m["asset_id"],m["region"],m["variant"]
        assert m["state"] not in group[k]
        group[k][m["state"]]=row
    sorting,splits=[],[]
    for c in cases:
        a,r=c["asset_id"],c["region"]
        base={x["record"]["state"]:diagnostic(x) for x in index[a,r]}
        ref=effects(base)
        metrics=next(x for x in source["NWB"]+source["brainstem"] if x.get("asset_id",x.get("asset"))==a)
        actual_method=metrics.get("quality_method",c["method"])
        c={**c,"method":actual_method,"fallback_reason":metrics.get("fallback_reason")}
        other=group[a,r,actual_method]
        orig_map={x["record"]["state"]:x["record"] for x in index[a,r]}
        if other:
            for state,row in other.items():
                m=row["record"]
                assert m["intervals_s"]==orig_map[state]["intervals_s"]
                assert set(m["unit_ids"])<=set(orig_map[state]["unit_ids"])
            assert len({tuple(v["record"]["unit_ids"]) for v in other.values()})==1
            adjusted={state:diagnostic(row) for state,row in other.items()}
            n=len(next(iter(other.values()))["record"]["unit_ids"])
            case=dict(c,alternate_units=n,reference_states=base,alternate_states=adjusted,
                      reference=ref,alternate=effects(adjusted),**impact(ref,effects(adjusted)))
        else:
            assert metrics.get("original_units")==metrics.get("subset_units")
            case=dict(c,alternate_units=metrics["subset_units"],reference_states=base,alternate_states=None,
                      reference=ref,alternate=None,assessment="SUBSET_IDENTICAL")
        sorting.append(case)
        if c.get("temporal_split"):
            first,last=(group[a,r,"within_state_"+t] for t in ("first","last"))
            for state in ("Wake","NREM"):
                u,v=first[state]["record"],last[state]["record"]
                assert u["unit_ids"]==v["unit_ids"]==orig_map[state]["unit_ids"]
                x,y=u["intervals_s"],v["intervals_s"]
                x,y=([x] if isinstance(x[0],(int,float)) else x),([y] if isinstance(y[0],(int,float)) else y)
                assert x[0][1]<=y[0][0] and u["parent_segment_s"]==v["parent_segment_s"]
            f,g=({k:diagnostic(v) for k,v in u.items()} for u in (first,last))
            assert all(abs(s["seconds"]-120)<1e-8 for d in (f,g) for s in d.values())
            splits.append(dict(c,first=f,last=g,first_effect=effects(f),last_effect=effects(g),
                               **impact(effects(f),effects(g))))
    write(ROOT/(OUT+"sorting_checks/summary.json"),dict(provenance=receipt,descriptive_material_delta=MATERIAL_DELTA,cases=sorting))
    write(ROOT/(OUT+"temporal_checks/within_state.json"),dict(provenance=receipt,cases=splits))
    return sorting


def selected_models(m):
    ds,r,context=m["dandiset"],m["region"],m["context"]
    if ds=="000041" and context=="HOME":
        return "rat|Frontal|000041|HOME"
    if ds=="000061" and r=="CA1" and context in ("PRE","POST"):
        return f"rat|CA1|000061|{context}"
    if (ds,r,context) in {("000552","CA1","HOME"),("000056","ADN","REST"),("000056","PoSub","REST"),
                         ("000939","PoSub","REST"),("zenodo:19462601","PAG","REST"),("zenodo:19462601","LDT","REST"),
                         ("hc-23","MEC II/III","OVERNIGHT")}:
        return f"{'rat' if ds=='hc-23' else 'mouse'}|{r}|{ds}|{context}"
    return None


def models():
    _,_,_,atlas=build()
    reference={c["id"]:c for cs in atlas.values() for c in cs}
    rows=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    for file in FIT_FILES:
        entries=read(file)["records"]
        for entry in entries if isinstance(entries,list) else [entries]:
            m=entry["record"]
            k=selected_models(m)
            if not k or m["variant"]!="reference":
                continue
            for spec,f in (("reference",entry["fit"]),("AR10",entry["long_reach"])):
                if f:
                    rows[k][spec][m["state"]].append((diagnostic({"fit":f}),len(m["unit_ids"])))
    result=[]
    for k,by_model in sorted(rows.items()):
        summary={}
        for spec,state_rows in by_model.items():
            for state,measure in state_rows.items():
                fits=[d for d,_ in measure]
                med=lambda key:median([x[key] for x in fits if x.get(key) is not None])
                summary.setdefault(spec,{})[state]=dict(
                    fits=len(fits),converged=sum(f["convergence"]>0 for f in fits),
                    min_ar_root=min(f["ar_root"] for f in fits),state_seconds_median=med("seconds"),
                    units_median=median(n for _,n in measure),acf_rmse_median=med("acf_rmse"),
                    residual_acf_median=med("max_abs_residual_acf"),excess_kurtosis_median=med("excess_kurtosis"),
                    residual_acf_warnings=sum(f["max_abs_residual_acf"]>.1 for f in fits),
                    kurtosis_warnings=sum(abs(f["excess_kurtosis"])>1 for f in fits))
        result.append(dict(id=k,animals=reference[k]["animals"],state_fit_summaries=summary,
                           paired={s:reference[k]["effects"][s] for s in by_model}))
    v1=read(V1_FILE)
    replay=read("results/validation/matlab-v1.json")["records"]
    prepared=read("manifests/mouse-preparation.json")["sessions"]
    unit_counts={r["session"]:r["units"] for r in prepared}
    result.append(dict(
        id="mouse|V1|PNAS-V1|PUBLISHED",animals=None,sessions=13,diagnostic_sample_session="YMV02",
        original_numerical={"fit_count":len(v1["primary"]),
                            "min_ar_root":min(r["fitted_min_root"] for r in v1["primary"]),
                            "successful_traced_solver":sum(r["solver_flag"]>0 for r in v1["optimizer_replay"]),
                            "units_median":median(unit_counts.values()),"diagnostic_units":unit_counts["YMV02"],
                            "state_seconds_median":{s:median(r["samples"]*.04 for r in v1["primary"] if r["state"]==s)
                                                    for s in ("WAKE","NREM")}},
        state_fit_summaries={s:{r["record"]["state"]:diagnostic({"fit":r[key]}) for r in replay}
                             for s,key in (("reference","fit"),("AR10","long_reach"))},
        paired={"reference":reference["mouse|V1|PNAS-V1|PUBLISHED"]["effects"]["reference"]},
        note="13 场仅有完整参考估计。AR(10)、残差和 ACF 对照限 YMV02；会话并非已知独立动物。"))
    write(ROOT/(OUT+"model_checks/adequacy.json"),dict(
        estimator=dict(beta=2,bin_seconds=.04,orders=[5,10],method="Yule-Walker",
                       categories=["NUMERICAL_VALIDITY","MODEL_ADEQUACY","BIOLOGICAL_INTERPRETATION"],
                       warnings={"abs_residual_acf":.1,"abs_kurtosis":1}),cohorts=result))
    return reference,result


def classify(c):
    species,region,ds,ctx=c["id"].split("|")
    n=c["animals"]
    if ds=="hc-23":
        return "UNRESOLVED","REM−NREM 在参考及 AR(10) 中偏负。8 场原群体静默分段的连续采集证据仍缺。"
    if region in ("LDT","DR","PRNc"):
        return "SENSITIVE","单元子集或 AR 阶数改变部分配对差值的方向。核团级样本小或共享原始动物。"
    if ds=="PNAS-V1":
        return "EXPLORATORY","13 个会话的 NW 方向一致。独立动物数和 REM 覆盖尚无可核验依据。"
    if ds=="hc-24":
        return "EXPLORATORY","仅 NREM 描述估计。缺少可配对状态。"
    if ctx in ("REST_PAPER","REST_QUIET_EXPLORATORY") or n<5:
        return "EXPLORATORY","样本较小或分期、行为场景属于补充比较。共享记录分别显示。"
    if (ds,region) in {("000041","Frontal"),("000552","CA1"),("000056","ADN"),("000939","PoSub"),("zenodo:19462601","PAG")}:
        return "ROBUST_PRELIMINARY","主要队列方向在已测试的规格中占优，动物级和会话级敏感性分别报告。"
    return "SENSITIVE","当前核团的方向一致性有限，分区或规格影响描述。"


def main():
    checks=spot_checks()
    reference,diagnostics=models()
    fits={d["id"]:d for d in diagnostics}
    evidence=[]
    for k,c in sorted(reference.items()):
        label,reason=classify(c)
        sample=[s for s in checks if k=="mouse|"+s["region"]+"|"+s["dandiset"]+"|"+s["context"]]
        evidence.append(dict(id=k,category=label,reason=reason,animals=c["animals"],sessions=c["sessions"],
             states=c["states"],paired=c["effects"]["reference"],ar10=c["effects"].get("AR10"),
             sensitivity=[{"session":s["session"],"assessment":s["assessment"],"changes":s.get("direction_changes",[])}
                          for s in sample],
             evidence="model_checks/adequacy.json" if k in fits else "independence.json",
             model_warnings="参考 AR 拟合的残差峰度与 ACF 偏离仅作模型诊断。"))
    evidence.append(dict(id="mouse|RSC|000987|REST",category="UNRESOLVED",
                         reason="现有材料只有自然慢波睡眠佐证，缺少同轴配对的静息 Wake 或 REM。",
                         animals=None,sessions=0,states=[],paired={},ar10=None,evidence="reports/22_phase10_dataset_audit.md"))
    write(ROOT/(OUT+"evidence_assessment.json"),evidence)


if __name__=="__main__":
    main()
