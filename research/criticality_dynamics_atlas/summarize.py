"""Animal-level summaries of matched observations and estimator controls."""
from collections import defaultdict
from hashlib import sha256
import argparse
import csv
import gzip
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.stats import rankdata

from inputs import key,write
from metrics import ar_modes

STATES=("Wake","NREM","REM")
PAIRS=(("NREM","Wake"),("REM","NREM"),("REM","Wake"))
MEASURES={
    "d2":None,"d2_AR10":None,"tau_AR_longest":"A","tau_population":"B","tau_integrated":"B",
    "neuron_tau_median":"B","neuron_tau_same_units":"B","m_MR":"B","m_short":"B","m_half":"B",
    "tau_MR":"B","pr_raw_n":"C","pr_corr_n":"C","pr_fixed10":"C","top_fraction":"C",
    "cross_fraction0":"C","cross_fraction_integrated":"C","pair_r":"C","coupling":"C","fano":"B",
    "neuron_fano_median":"C","silent_fraction":"B","off_duration":"B","off_hz":"B",
    "slope":"B","knee_hz":"B","persist1s":"B","persist4s":"B","variance":"B",
    "neuron_rate_tau_rho":"B","low_high":"B",
    "mean_rate":"B","rate_per_unit":"B"}
PRIMARY_ASSOC={
    ("d2","tau_population"),("d2","tau_integrated"),("d2","tau_AR_longest"),
    ("d2","neuron_tau_median"),("d2","m_MR"),("d2","pr_corr_n"),("d2","cross_fraction0"),
    ("d2","fano"),("d2","silent_fraction"),("d2","off_duration"),("d2","slope"),
    ("d2","knee_hz"),("d2","persist1s"),("d2","persist4s"),("tau_population","m_MR"),
    ("tau_integrated","pr_corr_n"),("m_MR","cross_fraction0"),
    ("neuron_tau_median","m_MR"),("tau_population","knee_hz"),("m_MR","low_high")}
PIN="f26c782dca3b1c9ab18cc92c6449591c44e0d507"


def read(p):
    return json.loads(p.read_bytes())


def fit_ok(x,p):
    return (bool(x) and x["d2"] is not None and x["d2"]>=0 and x["exit_status"]==0
            and x["min_root"] is not None and x["min_root"]>1
            and x["innovation_variance"] is not None and x["innovation_variance"]>0 and (p==1 or x["flag"] is not None and x["flag"]>0))


def read_fits(p,details):
    raw=p.read_bytes()
    assert sha256(raw).hexdigest()==details["digest"].removeprefix("sha256:")
    fits,records={},[]
    with ZipFile(p) as z:
        assert not z.testzip()
        plan={x["name"]:x for x in json.loads(z.read("plan.json"))}
        for name,entry in plan.items():
            result=json.loads(z.read(name+".fits.json"))
            assert result["plan"]==entry and result["toolkit_commit"]==PIN
            assert str(result["run_id"])==str(details["id"]) and result["commit"]==details["commit"]
            assert "R2024b" in result["matlab"] and len(result["rows"])==entry["records"]
            width,p1,p2=entry["spec"]
            for r in result["rows"]:
                assert r["id"] not in fits
                a,b=r["fit"],r["alternate"]
                val,alt=fit_ok(a,p1),fit_ok(b,p2)
                fits[r["id"]]=dict(d2=a["d2"] if val else None,d2_AR10=b["d2"] if alt else None,
                    valid_fit=val,fit=a,alternate=b or None,ar_order=p1,alt_order=p2,
                    modes=ar_modes(a["ar_coefficients"],width) if val else [],
                    tau_AR_longest=modes[0]["tau"] if (modes:=ar_modes(a["ar_coefficients"],width) if val else []) else None)
            records.append(dict(spec=entry["spec"],name=name,digest=entry["input_sha256"],n=len(result["rows"])))
    return fits,dict(details,batches=records,matlab=result["matlab"])


def valid_values(x):
    return [v for v in x if v is not None and np.isfinite(v)]


def avg(x):
    q=valid_values(x)
    return float(np.mean(q)) if q else None


def ci(values,fn=np.mean):
    x=np.asarray(values)
    if len(x)<3:
        return None
    rng=np.random.default_rng(44211)
    return np.quantile(fn(x[rng.integers(0,len(x),(1800,len(x)))],axis=1),[.025,.975]).tolist()


def group_id(r):
    m=r["meta"]
    return (r["layer"],m["dandiset"],m["region"],m["context"])


def groups(rows):
    clusters=defaultdict(list)
    for r in rows:
        r["rate_per_unit"]=r.get("mean_rate",0)/r["units"] if r.get("mean_rate") is not None else None
        r["pr_fixed10"]=r.get("fixed_n",{}).get("10",[None,None])[1]
        cells=[(u[0],u[2]) for u in r.get("neurons",[]) if u[2] is not None]
        if len(cells)>=5:
            rate,tau=rankdata(cells,axis=0).T
            rho=r_value(rate,tau)
            r["neuron_rate_tau_rho"]=float(rho) if np.isfinite(rho) else None
        clusters[group_id(r)].append(r)
    return clusters


def state_pairs(rows):
    by_session=defaultdict(dict)
    for r in rows:
        m=r["meta"]
        by_session[m["subject"],m["session"]][m["state"]]=r
    deltas=[]
    states=defaultdict(list)
    for (animal,session),st in by_session.items():
        for n,r in st.items():
            states[animal,n].append(r)
        for hi,lo in PAIRS:
            if hi not in st or lo not in st:
                continue
            a,b=st[hi],st[lo]
            row=dict(animal=animal,session=session,contrast=hi+"_"+lo)
            for k in MEASURES:
                x,y=a.get(k),b.get(k)
                row[k]=x-y if x is not None and y is not None else None
            row["paired_rate_log"]=float(np.log(a["rate_per_unit"]/b["rate_per_unit"])) if a["rate_per_unit"] and b["rate_per_unit"] else None
            row["paired_duration_log"]=float(np.log(a["duration"]/b["duration"]))
            if a.get("neurons") and b.get("neurons"):
                aa={u:x for u,x in zip(a["meta"]["unit_ids"],a["neurons"],strict=True)}
                bb={u:x for u,x in zip(b["meta"]["unit_ids"],b["neurons"],strict=True)}
                assert aa.keys()==bb.keys()
                dd=[u[2]-bb[i][2] for i,u in aa.items() if u[2] is not None and bb[i][2] is not None]
                row["paired_cells"]=len(dd)
                row["neuron_tau_same_units"]=float(np.median(dd)) if dd else None
            deltas.append(row)
    return deltas,states


def aggregate_deltas(deltas,measures):
    groups=defaultdict(list)
    for d in deltas:
        groups[d["contrast"],d["animal"]].append(d)
    return [dict(contrast=con,animal=animal,session_n=len(ds),
                 **{k:avg(r.get(k) for r in ds) for k in measures})
            for (con,animal),ds in sorted(groups.items())]


def r_value(u,v):
    u=u-u.mean(axis=-1,keepdims=True)
    v=v-v.mean(axis=-1,keepdims=True)
    den=np.sqrt((u*u).sum(-1)*(v*v).sum(-1))
    return np.divide((u*v).sum(-1),den,where=den>0,out=np.full(den.shape,np.nan))


def corr(xs,ys,labels=None,adjust=None):
    a=np.asarray(xs,dtype=float)
    b=np.asarray(ys,dtype=float)
    good=np.isfinite(a)&np.isfinite(b)
    if adjust is not None:
        z=np.asarray(adjust)[good]
        if not np.isfinite(z).all():
            return None
    a,b=a[good],b[good]
    ids=np.arange(len(a)) if labels is None else np.asarray(labels)[good]
    unique=np.unique(ids)
    if len(unique)<5 or a.std()==0 or b.std()==0:
        return None
    ar,br=rankdata(a),rankdata(b)
    if adjust is not None:
        z=z[:,np.ptp(z,axis=0)>1e-8]
        if len(unique)<max(12,5*(1+z.shape[1])) or not z.shape[1]:
            return None
        z=np.c_[np.ones(len(z)),rankdata(z,axis=0)]
        ar=ar-z@np.linalg.lstsq(z,ar,rcond=None)[0]
        br=br-z@np.linalg.lstsq(z,br,rcond=None)[0]
    rho=float(r_value(ar,br))
    cluster=len(unique)!=len(a)
    rng=np.random.default_rng(991)
    if cluster:
        by=[np.flatnonzero(ids==v) for v in unique]
        boot=[np.concatenate([by[i] for i in rng.integers(0,len(unique),len(unique))]) for _ in range(600)]
        sim=np.asarray([r_value(rankdata(a[z]),rankdata(b[z])) for z in boot])
        pv=None
    else:
        boot=rng.integers(0,len(a),size=(1200,len(a)))
        sim=r_value(ar[boot],br[boot]) if adjust is not None else r_value(
            rankdata(a[boot],axis=1),rankdata(b[boot],axis=1))
        # Standard rank-permutation for unadjusted effects; adjusted rank residuals reported descriptively.
        pv=float((1+np.count_nonzero(np.abs(r_value(ar,np.asarray([rng.permutation(br) for _ in range(999)])))>=abs(rho)))/1000) if adjust is None else None
    set_=sim[np.isfinite(sim)]
    return dict(rho=rho,n_animals=len(unique),n_points=len(a),
                ci=np.quantile(set_,[.025,.975]).tolist() if len(set_) else None,p=pv)


def association(state_rows,delta_rows):
    out=[]
    for x,y in sorted(PRIMARY_ASSOC):
        for state in STATES:
            q=[r for r in state_rows if r["state"]==state]
            stats=corr([r[x] for r in q],[r[y] for r in q])
            if stats:
                out.append(dict(x=x,y=y,level="state:"+state,type=MEASURES[y] if x=="d2" else "shared_data",**stats))
        for c in ["NREM_Wake","REM_NREM","REM_Wake"]:
            q=[r for r in delta_rows if r["contrast"]==c]
            for mode in ("delta","adjusted_delta"):
                z=[[r["paired_rate_log"],r["paired_duration_log"]] for r in q]
                z=np.asarray(z,dtype=float) if mode=="adjusted_delta" and q else None
                stats=corr([r[x] for r in q],[r[y] for r in q],adjust=z)
                if stats:
                    out.append(dict(x=x,y=y,level=mode+":"+c,type=MEASURES[y] if x=="d2" else "shared_data",**stats))
        both=[r for r in state_rows if r[x] is not None and r[y] is not None]
        animal=[r["animal"] for r in both]
        if len(set(animal))>=5:
            a=np.asarray([r[x] for r in both]); b=np.asarray([r[y] for r in both])
            levels=[("raw",a.copy(),b.copy()),("within_animal",a.copy(),b.copy())]
            for v in set(animal):
                idx=np.asarray(animal)==v
                if idx.sum()<2:
                    levels[1][1][idx]=np.nan
                    levels[1][2][idx]=np.nan
                else:
                    levels[1][1][idx]-=a[idx].mean()
                    levels[1][2][idx]-=b[idx].mean()
            for mode,u,w in levels:
                stats=corr(u,w,animal)
                if stats:
                    out.append(dict(x=x,y=y,level=mode,type=MEASURES[y] if x=="d2" else "shared_data",**stats))
    return out


def inference(rows):
    effects,states_out,assoc=[],[],[]
    for k,rs in groups(rows).items():
        d,ss=state_pairs(rs)
        measure=[*MEASURES,"paired_cells","paired_rate_log","paired_duration_log"]
        q=aggregate_deltas(d,measure)
        state_rows=[]
        for (animal,state),r in ss.items():
            state_rows.append(dict(group=k,animal=animal,state=state,session_n=len(r),
                 **{z:avg(x.get(z) for x in r) for z in MEASURES}))
        for contrast in ["NREM_Wake","REM_NREM","REM_Wake"]:
            ds=[r for r in q if r["contrast"]==contrast]
            for z in MEASURES:
                vals={r["animal"]:r[z] for r in ds if r[z] is not None}
                if vals:
                    xs=list(vals.values())
                    effects.append(dict(group=k,contrast=contrast,metric=z,animals=vals,
                        n=len(xs),median=float(np.median(xs)),mean=float(np.mean(xs)),ci=ci(xs),
                        positive=int(np.sum(np.asarray(xs)>0))))
        assoc.extend(dict(group=k,**r) for r in association(state_rows,q))
        states_out.extend(state_rows)
    family=sorted((r for r in assoc if r["p"] is not None),key=lambda r:r["p"])
    bound=1.
    for j in range(len(family)-1,-1,-1):
        bound=min(bound,len(family)*family[j]["p"]/(j+1))
        family[j]["q_BH"]=bound
    return effects,states_out,assoc


def simultaneous(rows):
    groups=defaultdict(list)
    for r in rows:
        m=r["meta"]
        if r["layer"] in {"unit_matched","NWB_matched"}:
            groups[(r["layer"],m["dandiset"],m["subject"],m["asset_id"],m["context"],m["state"])].append(r)
    out=[]
    for k,rows in groups.items():
        for i,a in enumerate(rows):
            for b in rows[i+1:]:
                u,v=sorted((a,b),key=lambda x:x["meta"]["region"])
                if u["meta"]["intervals_s"]!=v["meta"]["intervals_s"]:
                    continue
                diffs={field:u[field]-v[field] for field in MEASURES
                       if u.get(field) is not None and v.get(field) is not None}
                out.append(dict(key=k,session=u["meta"]["session"],region_a=u["meta"]["region"],
                                region_b=v["meta"]["region"],**diffs))
    return out


def multiscale_profile(rows):
    """Summarize the shape of the paired lag profile without fitting to d₂."""
    bins=(2,5,10,25,50,100)
    measures=["d2",*(f"lag_{k}" for k in bins)]
    pairs=defaultdict(list)
    for r in rows:
        m=r["meta"]
        if r.get("acf") and r.get("d2") is not None:
            values={f"lag_{k}":r["acf"][k] for k in bins}
            pairs[group_id(r)].append(dict(r,**values))
    out=[]
    for group,rs in sorted(pairs.items()):
        by_session=defaultdict(dict)
        for r in rs:
            m=r["meta"]
            by_session[m["subject"],m["session"]][m["state"]]=r
        ds=[]
        for (animal,_),st in by_session.items():
            for a,b in PAIRS:
                if a in st and b in st and all(st[a][k] is not None and st[b][k] is not None for k in measures):
                    ds.append(dict(animal=animal,contrast=a+"_"+b,**{k:st[a][k]-st[b][k] for k in measures}))
        for con in sorted({d["contrast"] for d in ds}):
            ar=aggregate_deltas([d for d in ds if d["contrast"]==con],measures)
            if len(ar)<10:
                continue
            x=np.asarray([[r[k] for k in measures[1:]] for r in ar])
            sd=x.std(0)
            if not np.all(sd>0):
                continue
            u,s,v=np.linalg.svd((x-x.mean(0))/sd,full_matrices=False)
            sign=np.sign(v[0,np.argmax(np.abs(v[0]))])
            result=corr([r["d2"] for r in ar],u[:,0]*s[0]*sign)
            if result:
                result.pop("p")
                out.append(dict(group=group,contrast=con,lag_s=[k*.04 for k in bins],
                    loadings=(sign*v[0]).tolist(),variance_explained=float(s[0]**2/(s@s)),**result))
    return out


def region_interactions(sim):
    by_rec=defaultdict(dict)
    for r in sim:
        layer,dataset,animal,asset,ctx,state=r["key"]
        g=(layer,dataset,r["region_a"],r["region_b"],ctx)
        by_rec[g,animal,asset][state]=r
    effects=defaultdict(list)
    for (g,animal,_),r in by_rec.items():
        for high,low in PAIRS:
            if high in r and low in r:
                for k in MEASURES:
                    if k in r[high] and k in r[low]:
                        effects[g,high+"_"+low,k,animal].append(r[high][k]-r[low][k])
    grouped=defaultdict(dict)
    for (g,c,k,animal),v in effects.items():
        grouped[g,c,k][animal]=avg(v)
    return [dict(group=g,contrast=c,metric=k,animals=a,n=len(a),mean=avg(a.values()),ci=ci(list(a.values())))
            for (g,c,k),a in sorted(grouped.items())]


def surrogate_effects(rows,sens):
    anchors={r["id"]:r for r in rows if r.get("id")}
    by=defaultdict(list)
    for r in sens:
        if r["kind"] not in {"iaaft","shift","permute"}:
            continue
        name=r["id"].rsplit("_",2)[0]
        original=anchors[name]
        states=r["meta"]
        if r["kind"]=="iaaft" and r["acf_rmse"]>.025:
            continue
        t=group_id(original)
        q={v:r[v]-original[v] for v in ("d2","tau_integrated","fano","silent_fraction","off_duration",
                "pr_raw_n","cross_fraction0","burst_p95") if r.get(v) is not None and original.get(v) is not None}
        by[t,states["subject"],states["session"],states["state"],r["kind"]].append(q)
    rows=[]
    for (g,a,s,state,k),rep in by.items():
        rows.append(dict(group=g,animal=a,session=s,state=state,kind=k,n=len(rep),
                    **{x:avg(r.get(x) for r in rep) for x in set().union(*[set(r) for r in rep])}))
    return rows


def compress(p,data):
    raw=json.dumps(data,ensure_ascii=False,allow_nan=False,separators=(",",":")).encode()
    p.write_bytes(gzip.compress(raw,mtime=0))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache",type=Path)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--artifact",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--native-fits",type=Path)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    fits,receipt=read_fits(a.artifact,read(a.receipt))
    pop=read(a.cache/"pop/population.json")
    with ZipFile(a.artifact) as artifact:
        names=("unit.json","sensitivities.json","theory.json")
        have=all(n in artifact.namelist() for n in names)
        unit,sens,theory=(json.loads(artifact.read(n)) if have else read(a.cache/"matched"/n)
                          for n in names)
        receipt["metric_origin"]="same_runner" if have else "independent_replay"
        if not have:
            local={x["name"]:x for x in read(a.cache/"matched/plan.json")}
            assert all(local[b["name"]]["records"]==b["n"] for b in receipt["batches"])

    native=read(a.cache/"native.json") if (a.cache/"native.json").exists() else []
    if a.native_fits:
        fitted = {key(r["meta"]): r for r in read(a.native_fits)["rows"]
                  if r["layer"] == "PAG_1ms" and r["width_ms"] == 40}
        for r in native:
            if r["layer"] != "PAG_1ms":
                continue
            f = fitted.pop(key(r["meta"]))
            np.testing.assert_allclose(f["meta"]["intervals_s"], r["meta"]["intervals_s"],
                                       rtol=0, atol=6e-9)
            assert f["fit"]["duration"] == r["duration"] and f["units"] == r["units"]
            assert f["source_sha256"] == r["raw_receipt"]["fixture_sha256"]
            modes = ar_modes(f["fit"]["ar_coefficients"], .04) if f["d2"] is not None else []
            r.update(d2=f["d2"], d2_AR10=f["d2_alt"], modes=modes,
                     tau_AR_longest=modes[0]["tau"] if modes else None,
                     fit=f["fit"], alternate=f["alternate"])
        assert not fitted
    for group in (unit,sens,theory):
        for r in group:
            r.update(fits.pop(r["id"]))
    assert not fits
    for t in [*pop,*unit,*native]:
        if "valid_fit" not in t:
            t["valid_fit"]=t.get("d2") is not None
    rows=pop+unit+native
    e,s,ass=inference(rows)
    inventory=read(a.cache/"pop/inventory.json")
    usable={key(r["meta"]) for r in unit+native}
    for item in inventory:
        item["usable_for_neuron_timescale"]=item["usable_for_dimensionality"]=key(item) in usable
    for name,doc in (("inventory.json.gz",inventory),("observations.json.gz",rows),
                     ("surrogates.json.gz",sens),("continuous_windows.json.gz",read(a.cache/"pop/windows.json"))):
        compress(a.out/name,doc)
    sim=simultaneous(rows)
    for name,doc in (("paired.json",e),("animal_states.json",s),("associations.json",ass),("null_models.json",theory),
                     ("simultaneous.json",sim),("region_interactions.json",region_interactions(sim)),
                     ("acf_profile_associations.json",multiscale_profile(rows)),("surrogate_effects.json",surrogate_effects(unit,sens)),
                     ("matlab_receipt.json",receipt)):
        write(a.out/name,doc)
    fields=["layer","dandiset","subject","session","region","context","state","units","duration",
            *MEASURES,"valid_fit","neuron_tau_fraction","block_mean_cv"]
    with (a.out/"atlas.csv").open("w",newline="") as dest:
        wr=csv.DictWriter(dest,fieldnames=fields,extrasaction="ignore",lineterminator="\n")
        wr.writeheader()
        for r in rows:
            wr.writerow(r|r["meta"])
    print(f"Animal summaries {len(e)}, association rows {len(ass)}, observations {len(rows)}.")


if __name__=="__main__":
    main()
