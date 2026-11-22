"""Calculate population and unit statistics and prepare the independent fitting stage."""
import argparse
from collections import defaultdict
from hashlib import sha256
from pathlib import Path

import numpy as np

from inputs import ROOT,FITS,DIR,export,key,reference,source,write
from metrics import ar_modes,iaaft,population_metrics,shift,unit_metrics
from theory import prepare as theoretical

STATES=("Wake","NREM","REM")
SEEDS=(4101,4102,4103)
DT=.04


def seed_for(k):
    return int.from_bytes(sha256("|".join(k).encode()).digest()[:4],"little")


def wrap(m):
    return {k:m.get(k) for k in ("dandiset","subject","session","region","context","state",
            "asset_id","unit_ids","state_path","intervals_s","version","digest")}


def key_id(m):
    return "_".join(key(m))


def fold(blocks,n):
    return [x[:len(x)//n*n].reshape(-1,n,*x.shape[1:]).sum(1)
            for x in blocks if len(x)>=n]


def selected(blocks,dt):
    # Long bins are admissible only wholly inside a single continuous segment.
    for sec in (60,120,240):
        n=round(sec/dt)
        eligible=[(i,j) for i,x in enumerate(blocks) for j in range(0,len(x)-n+1,n)]
        for i in sorted({0,(len(eligible)-1)//2,len(eligible)-1}) if eligible else []:
            b,start=eligible[i]
            yield sec,b,start,blocks[b][start:start+n]


def main_population():
    out,inventory,windows=[],[],[]
    for label,fp in FITS.items():
        fits=reference(label)
        seen=set()
        for m,blocks,provenance in source(label):
            k=key(m)
            if k not in fits or m["variant"]!="reference":
                continue
            assert k not in seen,k
            seen.add(k)
            row=fits[k]
            fit=row["fit"]
            assert len(m["unit_ids"])==len(row["record"]["unit_ids"])
            secs=sum(map(len,blocks))*DT
            assert abs(secs-fit["duration_s"])<.002,(k,secs,fit["duration_s"])
            assert np.isclose(sum(v.sum() for v in blocks),fit["rate_hz"]*secs,atol=.01),(k,fit)
            quality=bool(fit["exit_status"]==0 and fit["solver_flag"]>0 and fit["min_ar_root"]>1)
            unit_available=m["dandiset"].isdigit() or m["dandiset"].startswith("zenodo")
            meta=wrap(m)
            inv=dict(**meta,fixture=provenance,fit_path=fp,duration=secs,unit_count=len(meta["unit_ids"]),
                     species="rat" if m["dandiset"] in {"000041","000061","hc-23","hc-24"} else "mouse",
                     data_modality="sorted spike counts",spike_data_available=unit_available,
                     population_count_available=True,d2_available=quality,annotation_source=m["state_path"],
                     usable_for_neuron_timescale=False,usable_for_population_timescale=True,
                     usable_for_dimensionality=False,usable_for_MR=True,usable_for_ON_OFF=True)
            inventory.append(inv)
            if not quality:
                continue
            # Covariance estimators use a window's own mean and need valid pairs at the target lag.
            usable=[x for x in blocks if len(x)>=25]
            q=population_metrics(usable,DT)
            modes=ar_modes(fit["ar_coefficients"],DT)
            alt=row["long_reach"]
            out.append(dict(meta=meta,layer="historical_population",duration=secs,units=len(meta["unit_ids"]),
                            d2=fit["d2"],d2_AR10=alt.get("d2") if alt else None,modes=modes,
                            tau_AR_longest=modes[0]["tau"],uncertain_clock=m["dandiset"]=="hc-23",
                            valid_fit=True,**q))
            wins=[]
            for sec,i,start,x in selected(blocks,DT):
                cov=population_metrics([x],DT,False)
                if cov:
                    wins.append(dict(window_s=sec,segment_index=i,start_bin=start,**cov))
            windows.append(dict(meta=meta,windows=wins))
        assert seen==set(fits), (label,len(seen),len(fits))
        print("historical",label,len(seen),flush=True)
    return out,inventory,windows


def count_fit_meta(m,id_,variant,dt):
    return {**{k:m[k] for k in ("dandiset","subject","session","state","region","context")},
            "id":id_,"variant":variant,"bin_seconds":dt}


def focus(m):
    return m["region"] in {"PoSub","ADN","CA1","BLA","PAG","mPFC"}


def control_fits(fetch):
    table,auras=[],[]
    groups={40:[],80:[],200:[]}
    seen=set()
    for tag in ("unit_phase8","unit_phase10"):
        for m,x,provenance in source(tag,fetch):
            k=key(m)
            if m["variant"]!="matched_controls":
                continue
            assert k not in seen
            seen.add(k)
            blocks=[v.sum(1) for v in x]
            row=wrap(m)
            row["fixture"]=provenance
            sid=key_id(m)
            s=seed_for(k[:-1])
            table.append(dict(meta=row,layer="unit_matched",id=sid,duration=180,units=len(m["unit_ids"]),
                              **population_metrics(blocks,DT),**unit_metrics(x,DT,s,True)))
            groups[40].append((count_fit_meta(m,sid,"reference",DT),blocks))
            for width in (80,200):
                if not focus(m):
                    continue
                dt=width/1000
                b=fold(blocks,round(dt/DT))
                j=f"{sid}_{width}ms"
                auras.append(dict(meta=row,id=j,kind="multiscale",bin_seconds=dt,
                                  **population_metrics(b,dt)))
                groups[width].append((count_fit_meta(m,j,"reference",dt),b))
            if not focus(m):
                continue
            for seed in SEEDS:
                rng=np.random.default_rng(seed+s)
                shifted=[shift(v,rng) for v in x]
                variants=[("shift",[v.sum(1) for v in shifted],shifted),
                          ("permute",[rng.permutation(v) for v in blocks],None),
                          ("iaaft",[iaaft(v,rng) for v in blocks],None)]
                for name,b,matrix in variants:
                    j=f"{sid}_{name}_{seed}"
                    pm=population_metrics(b,DT,False)
                    if not pm:
                        continue
                    out=dict(meta=row,seed=seed,id=j,kind=name,**pm)
                    if matrix is not None:
                        c=unit_metrics(matrix,DT,s,False)
                        ac=[(v-v.mean(0)) for v in matrix]
                        auto=sum(np.square(v).sum() for v in ac)
                        full=sum(np.square(v.sum(1)).sum() for v in ac)
                        out.update(c,cross_fraction0=float(1-auto/full) if full>0 else None)
                        from metrics import lag_cov
                        auto=lag_cov(matrix,50).sum(1)
                        total=lag_cov(b,50)
                        out.update(short_cov_auto=auto.tolist(),short_cov_cross=(total-auto).tolist(),
                                   short_cov_total=total.tolist())
                    if name=="iaaft":
                        from metrics import lag_cov
                        ac,b0=lag_cov(b,20),lag_cov(blocks,20)
                        out["acf_rmse"]=float(np.sqrt(np.mean((ac/ac[0]-b0/b0[0])**2)))
                    auras.append(out)
                    groups[40].append((count_fit_meta(m,j,name,DT),b))
            # Conditional activity: exclude silent runs >=80 ms, keep separate valid active pieces.
            if m["dandiset"]=="000041":
                seg=[]
                for v in blocks:
                    starts=np.flatnonzero(np.diff(np.r_[False,v==0,False]))
                    forbidden=[(a,b) for a,b in starts.reshape(-1,2) if b-a>=2]
                    ends=[0]+[b for _,b in forbidden]
                    starts=[a for a,_ in forbidden]+[len(v)]
                    seg.extend(v[a:b] for a,b in zip(ends,starts,strict=True) if b-a>=50)
                if sum(map(len,seg))*DT>=60:
                    j=sid+"_active_only"
                    auras.append(dict(meta=row,id=j,kind="active_only",duration=sum(map(len,seg))*DT,
                                      **population_metrics(seg,DT,False)))
                    groups[40].append((count_fit_meta(m,j,"active_only",DT),seg))
        print("unit archive",tag,len(table),flush=True)
    for k in {k[:-1] for k in seen}:
        paired=[x for x in table if key(x["meta"])[:-1]==k]
        assert len({tuple(r["meta"]["unit_ids"]) for r in paired})==1,k
    return table,auras,groups


def cloud(out,fetch):
    rows,sur,groups=control_fits(fetch)
    t,ts=theoretical()
    write(out/"unit.json",rows)
    write(out/"sensitivities.json",sur)
    write(out/"theory.json",t)
    export(out,[("bins40",.04,5,10,groups[40]),("bins80",.08,5,2,groups[80]),
                ("bins200",.20,5,1,groups[200]),("null_ar1",.04,1,1,ts[1]),("null_ar5",.04,5,10,ts[5])])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=["population","fitting"])
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--fetch",action="store_true")
    a=p.parse_args()
    if a.mode=="population":
        pop,inv,windows=main_population()
        write(a.out/"population.json",pop)
        write(a.out/"inventory.json",inv)
        write(a.out/"windows.json",windows)
    else:
        cloud(a.out,a.fetch)


if __name__=="__main__":
    main()
