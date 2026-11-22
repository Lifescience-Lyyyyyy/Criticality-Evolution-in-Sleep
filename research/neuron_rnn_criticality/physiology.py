"""Reconstruct neuron-level recordings and audit author UP/DOWN dynamics."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

from run_jercog import COMMIT, DT, SEEDS, CONDITIONS, digest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent/"criticality_dynamics_atlas"))
from inputs import write

BIN = .01
SIZE = 5000

# Figure 7H gamma fits (shape; scale in seconds). See doi:10.7554/eLife.22425.013.
REF = {"up_mean": 1.1 * .7, "down_mean": 2.3 * .4,
       "up_cv": 1 / np.sqrt(1.1), "down_cv": 1 / np.sqrt(2.3),
       "within_up_isi_cv_E": .74, "within_up_isi_cv_I": .94}
METHOD = "Two-state Poisson HMM/Viterbi on the author 90E+10I countdet, 10 ms; primary removes interior events <50 ms."


def read_record(path):
    m=json.loads((path/"manifest.json").read_text())
    assert m["model_commit"]==COMMIT and m["integration_step_seconds"]==DT
    assert m["neuron_counts"]=={"E":4000,"I":1000}
    native={}
    counts={suffix:path/name for name in m["files"]
            for suffix in ("rasterE","rasterI","countE","countI","countdet","Ada","kicks")
            if name.endswith(suffix+".bin")}
    assert len(counts)==7 and all(digest(p)==m["files"][p.name] for p in counts.values())
    ns=round(m["duration_seconds"]/BIN)
    for label,n in m["neuron_counts"].items():
        a=np.fromfile(counts[f"raster{label}"],dtype="<i4").reshape(-1,2)
        c=np.fromfile(counts[f"count{label}"],dtype="<i4")
        assert len(c)==ns and (a>=0).all() and a[:,1].max()<n and a[:,0].max()<ns*200
        assert np.all(np.diff(a[:,0])>=0)
        # Time interval convention is [left, right), 200 native steps per 10 ms.
        np.testing.assert_array_equal(np.bincount(a[:,0]//200,minlength=ns),c)
        native[label]=a
        native["rate"+label]=c/(n*BIN)
    native["detect"]=np.fromfile(counts["countdet"],dtype="<i4")
    native["adaptation"]=np.fromfile(counts["Ada"],dtype="<i4")*.001
    native["kicks"]=np.fromfile(counts["kicks"],dtype="<f8").reshape(-1,2)
    det=np.concatenate([native["E"][(native["E"][:,1]>=991)&(native["E"][:,1]<=1080),0],
                        native["I"][(native["I"][:,1]>=99)&(native["I"][:,1]<=108),0]])
    assert ns==len(native["detect"])==len(native["adaptation"])
    np.testing.assert_array_equal(np.bincount(det//200,minlength=ns),native["detect"])
    return m,native


def hmm(count):
    """Fit a two-state Poisson HMM by EM, and decode the most likely path."""
    rate=np.array([max(.01,np.quantile(count,.15)),max(1.,np.quantile(count,.85))])
    p=np.array([[.985,.015],[.015,.985]])
    e=np.empty((len(count),2))
    a=np.empty_like(e)
    b=np.empty_like(e)
    norm=np.empty(len(count))
    previous=-np.inf
    for _ in range(45):
        loge=count[:,None]*np.log(rate)-rate
        maxima=loge.max(axis=1)
        e=np.exp(loge-maxima[:,None])
        a[0]=e[0]*.5
        norm[0]=a[0].sum()
        a[0]/=norm[0]
        for j in range(1,len(e)):
            a[j]=(a[j-1]@p)*e[j]
            norm[j]=a[j].sum()
            a[j]/=norm[j]
        likelihood=float(np.log(norm).sum()+maxima.sum())
        if abs(likelihood-previous)<.01:
            break
        previous=likelihood
        b[-1]=1.
        for j in range(len(e)-2,-1,-1):
            b[j]=p@(e[j+1]*b[j+1])/norm[j+1]
        g=a*b
        xi=(a[:-1,:,None]*p*(e[1:]*b[1:])[:,None,:]/norm[1:,None,None]).sum(0)
        p=np.clip(xi / xi.sum(1,keepdims=True),1e-8,1.-1e-8)
        rate=np.maximum((g.T@count)/g.sum(0),1e-6)
    le=count[:,None]*np.log(rate)-rate
    v=le[0].copy()
    link=np.empty((len(le),2),np.int8)
    logp=np.log(p)
    for j in range(1,len(le)):
        u=v[:,None]+logp
        link[j]=np.argmax(u,axis=0)
        v=np.max(u,axis=0)+le[j]
    labels=np.empty(len(le),np.int8)
    labels[-1]=int(np.argmax(v))
    for j in range(len(le)-1,0,-1):
        labels[j-1]=link[j,labels[j]]
    if rate[1]<rate[0]:
        labels=1-labels
        rate=rate[::-1]
    return labels.astype(bool), rate.tolist()


def periods(state):
    edge=np.r_[0,np.flatnonzero(state[1:]!=state[:-1])+1,len(state)]
    return [(bool(state[a]),int(a),int(b)) for a,b in zip(edge[:-1],edge[1:])]


def merge_short(state,minimum):
    out=state.copy()
    while True:
        q=[(b-a,a,b) for _,a,b in periods(out)[1:-1] if b-a<minimum]
        if not q:
            return out
        _,a,b=min(q)
        out[a:b]=out[a-1]


def episode_stats(state,native,burn):
    q=periods(state)
    complete=q[1:-1]
    by={s:np.asarray([b-a for k,a,b in complete if k==s])*BIN for s in (True,False)}
    res={"up_occupancy":float(np.mean(state)), "analysis_seconds":len(state)*BIN}
    for label,s in (("up",True),("down",False)):
        dur=by[s]
        res.update({f"{label}_n":len(dur),f"{label}_mean":float(dur.mean()) if len(dur) else None,
                    f"{label}_cv":float(dur.std()/dur.mean()) if len(dur)>1 else None,
                    f"{label}_durations":dur.tolist()})
    for s,field in ((False,"corr_preceding_down_up"),(True,"corr_up_following_down")):
        x,y=zip(*[(b-a,d-c) for (k,a,b),(_,c,d) in zip(complete,complete[1:]) if k==s]) if len(complete)>2 else ([],[])
        res[field]=float(np.corrcoef(x,y)[0,1]) if len(x)>3 and np.std(x)>0 and np.std(y)>0 else None
    res["aligned_eligible_up"]=len([k for k,a,b in complete if k and b-a>=80])
    for typ in ("E","I"):
        rate=native["rate"+typ][burn:burn+len(state)]
        use=[(a,b) for k,a,b in complete if k and b-a>=80]
        onset=[rate[a-10:a+35] for a,b in use if a>=10]
        offset=[rate[b-35:b+10] for a,b in use if b+10<=len(rate)]
        early=[rate[a+15:a+30].mean() for a,b in use]
        late=[rate[b-30:b-15].mean() for a,b in use]
        res[typ]={"rate_hz":float(rate.mean()),"up_hz":float(rate[state].mean()) if state.any() else None,
                  "down_hz":float(rate[~state].mean()) if not state.all() else None,
                  "up_relative_decay":float(1-np.mean(late)/np.mean(early)) if early and np.mean(early)>0 else None,
                  "onset_n":len(onset),"offset_n":len(offset),
                  "onset_hz":np.mean(onset,axis=0).tolist() if onset else [],
                  "offset_hz":np.mean(offset,axis=0).tolist() if offset else []}
    return res


def ragged(native):
    times=np.concatenate([native["E"][:,0],native["I"][:,0]])
    ids=np.r_[native["E"][:,1],native["I"][:,1]+4000]
    order=np.argsort(ids,kind="stable")
    return ids,times,np.r_[0,np.cumsum(np.bincount(ids,minlength=SIZE))],times[order]*DT


def interspike_cv(times,offset,state,burn):
    # Only consecutive spikes within the same UP episode contribute.
    episode=np.full(len(state),-1)
    for j,(s,a,b) in enumerate(periods(state)):
        if s:
            episode[a:b]=j
    metrics=[]
    for start,end in zip(offset[:-1],offset[1:]):
        t=times[start:end]
        k=np.floor(t/BIN+1e-8).astype(int)-burn
        ix=(k>=0)&(k<len(episode))
        t,k=t[ix],k[ix]
        u=episode[k]
        d=np.diff(t)[(u[1:]==u[:-1])&(u[1:]>=0)]
        metrics.append(float(d.std()/d.mean()) if len(d)>=10 and d.mean()>0 else None)
    return metrics


def compare(path,out):
    m,native=read_record(path)
    burn=round(m["burn_in_seconds"]/BIN)
    state,lam=hmm(native["detect"][burn:])
    raw=episode_stats(state,native,burn)
    state=merge_short(state,5)
    main=episode_stats(state,native,burn)
    alternate=episode_stats(merge_short(state,10),native,burn)
    ids,steps,offset,times=ragged(native)
    cv=interspike_cv(times,offset,state,burn)
    for typ,data in (("E",cv[:4000]),("I",cv[4000:])):
        q=[v for v in data if v is not None]
        main[typ]["within_up_isi_cv_median"]=float(np.median(q)) if q else None
        main[typ]["within_up_isi_cv_n"]=len(q)
        main[typ]["within_up_isi_cv_mean"]=float(np.mean(q)) if q else None
    target=out/path.name
    target.mkdir(parents=True,exist_ok=True)
    spikes=target/"simulation_record.npz"
    np.savez_compressed(spikes,neuron_ids=np.arange(SIZE,dtype=np.int32),
        neuron_types=np.repeat(["E","I"],[4000,1000]),spike_times_seconds=times,
        spike_times_offsets=offset)
    states=target/"states.npz"
    np.savez_compressed(states,up_mask=state,detector_counts=native["detect"],
        rate_E_hz=native["rateE"],rate_I_hz=native["rateI"],adaptation=native["adaptation"],
        kick_seconds=native["kicks"][:,0]*DT)
    main["gate_checks"]={
        "alternation":min(main["up_n"],main["down_n"])>=20 and .05<=main["up_occupancy"]<=.95,
        "rate_separation":all(main[t]["up_hz"] is not None and main[t]["down_hz"] is not None
                              and main[t]["up_hz"]>3*main[t]["down_hz"] for t in "EI"),
        "published_duration_range": all(main[k] is not None and .4<=main[k]/REF[k]<=2.5
                                        for k in ("up_mean","down_mean","up_cv","down_cv")) if m["condition"]=="figure7" else None,
    }
    m["network_metadata"]["actual_recorded_cells"]=SIZE
    m["provenance"].update(simulation_record_sha256=digest(spikes),states_sha256=digest(states),
                           raw_folder=path.name)
    write(target/"manifest.json",m)
    raster=[]
    for typ,n,off in (("E",4000,0),("I",1000,4000)):
        pick=np.random.default_rng(614).choice(n,65,False)
        native_spikes=native[typ]
        keep=np.isin(native_spikes[:,1],pick)&(native_spikes[:,0]*DT>=12)&(native_spikes[:,0]*DT<18)
        index={k:i+off for i,k in enumerate(pick)}
        raster.append([(int(k)*DT,index[int(n)]) for k,n in native_spikes[keep]])
    preview={"raster_E":raster[0],"raster_I":[(x,y-4000) for x,y in raster[1]],
             "time_s":(np.arange(1200,1800)*BIN).tolist(),
             "up":state[1200-burn:1800-burn].astype(int).tolist(),
             "rate_E":native["rateE"][1200:1800].tolist(),"rate_I":native["rateI"][1200:1800].tolist()}
    return dict(manifest=m,primary=main,unsmoothed=raw,min100ms=alternate,poisson_counts_per_bin=lam,
                preview=preview,rate_quantiles={t:np.quantile(np.diff(offset)[i:j]/m["duration_seconds"],
                                [.10,.25,.50,.75,.90]).tolist()
                                for t,i,j in (("E",0,4000),("I",4000,5000))})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("raw",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    runs=[compare(a.raw/f"{c}-c{k}-n{n}",a.out) for k,n in SEEDS for c in CONDITIONS]
    checked=[x["primary"] for x in runs if x["manifest"]["condition"]=="figure7"]
    qualifier=sum(x["I"]["up_relative_decay"] is not None and x["E"]["up_relative_decay"] is not None
                  and x["I"]["up_relative_decay"]>x["E"]["up_relative_decay"] for x in checked)
    pass_conditions={}
    for cond in CONDITIONS:
        pass_conditions[cond]=all(all(v is not False for v in x["primary"]["gate_checks"].values())
                                  for x in runs if x["manifest"]["condition"]==cond)
    pass_conditions["figure7"] &= qualifier>=2
    write(a.out/"physiology.json",dict(reference=REF,method=METHOD,runs=runs,passes=pass_conditions,
                                      i_decreases_more_seeds=qualifier))
    for x in runs:
        m=x["manifest"]
        s=x["primary"]
        print(m["condition"],m["simulation_seed"],*((round(s[k],3) if s[k] is not None else None) for k in ("up_mean","down_mean","up_cv","down_cv")),
              s["gate_checks"],flush=True)
    print("Physiology gates",pass_conditions,flush=True)


if __name__=="__main__":
    main()
