"""Produce neuron dynamics, episode-conditioned surrogates and compact MATLAB inputs."""
import argparse
from hashlib import sha256
import gzip
import json
from pathlib import Path

import numpy as np

from physiology import SIZE, periods, digest, write
from inputs import export,source
from metrics import exp_fit, integrated, lag_cov, population_metrics

BIN=.04
MAX_LAG=100
SHIFT_SEEDS=(72801,72802,72803)
PARAMETER_SAMPLE_SEED=99017


def observed_counts(spikes,manifest):
    with np.load(spikes,allow_pickle=False) as n:
        times=n["spike_times_seconds"]
        off=n["spike_times_offsets"]
        ids=n["neuron_ids"]
        types=n["neuron_types"]
    seconds=manifest["duration_seconds"]
    burn=manifest["burn_in_seconds"]
    assert len(ids)==SIZE and np.array_equal(ids,np.arange(SIZE))
    assert np.count_nonzero(types=="E")==4000 and len(off)==SIZE+1 and off[-1]==len(times)
    t=np.floor(times/BIN+1e-8).astype(np.int32)
    lo,hi=round(burn/BIN),round(seconds/BIN)
    if not np.all((t>=0)&(t<hi)):
        raise ValueError("Unexpected spike time outside the simulation interval.")
    uid=np.repeat(ids,np.diff(off))
    use=t>=lo
    merged=np.bincount((t[use]-lo)*SIZE+uid[use],minlength=(hi-lo)*SIZE).reshape(-1,SIZE)
    assert merged.max()<65536
    return merged.astype(np.uint16)


def components(x,lag=MAX_LAG):
    """Common denominators, centered once per continuous segment."""
    cells=np.concatenate([lag_cov([x[:,i:i+100]],lag) for i in range(0,x.shape[1],100)],axis=1)
    total=lag_cov([x.sum(1,dtype=np.int64)],lag)
    within=cells.sum(1)
    between=total-within
    # Independent direct products of centered vectors for selected lags.
    q=x.astype(float)-x.mean(0)
    y=q.sum(1)
    for k in (0,1,min(10,lag)):
        a,b=q[:len(q)-k],q[k:]
        lhs=np.einsum("ij,ij->",a,b)/(len(q)-k)
        rhs=y[:len(q)-k]@y[k:]/(len(q)-k)
        np.testing.assert_allclose([lhs,rhs],[within[k],total[k]],atol=1e-7,rtol=2e-9)
    return cells,within,between,total


def profile(auto,cross,total):
    return dict(auto=auto.tolist(),cross=cross.tolist(),total=total.tolist(),
                cross_fraction0=float(cross[0]/total[0]) if total[0]>0 else None,
                cross_fraction_at_200ms=float(cross[5]/total[5]) if abs(total[5])>1e-8 else None)


def cell_profiles(x,cells):
    rates=x.mean(0)/BIN
    acf=np.divide(cells,cells[0],where=cells[0]>0,out=np.full_like(cells,np.nan))
    tau=np.full((len(rates),2),np.nan)
    for j in np.flatnonzero((x.sum(0)>=100)&(cells[0]>0)):
        f=exp_fit(acf[:,j],BIN,2*BIN,.8,True)
        if f.get("valid"):
            tau[j,0]=f["tau"]
        tau[j,1]=integrated(acf[:,j],BIN)[0]
    return rates,acf,tau


def pure_episodes(up):
    p=up.reshape(-1,4)
    # Mixed bins are flagged separately and remain unchanged.
    labels=np.where(np.all(p,axis=1),1,np.where(np.any(p,axis=1),2,0))
    edge=np.r_[0,1+np.flatnonzero(np.diff(labels)),len(labels)]
    return [(int(labels[a]),int(a),int(b)) for a,b in zip(edge[:-1],edge[1:])]


def shifted(x,episodes,seed):
    rng=np.random.default_rng(seed)
    y=x.copy()
    eligible=0
    err=0.
    for st,a,b in episodes:
        n=b-a
        if st==2 or n<8:
            continue
        v=x[a:b]
        steps=rng.integers(int(np.ceil(n/4)),int(np.floor(3*n/4))+1,size=x.shape[1])
        rolled=np.take_along_axis(v,(np.arange(n)[:,None]-steps)%n,axis=0)
        np.testing.assert_array_equal(v.sum(0),rolled.sum(0))
        for k in (0,1,2,min(5,n-1)):
            old=np.sum(v.astype(np.int64)*np.roll(v,-k,axis=0),axis=0)
            new=np.sum(rolled.astype(np.int64)*np.roll(rolled,-k,axis=0),axis=0)
            err=max(err,float(np.max(np.abs(new-old))))
        y[a:b]=rolled
        eligible+=n
    assert err==0.
    return y,dict(eligible_bins=eligible,total_bins=len(x),max_cyclic_autocorrelation_error=err)


def cells_of_size(n):
    pick=np.random.default_rng(PARAMETER_SAMPLE_SEED)
    return np.r_[pick.choice(4000,4*n//5,False),4000+pick.choice(1000,n//5,False)]


def slim(q):
    return {k:v for k,v in q.items() if k not in {"freq","psd","mr_coeff"}}


def entry(out, name, condition, session, category, x, **meta):
    """Register only population sums, with the exact segment boundaries implicit in block lengths."""
    blocks=[np.asarray(a,dtype=np.int64).tolist() for a in x]
    out.append(dict(id=name,session=str(session),region=meta.pop("region","cortex"),
        context=condition,state=meta.pop("state","Simulated"),variant="reference",bin_seconds=BIN,
        category=category,counts=blocks,duration_seconds=sum(map(len,blocks))*BIN,**meta))
    return name


def run_one(member,local,cases):
    m=member["manifest"]
    k=f"{m['condition']}-c{m['connectivity_seed']}-n{m['simulation_seed']}"
    a=local/k
    spikes=a/"simulation_record.npz"
    states=a/"states.npz"
    assert digest(spikes)==m["provenance"]["simulation_record_sha256"]
    assert digest(states)==m["provenance"]["states_sha256"]
    x=observed_counts(spikes,m)
    with np.load(states,allow_pickle=False) as p:
        state=p["up_mask"]
        burn=round(m["burn_in_seconds"]/.01)
        e=p["rate_E_hz"][burn:]*40
        i=p["rate_I_hz"][burn:]*10
        np.testing.assert_allclose(x[:,:4000].sum(1),e.reshape(-1,4).sum(1),rtol=1e-12)
        np.testing.assert_allclose(x[:,4000:].sum(1),i.reshape(-1,4).sum(1),rtol=1e-12)
    phases=pure_episodes(state)
    c,u,v,t=components(x)
    hz,acf,tau=cell_profiles(x,c)
    matrix=a/"neuron_dynamics_40ms.npz"
    np.savez_compressed(matrix,neuron_ids=np.arange(SIZE),counts=x,acf=acf.astype(np.float32),
                        firing_rate_hz=hz,tau_s=tau[:,0],integrated_tau_s=tau[:,1])
    ctx=m["condition"]
    meta={"model":"JercogEtAl2017","connectivity_seed":m["connectivity_seed"],"noise_seed":m["simulation_seed"],
          "parameter_case":m["parameter_config"],"provenance_sha256":digest(spikes)}
    base=k+"_original"
    q=slim(population_metrics([x.sum(1)],BIN))
    row=dict(id=base, manifest=m,population=q,decomposition=profile(u,v,t),controls=[],
             neuron_dynamics_sha256=digest(matrix),n_units=SIZE,individual=[])
    for j in range(SIZE):
        row["individual"].append([float(hz[j]),float(tau[j,0]),float(tau[j,1])])
    entry(cases,base,ctx,k,"original",[x.sum(1)],n_units=SIZE,**meta)
    for n in (100,500,1000):
        pick=cells_of_size(n)
        z=x[:,pick].sum(1)
        id_=entry(cases,k+f"_sample{n}",ctx,k,"sampling",[z],n_units=n,**meta)
        row["controls"].append(dict(id=id_,n_units=n,kind="sampling",population=slim(population_metrics([z],BIN))))
    for i,z in enumerate(np.split(x.sum(1),3)):
        id_=entry(cases,k+f"_t{i}",ctx,k,"duration",[z],n_units=SIZE,**meta)
        row["controls"].append(dict(id=id_,n_units=SIZE,kind="duration",start_s=30*i,
                                    population=slim(population_metrics([z],BIN))))
    for seed in SHIFT_SEEDS:
        y,info=shifted(x,phases,seed)
        _,au,cr,to=components(y,lag=50)
        id_=entry(cases,k+f"_shift{seed}",ctx,k,"within_episode_shift",[y.sum(1)],n_units=SIZE,
                  surrogate_seed=seed,**meta)
        row["controls"].append(dict(id=id_,n_units=SIZE,kind="shift",shift=info,
            decomposition=profile(au,cr,to),population=slim(population_metrics([y.sum(1)],BIN))))
    # The low-N comparison shares the first 60 seconds with empirical windows.
    pick=cells_of_size(10)
    y=x[:1500,pick].sum(1)
    id_=entry(cases,k+"_N10_60s",ctx,k,"empirical_match",[y],n_units=10,**meta)
    row["controls"].append(dict(id=id_,n_units=10,kind="empirical_match",
                                population=slim(population_metrics([y],BIN))))
    print(k,len(x),m["duration_seconds"],len(phases),flush=True)
    return row


def empirical(cases):
    rows=[]
    selected={}
    for m,blocks,src in source("unit_phase8"):
        if m.get("dandiset")!="000041" or m["region"]!="mPFC" or m["variant"]!="matched_controls":
            continue
        ids=m["unit_ids"]
        if len(ids)<10:
            continue
        key=(m["subject"],m["session"],m["context"])
        if key not in selected:
            seed=int.from_bytes(sha256("/".join(key).encode()).digest()[:4],"little")
            selected[key]=(tuple(ids),np.sort(np.random.default_rng(seed).choice(len(ids),10,False)))
        old,pick=selected[key]
        assert tuple(ids)==old
        x=blocks[0][:,pick]
        c,u,v,t=components(x)
        rates,_,tau=cell_profiles(x,c)
        case_id=f"emp_{m['session']}_{m['state']}"
        state=entry(cases,case_id,m["context"],m["session"],"empirical",[x.sum(1)],n_units=10,
                    region="mPFC",state=m["state"],animal=m["subject"],model="DANDI-000041",
                    provenance_sha256=src["sha256"])
        rows.append(dict(id=state,meta=m,selected_ids=[ids[j] for j in pick],source=src,
                         population=slim(population_metrics([x.sum(1)],BIN)),
                         decomposition=profile(u,v,t),
                         individual=[[float(value),float(tau[i,0]),float(tau[i,1])] for i,value in enumerate(rates)]))
    assert rows
    return rows


def epoch_qc(native,dst):
    """Test how much population memory remains after centering each episode."""
    phys=json.loads((native/"physiology.json").read_text())
    out=[]
    for run in phys["runs"]:
        m=run["manifest"]
        k=f"{m['condition']}-c{m['connectivity_seed']}-n{m['simulation_seed']}"
        root=native/k
        matrix=root/"neuron_dynamics_40ms.npz"
        with np.load(matrix,allow_pickle=False) as source:
            counts=source["counts"]
        with np.load(root/"states.npz",allow_pickle=False) as features:
            epochs=pure_episodes(features["up_mask"])
        for label,val in (("up",1),("down",0)):
            blocks=[counts[a:b] for st,a,b in epochs if st==val and b-a>=8]
            if not blocks:
                continue
            lags=min(25,max(map(len,blocks))-1)
            total=lag_cov([x.sum(1) for x in blocks],lags)
            auto=np.zeros_like(total)
            for i in range(0,SIZE,100):
                auto+=lag_cov([x[:,i:i+100] for x in blocks],lags).sum(1)
            dz=total-auto
            assert np.allclose(auto+dz,total)
            acf=total/total[0] if total[0]>0 else np.full_like(total,np.nan)
            t=integrated(acf,BIN)[0] if total[0]>0 else None
            exposure=[sum(max(0,len(x)-j) for x in blocks) for j in range(lags+1)]
            out.append(dict(condition=m["condition"],connectivity_seed=m["connectivity_seed"],
                episode=label,n_episodes=len(blocks),lag_seconds=(np.arange(lags+1)*BIN).tolist(),
                valid_pairs=exposure,acf=acf.tolist(),tau_integrated_s=t,
                covariance=profile(auto,dz,total)))
    write(dst,out)
    return out


def pack_input(bundle_path,dst):
    with gzip.open(bundle_path,"rt",encoding="utf-8") as f:
        doc=json.load(f)
    cases=doc["cases"]
    gates=doc["physiology_pass"]
    if not gates or any(passed is not True for passed in gates.values()):
        raise ValueError(f"MATLAB input requires the passed gates for the included cases: {gates}")
    assert len(set(z["id"] for z in cases))==len(cases)
    data=[({k:v for k,v in z.items() if k!="counts"},[np.asarray(c,np.float64) for c in z["counts"]])
          for z in cases]
    export(dst,[("snn_matlab",BIN,5,10,data)])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("analysis","fixture","epochs"))
    p.add_argument("source",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    if a.mode=="fixture":
        pack_input(a.source,a.out)
        return
    if a.mode=="epochs":
        epoch_qc(a.source,a.out)
        return
    phys=json.loads((a.source/"physiology.json").read_text())
    eligible=[q for q in phys["runs"] if phys["passes"][q["manifest"]["condition"]]]
    if not any(q["manifest"]["condition"]=="figure7" for q in eligible):
        raise SystemExit("Figure 7 physiology gate did not pass. Retain and investigate diagnostics.")
    cases=[]
    rs=[run_one(z,a.source,cases) for z in eligible]
    emp=empirical(cases)
    a.out.mkdir(parents=True,exist_ok=True)
    write(a.out/"analysis.json",dict(runs=rs,empirical=emp))
    gates={condition:True for condition, passed in phys["passes"].items() if passed}
    raw=json.dumps(dict(physiology_pass=gates,cases=cases),separators=(",",":")).encode()
    (a.out/"matlab_counts.json.gz").write_bytes(gzip.compress(raw,mtime=0))
    print("Official MATLAB input rows",len(cases))


if __name__=="__main__":
    main()
