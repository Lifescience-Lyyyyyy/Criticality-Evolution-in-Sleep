"""Audit the author recordings before staging their matched-60-s Sooter fits."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.signal import welch
import yaml

from analyse_snn import BIN, cell_profiles, components, entry, profile, slim
from metrics import population_metrics
from physiology import digest, episode_stats, hmm, merge_short, write
from run_bugnon import CASES, REPLICATE

STABLE={"wake":"02_noise","sleep":"03_sleep"}


def metadata(root,kind):
    for file in sorted(root.rglob(kind+"*.yml")):
        m=yaml.full_load(file.read_text())
        if isinstance(m,dict) and m.get("type")==kind and str(m.get("layer_name","")).startswith("Vp"):
            yield m


def source_files(root,recorder):
    for name in recorder["filenames"]:
        files=list(root.rglob(Path(name).name))
        if len(files)!=1:
            raise ValueError(f"Expected one recorder file for {name}, found {files}")
        yield files[0]


def events(root,r):
    cols=r["colnames"]
    assert cols==["gid","time"]
    rows=[np.loadtxt(p,ndmin=2) for p in source_files(root,r) if p.stat().st_size]
    return np.concatenate([a for a in rows if a.size],axis=0) if any(a.size for a in rows) else np.empty((0,2))


def isi_cv(times, cells, n):
    """CV of ISIs, with at least ten spikes in the stable session."""
    same=cells[1:]==cells[:-1]
    intervals=np.diff(times)[same]
    who=cells[1:][same]
    counts=np.bincount(who,minlength=n)
    mu=np.bincount(who,weights=intervals,minlength=n)/np.maximum(counts,1)
    second=np.bincount(who,weights=intervals**2,minlength=n)/np.maximum(counts,1)
    return np.divide(np.sqrt(np.maximum(second-mu**2,0)),mu,
                     out=np.full(n,np.nan),where=(counts>=9)&(mu>0))


def cortical_vm(root,t0,t1):
    # Retain each author's cortical recorder as a separate voltage series.
    out={}
    for r in metadata(root,"multimeter"):
        cols=r["colnames"]
        if "V_m" not in cols:
            continue
        sums=np.zeros(round((t1-t0)/10))
        counts=np.zeros_like(sums)
        xv=cols.index("V_m")
        for file in source_files(root,r):
            with file.open() as stream:
                while True:
                    a=np.loadtxt(stream,usecols=(1,xv),max_rows=400000,ndmin=2)
                    if not len(a):
                        break
                    use=(a[:,0]>=t0)&(a[:,0]<t1)
                    j=np.floor((a[use,0]-t0)/10).astype(int)
                    sums+=np.bincount(j,weights=a[use,1],minlength=len(sums))
                    counts+=np.bincount(j,minlength=len(counts))
        ok=counts>0
        if ok.mean()<.98:
            continue
        v=np.interp(np.arange(len(sums)),np.flatnonzero(ok),sums[ok]/counts[ok])
        f,p=welch(v,fs=100,nperseg=800)
        band=(f>=.5)&(f<=2)
        out[r["label"]]={"band_power_mv2":float(p[band].sum()*(f[1]-f[0])),
                         "sample_fraction":float(ok.mean()),"f":f.tolist(),"psd":p.tolist()}
    return out


def audit(root,authors,condition):
    rec=json.loads((root/"run_receipt.json").read_text())
    assert rec["model_commit"]==REPLICATE and rec["condition"]==condition and rec["resolution_ms"]==.2
    assert digest(authors/rec["parameter_source"])==rec["parameter_sha256"]
    times_path,=root.rglob("session_times.yml")
    sessions=yaml.full_load(times_path.read_text())
    start,stop=sessions[STABLE[condition]]
    assert stop-start==(60000 if condition=="wake" else 120000)
    reference=authors/"simulations"/CASES[condition]/"analysis"
    with (reference/"df_rate_stats.csv").open() as f:
        refs={r["metadata_filename"]:r for r in csv.DictReader(f) if r["session"]==STABLE[condition]}
    with (reference/"df_cv_stats.csv").open() as f:
        cv_refs={r["metadata_filename"]:r for r in csv.DictReader(f)
                 if r["session"]==STABLE[condition] and r["min_n_spikes"]=="10"
                 and r["duration"]=="0" and r["t_min"]=="0"}
    stats=[]
    corpus=[]
    for r in metadata(root,"spike_detector"):
        a=events(root,r)
        ids=np.asarray(sorted(r["gids"]),dtype=np.int32)
        which=np.searchsorted(ids,a[:,0]).astype(np.int32)
        assert np.all(which<len(ids)) and np.array_equal(a[:,0],ids[which])
        use=(a[:,1]>=start)&(a[:,1]<stop)
        rates=np.bincount(which[use],minlength=len(ids))*1000/(stop-start)
        ref=refs[r["label"]]
        assert int(float(ref["count"]))==len(ids)
        typ=r["population_name"].rsplit("_",1)[-1]
        stats.append(dict(label=r["label"],layer=r["layer_name"],type=typ,n=len(ids),
            mean_rate=float(rates.mean()),reference_mean=float(ref["mean"]),
            quantiles=np.quantile(rates,[.25,.5,.75]).tolist(),
            reference_quantiles=[float(ref[c]) for c in ("25%","50%","75%")],
            zero_rate_fraction=float(np.mean(rates==0))))
        corpus.append((ids,typ,r["layer_name"],a[use,1]-start,which[use]))
    assert stats and {x["type"] for x in stats}=={"PC","BC","MC"}
    id_all=np.concatenate([x[0] for x in corpus])
    assert len(np.unique(id_all))==len(id_all)
    types=np.concatenate([np.repeat(x[1],len(x[0])) for x in corpus])
    layers=np.concatenate([np.repeat(x[2],len(x[0])) for x in corpus])
    spike_time=np.concatenate([x[3] for x in corpus])*.001
    start_cols=np.r_[0,np.cumsum([len(x[0]) for x in corpus])]
    lookup=np.concatenate([x[4]+n for n,x in zip(start_cols,corpus)])
    idx=np.lexsort((spike_time,lookup))
    cvs=isi_cv(spike_time[idx],lookup[idx],len(id_all))
    for j,info in enumerate(stats):
        q=cvs[start_cols[j]:start_cols[j+1]]
        q=q[np.isfinite(q)]
        ref=cv_refs.get(info["label"],{})
        info.update(isi_cv_n=len(q),isi_cv_median=float(np.median(q)) if len(q) else None,
                    reference_isi_cv_n=int(float(ref["count"])) if ref.get("count") else None,
                    reference_isi_cv_median=float(ref["50%"]) if ref.get("50%") else None)
    receipt_ref=root/"cortical_spike_record.npz"
    np.savez_compressed(receipt_ref,neuron_ids=id_all,neuron_types=types,neuron_layers=layers,
                        spike_times_seconds=spike_time[idx],
                        spike_times_offsets=np.r_[0,np.cumsum(np.bincount(lookup,minlength=len(id_all)))])
    native_step=np.rint(spike_time/.0002).astype(int)
    np.testing.assert_allclose(spike_time,native_step*.0002,atol=1e-8,rtol=0)
    short=native_step<300_000
    k=native_step[short]//200
    assert np.all(k<1500)
    counts=np.bincount(k*len(id_all)+lookup[short],minlength=1500*len(id_all)).reshape(1500,len(id_all))
    assert counts.max()<=np.iinfo(np.uint16).max and counts.sum()==short.sum()
    count_path=root/"cortical_counts_40ms.npy"
    np.save(count_path,counts.astype(np.uint16),allow_pickle=False)
    pop=counts.sum(1)
    cells,u,v,t=components(counts)
    _,_,tau=cell_profiles(counts,cells)
    by_type={}
    type_counts={}
    for typ in ("PC","BC","MC"):
        xx=counts[:,types==typ]
        type_counts[typ]=xx.sum(1)
        q=tau[types==typ,0]
        q=q[np.isfinite(q)]
        by_type[typ]={"n":xx.shape[1],"median_tau_s":float(np.median(q)) if len(q) else None,
                      "estimated":len(q),"rate_hz":float(xx.mean()/BIN),
                      "population_dynamics":slim(population_metrics([type_counts[typ]],BIN))}
    pc=type_counts["PC"].astype(float)
    inh=(type_counts["BC"]+type_counts["MC"]).astype(float)
    lag=np.arange(-12,13)
    # Positive lag reports I counts occurring later than E counts.
    corr=[]
    for k in lag:
        e,i=(pc[:-k],inh[k:]) if k>0 else ((pc[-k:],inh[:k]) if k<0 else (pc,inh))
        corr.append(float(np.corrcoef(e,i)[0,1]) if e.std() and i.std() else None)
    ei=dict(lag_s=(lag*BIN).tolist(),correlation=corr)
    k10=native_step//50
    bins=round((stop-start)/10)
    k10=np.minimum(k10,bins-1)
    q={f"rate{i}":np.bincount(k10[use],minlength=bins)*100/int(np.count_nonzero(t)) for i,t,use in
       (("E",types=="PC",types[lookup]=="PC"),("I",types!="PC",types[lookup]!="PC"))}
    state,_=hmm(np.bincount(k10,minlength=bins))
    state=merge_short(state,5)
    phases=episode_stats(state,q,0)
    ratios=[x["mean_rate"]/x["reference_mean"] for x in stats if x["reference_mean"]>0]
    rec.update(duration_seconds=(stop-start)/1000,burn_in_seconds=start/1000,
               neuron_count_recorded=len(id_all),neuron_types=by_type,
               spike_record_sha256=digest(receipt_ref), spike_count_sha256=digest(count_path))
    return dict(receipt=rec,populations=stats,episode_proxy=phases,dynamics=slim(population_metrics([pop],BIN)),
                decomposition=profile(u,v,t),ei_synchrony=ei,voltage=cortical_vm(root,start,stop),
                median_absolute_log_rate_error=float(np.median(np.abs(np.log(ratios)))),rate_ratios=ratios),pop


def fixed_subsamples(types):
    """Keep the same cell indices and subset seeds in both author conditions."""
    for seed in (2401,2402,2403):
        rng=np.random.default_rng(seed)
        masks=[(types=="PC",8),(types!="PC",2)]
        large=[(types==kind,n) for kind,n in (("PC",395),("BC",79),("MC",26))]
        for n,groups in ((10,masks),(500,large)):
            cols=np.concatenate([rng.choice(np.flatnonzero(mask),size,replace=False)
                                 for mask,size in groups])
            yield seed,n,np.sort(cols)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--wake",type=Path,required=True)
    p.add_argument("--sleep",type=Path,required=True)
    p.add_argument("--authors",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(exist_ok=True,parents=True)
    rows={}
    series={}
    for label in ("wake","sleep"):
        rows[label],series[label]=audit(getattr(a,label),a.authors,label)
    with np.load(a.wake/"cortical_spike_record.npz") as wake, np.load(a.sleep/"cortical_spike_record.npz") as sleep:
        ids,types=wake["neuron_ids"],wake["neuron_types"]
        for field in ("neuron_ids","neuron_types","neuron_layers"):
            np.testing.assert_array_equal(wake[field],sleep[field])
    overlap=set(rows["wake"]["voltage"])&set(rows["sleep"]["voltage"])
    ratios={k:rows["sleep"]["voltage"][k]["band_power_mv2"]/rows["wake"]["voltage"][k]["band_power_mv2"]
            for k in sorted(overlap) if rows["wake"]["voltage"][k]["band_power_mv2"]>0}
    gates={"layer_rate":all(rows[s]["median_absolute_log_rate_error"]<=np.log(1.5) and
                          np.mean((np.asarray(rows[s]["rate_ratios"])>=.5)&
                                  (np.asarray(rows[s]["rate_ratios"])<=2))>=.70 for s in rows),
           "slow_modulation":bool(ratios) and np.median(list(ratios.values()))>=1.2,
           "sleep_bouts":min(rows["sleep"]["episode_proxy"][k] for k in ("up_n","down_n"))>=10}
    write(a.out/"bugnon_physiology.json",dict(conditions=rows,physiology_pass=gates,slow_band_voltage_ratios=ratios))
    print("Bugnon physiology gate:",gates,flush=True)
    if not all(gates.values()):
        raise SystemExit("Author-model data retained. No Sooter stage until its gate passes.")
    cases=[]
    subsets=list(fixed_subsamples(types))
    for cond,pop in series.items():
        m=rows[cond]["receipt"]
        run=dict(model="BugnonEtAl2024",connectivity_seed=m["global_seed"],noise_seeds=m["vp_seeds"],
                 provenance_sha256=m["spike_record_sha256"])
        entry(cases,"bugnon_"+cond,cond,cond,"bugnon_original",[pop],
              n_units=m["neuron_count_recorded"],**run)
        counts=np.load(getattr(a,cond)/"cortical_counts_40ms.npy",mmap_mode="r",allow_pickle=False)
        for seed,n,cols in subsets:
            name=f"bugnon_{cond}_n{n}_s{seed}"
            selection=dict(sample_seed=seed,sampled_gids=ids[cols].tolist(),
                           population_counts={kind:int(np.sum(types[cols]==kind))
                                              for kind in ("PC","BC","MC")})
            entry(cases,name,cond,name,"bugnon_subsample",[counts[:,cols].sum(1)],
                  n_units=n,**selection,**run)
    raw=json.dumps(dict(physiology_pass=gates,cases=cases),separators=(",",":")).encode()
    (a.out/"bugnon_counts.json.gz").write_bytes(gzip.compress(raw,mtime=0))


if __name__=="__main__":
    main()
