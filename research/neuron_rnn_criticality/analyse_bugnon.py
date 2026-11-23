"""Audit the author recordings before staging their matched-60-s Sooter fits."""
import argparse
import csv
import gzip
from hashlib import sha256
import json
from itertools import chain, islice
from pathlib import Path

import numpy as np
from scipy.signal import butter, find_peaks, sosfiltfilt, welch
import yaml

from analyse_snn import BIN, cell_profiles, components, entry, profile, slim
from describe_neurons import describe_cells
from metrics import population_metrics
from physiology import digest, episode_stats, hmm, merge_short, write
from run_bugnon import CASES, DENEST, NEST, REPLICATE

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


def lfp_statistics(ampa,gaba):
    # 10 ms recorded samples. Orient the current proxy like the authors' LFP.
    samples=np.arange(len(ampa))
    delayed=np.interp(samples-.6,samples,ampa)
    lfp=1.65*gaba[1:]-delayed[1:]
    f,p=welch(lfp,fs=100,nperseg=800)
    slow=(f>=.5)&(f<=2)
    filtered=sosfiltfilt(butter(4,[.5,4],btype="bandpass",fs=100,output="sos"),lfp)
    troughs=find_peaks(-filtered,height=0)[0]
    durations=[]
    for a,b in zip(troughs[:-1],troughs[1:]):
        if np.max(filtered[a:b+1])>0:
            durations.append((b-a)*.01)
    return dict(band_power=float(p[slow].sum()*(f[1]-f[0])),f=f.tolist(),psd=p.tolist(),
                candidate_waves=len(durations),
                candidate_wave_median_s=float(np.median(durations)) if durations else None)


def cortical_vm(root,t0,t1):
    # Mean Vm by recorder plus a global pyramidal-cell current proxy.
    out={}
    nbins=round((t1-t0)/10)
    currents=np.zeros((2,nbins))
    coverage=np.zeros(nbins)
    total_pc=valid_pc=0
    for r in metadata(root,"multimeter"):
        cols=r["colnames"]
        if "V_m" not in cols:
            continue
        is_pc=r["population_name"].endswith("_PC")
        total_pc+=int(is_pc)
        cols_read=(1,cols.index("V_m"))
        if is_pc:
            cols_read+=tuple(cols.index(k) for k in ("g_AMPA","g_GABA_A"))
        sums=np.zeros(nbins)
        counts=np.zeros(nbins)
        for file in source_files(root,r):
            with file.open() as stream:
                while (first:=next(stream,None)) is not None:
                    a=np.loadtxt(chain((first,),islice(stream,399999)),usecols=cols_read,ndmin=2)
                    use=(a[:,0]>=t0)&(a[:,0]<t1)
                    selected=a[use]
                    j=np.floor((selected[:,0]-t0)/10).astype(int)
                    vm=selected[:,1]
                    sums+=np.bincount(j,weights=vm,minlength=nbins)
                    count=np.bincount(j,minlength=nbins)
                    counts+=count
                    if is_pc:
                        # Vp_PC inherits E_AMPA=0 and E_GABA_A=-70 mV in both conditions.
                        currents[0]+=np.bincount(j,weights=selected[:,2]*(-vm),minlength=nbins)
                        currents[1]+=np.bincount(j,weights=selected[:,3]*(-70-vm),minlength=nbins)
                        coverage+=count
        ok=counts>0
        if ok.mean()<.98:
            continue
        valid_pc+=int(is_pc)
        v=np.interp(np.arange(nbins),np.flatnonzero(ok),sums[ok]/counts[ok])
        f,p=welch(v,fs=100,nperseg=800)
        band=(f>=.5)&(f<=2)
        out[r["label"]]=dict(band_power_mv2=float(p[band].sum()*(f[1]-f[0])),
                             sample_fraction=float(ok.mean()),f=f.tolist(),psd=p.tolist())
    lfp=dict(pyramidal_recorders=total_pc,covered_recorders=valid_pc,statistics=None)
    if total_pc and valid_pc==total_pc:
        ok=coverage>0
        means=[np.interp(np.arange(nbins),np.flatnonzero(ok),x[ok]/coverage[ok]) for x in currents]
        lfp["statistics"]=lfp_statistics(*means)
    return out,lfp


def archive_waves(path,start,stop):
    fields=("maxPeak_Amp","duration","halfSlopetoFirstPeakTime","halfSlopefromLastPeakTime")
    with path.open() as stream:
        rows=[[float(row[key]) for key in fields] for row in csv.DictReader(stream)
              if row["subnet"]=="full" and start<=float(row["startTime"])
              and float(row["endTime"])<=stop]
    return np.asarray(rows,dtype=float).reshape(-1,len(fields))


def archive_wave_stats(reference,start,stop,thresholds):
    path=reference/"df_sw.csv"
    waves=archive_waves(path,start,stop)
    results={}
    for rule,threshold in thresholds.items():
        above=waves[waves[:,0]>threshold]
        largest=above[above[:,0]>=np.quantile(above[:,0],.75)] if len(above) else above
        down_s=(largest[:,3]-largest[:,2])/1000
        results[rule]=dict(above_threshold=len(above),largest_quartile=len(largest),
                           median_wave_s=float(np.median(above[:,1])/1000) if len(above) else None,
                           median_down_half_slope_s=float(np.median(down_s)) if len(down_s) else None)
    return dict(source_sha256=digest(path),candidate_waves=len(waves),threshold_sensitivity=results,
                candidate_wave_median_s=float(np.median(waves[:,1])/1000) if len(waves) else None)


def parameter_evidence(author_path,runtime_path):
    author=yaml.full_load(author_path.read_text())
    runtime=yaml.full_load(runtime_path.read_text())
    for tree in (author,runtime):
        for key in ("input_dir","output_dir"):
            tree["simulation"]["params"].pop(key,None)
    if author!=runtime:
        raise ValueError(f"Executed parameter tree differs from the author archive: {runtime_path}")
    stages=runtime["simulation"]["params"]["sessions"]
    recorded={f"{i:02d}_{name}":runtime["session_models"][name]["params"].get("record",True)
              for i,name in enumerate(stages)}
    topology={k:runtime["network"][k] for k in ("layers","projection_models","synapse_models","topology")}
    topology["kernel"]=runtime["kernel"]
    normalized=json.dumps(topology,sort_keys=True,separators=(",",":")).encode()
    return dict(executed_parameter_sha256=digest(runtime_path),
                topology_and_kernel_spec_sha256=sha256(normalized).hexdigest(),session_recorded=recorded)


def phase_spike_counts(times,sessions):
    return {name:int(np.count_nonzero((times>=start)&(times<stop)))
            for name,(start,stop) in sessions.items()}


def paired_timescales(a,b,types):
    out={}
    for kind in ("PC","BC","MC"):
        estimates={}
        for i,name in enumerate(("exponential","integrated")):
            use=(types==kind)&np.isfinite(a[:,i])&np.isfinite(b[:,i])
            difference=b[use,i]-a[use,i]
            estimates[name]=dict(
                n_paired=len(difference),
                change_s_p10_p50_p90=np.percentile(difference,[10,50,90]).tolist() if len(difference) else None,
                median_absolute_change_s=float(np.median(np.abs(difference))) if len(difference) else None)
        out[kind]=estimates
    return out


def split_timescales(counts,types):
    first,last=np.split(counts,2)
    use=(first.sum(0)>=100)&(last.sum(0)>=100)
    if not use.any():
        empty=np.empty((0,2))
        return paired_timescales(empty,empty,types[use])
    a,b=(cell_profiles(piece,components(piece)[0])[2] for piece in (first[:,use],last[:,use]))
    return paired_timescales(a,b,types[use])


def audit(root,authors,condition,thresholds):
    rec=json.loads((root/"run_receipt.json").read_text())
    assert rec["model_commit"]==REPLICATE and rec["condition"]==condition and rec["resolution_ms"]==.2
    assert rec["denest_commit"]==DENEST and rec["nest_commit"]==NEST and rec["threads"]==5
    author_dir=authors/"simulations"/CASES[condition]
    param=author_dir/"parameter_tree.yml"
    assert rec["parameter_source"]==str(param.relative_to(authors)) and digest(param)==rec["parameter_sha256"]
    evidence=parameter_evidence(param,root/"parameter_tree.yml")
    recording=evidence["session_recorded"]
    times_path,=root.rglob("session_times.yml")
    sessions=yaml.full_load(times_path.read_text())
    author_sessions=yaml.full_load((author_dir/"session_times.yml").read_text())
    if sessions!=author_sessions or sessions.keys()!=recording.keys():
        raise ValueError(f"Session boundaries differ from the {condition} author archive: {sessions}")
    if not recording[STABLE[condition]]:
        raise ValueError(f"The {condition} stable session was not recorded")
    start,stop=sessions[STABLE[condition]]
    reference=author_dir/"analysis"
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
        times=a[:,1]
        use=(times>=start)&(times<stop)
        by_session=phase_spike_counts(times,sessions)
        assert by_session[STABLE[condition]]==np.count_nonzero(use)
        rates=np.bincount(which[use],minlength=len(ids))*1000/(stop-start)
        ref=refs[r["label"]]
        assert int(float(ref["count"]))==len(ids)
        typ=r["population_name"].rsplit("_",1)[-1]
        stats.append(dict(label=r["label"],layer=r["layer_name"],type=typ,n=len(ids),
            mean_rate=float(rates.mean()),reference_mean=float(ref["mean"]),
            quantiles=np.quantile(rates,[.25,.5,.75]).tolist(),
            reference_quantiles=[float(ref[c]) for c in ("25%","50%","75%")],
            zero_rate_fraction=float(np.mean(rates==0)),
            spike_counts_by_session={name:n if recording[name] else None for name,n in by_session.items()},
            spikes_in_unrecorded_intervals=sum(n for name,n in by_session.items() if not recording[name]),
            spikes_outside_sessions=len(times)-sum(by_session.values())))
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
    rates,acf,tau=cell_profiles(counts,cells)
    dynamics_path=root/"cortical_neuron_dynamics_40ms.npz"
    np.savez_compressed(dynamics_path,neuron_ids=id_all,neuron_types=types,neuron_layers=layers,
                        lags_seconds=np.arange(len(acf))*BIN,acf=acf.astype(np.float32),
                        firing_rate_hz=rates,tau_s=tau[:,0],integrated_tau_s=tau[:,1])
    by_type={}
    type_counts={}
    for typ in ("PC","BC","MC"):
        use=types==typ
        type_counts[typ]=counts[:,use].sum(1)
        by_type[typ]=describe_cells(rates[use],*tau[use].T)
        by_type[typ]["population_dynamics"]=slim(population_metrics([type_counts[typ]],BIN))
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
    rec.update(**evidence,session_intervals_ms=sessions,stable_session=STABLE[condition],
               duration_seconds=(stop-start)/1000,burn_in_seconds=start/1000,
               neuron_count_recorded=len(id_all),neuron_types=by_type,
               spike_record_sha256=digest(receipt_ref),spike_count_sha256=digest(count_path),
               neuron_dynamics_sha256=digest(dynamics_path))
    voltage,lfp=cortical_vm(root,start,stop)
    lfp["archive_full_network_candidates"]=archive_wave_stats(reference,start,stop,thresholds)
    return dict(receipt=rec,populations=stats,episode_proxy=phases,dynamics=slim(population_metrics([pop],BIN)),
                decomposition=profile(u,v,t),ei_synchrony=ei,voltage=voltage,lfp_proxy=lfp,
                split_half_timescales=split_timescales(counts,types),
                median_absolute_log_rate_error=float(np.median(np.abs(np.log(ratios)))),rate_ratios=ratios),pop,tau


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
    stimulus=a.authors/"simulations"/"fig0_wake_mode_stimuli"
    stimulus_times=yaml.full_load((stimulus/"session_times.yml").read_text())
    stimulus_start=stimulus_times["03_grating"][0]
    stimulus_stop=max(end for _,end in stimulus_times.values())
    stimulus_path=stimulus/"analysis"/"df_sw.csv"
    amplitudes=archive_waves(stimulus_path,stimulus_start,stimulus_stop)[:,0]
    if not len(amplitudes) or not np.all(np.isfinite(amplitudes)):
        raise ValueError(f"Missing finite full-network stimulus wave amplitudes: {stimulus_path}")
    thresholds=dict(fig4_p99=float(np.quantile(amplitudes,.99)),methods_maximum=float(amplitudes.max()))
    rows={}
    series={}
    timescales={}
    for label in ("wake","sleep"):
        rows[label],series[label],timescales[label]=audit(getattr(a,label),a.authors,label,thresholds)
    wake_receipt,sleep_receipt=(rows[k]["receipt"] for k in ("wake","sleep"))
    for field in ("global_seed","vp_seeds","topology_and_kernel_spec_sha256"):
        if wake_receipt[field]!=sleep_receipt[field]:
            raise ValueError(f"The two Bugnon conditions differ in {field}")
    with np.load(a.wake/"cortical_spike_record.npz") as wake, np.load(a.sleep/"cortical_spike_record.npz") as sleep:
        ids,types=wake["neuron_ids"],wake["neuron_types"]
        for field in ("neuron_ids","neuron_types","neuron_layers"):
            np.testing.assert_array_equal(wake[field],sleep[field])
        wake_rates=np.diff(wake["spike_times_offsets"])/rows["wake"]["receipt"]["duration_seconds"]
        sleep_rates=np.diff(sleep["spike_times_offsets"])/rows["sleep"]["receipt"]["duration_seconds"]
    rate_corr={}
    for kind,reference in (("PC",.90),("BC",.80),("MC",None)):
        mask=types==kind
        w,s=wake_rates[mask],sleep_rates[mask]
        pearson=float(np.corrcoef(w,s)[0,1]) if w.std() and s.std() else None
        rate_corr[kind]=dict(n=len(w),pearson=pearson,published_fig4a_r=reference)
    overlap=set(rows["wake"]["voltage"])&set(rows["sleep"]["voltage"])
    ratios={k:rows["sleep"]["voltage"][k]["band_power_mv2"]/rows["wake"]["voltage"][k]["band_power_mv2"]
            for k in sorted(overlap) if rows["wake"]["voltage"][k]["band_power_mv2"]>0}
    gates={"layer_rate":all(rows[s]["median_absolute_log_rate_error"]<=np.log(1.5) and
                          np.mean((np.asarray(rows[s]["rate_ratios"])>=.5)&
                                  (np.asarray(rows[s]["rate_ratios"])<=2))>=.70 for s in rows),
           "slow_modulation":bool(ratios) and np.median(list(ratios.values()))>=1.2,
           "sleep_bouts":min(rows["sleep"]["episode_proxy"][k] for k in ("up_n","down_n"))>=10}
    write(a.out/"bugnon_physiology.json",
          dict(conditions=rows,physiology_pass=gates,slow_band_voltage_ratios=ratios,
               paired_cell_rate_correlations=rate_corr,
               paired_timescale_change_sleep_minus_wake=paired_timescales(timescales["wake"],timescales["sleep"],types),
               archive_stimulus_thresholds=dict(source_sha256=digest(stimulus_path),
                                                stimulus_window_ms=[stimulus_start,stimulus_stop],
                                                candidate_waves=len(amplitudes),thresholds=thresholds)))
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
