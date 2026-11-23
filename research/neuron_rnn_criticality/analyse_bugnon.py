"""Audit the author recordings before staging their matched-60-s Sooter fits."""
import argparse
import csv
import gzip
from hashlib import sha256
import json
from itertools import chain, islice
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import butter, find_peaks, sosfiltfilt, welch
import yaml

from analyse_snn import BIN, cell_profiles, components, entry, profile, slim
from describe_neurons import describe_cells
from inputs import shrink
from metrics import population_metrics
from physiology import digest, episode_stats, hmm, merge_short, write
from run_bugnon import CASES, DENEST, NEST, REPLICATE
from run_bugnon_paper import COMMIT as PAPER_COMMIT, session_overrides

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


def lfp_proxy(ampa,gaba):
    # 10 ms recorded samples. Orient the current proxy like the authors' LFP.
    samples=np.arange(len(ampa))
    delayed=np.interp(samples-.6,samples,ampa)
    return 1.65*gaba[1:]-delayed[1:]


def positive_waves(lfp):
    filtered=sosfiltfilt(butter(4,[.5,4],btype="bandpass",fs=100,output="sos"),lfp)
    troughs=find_peaks(-filtered,height=0)[0]
    waves=[]
    for a,b in zip(troughs[:-1],troughs[1:]):
        peak=a+np.argmax(filtered[a:b+1])
        if filtered[peak]>0:
            waves.append((int(a),int(peak),int(b),float(filtered[peak])))
    return waves


def lfp_statistics(ampa,gaba):
    lfp=lfp_proxy(ampa,gaba)
    f,p=welch(lfp,fs=100,nperseg=800)
    slow=(f>=.5)&(f<=2)
    durations=[(b-a)*.01 for a,_,b,_ in positive_waves(lfp)]
    return dict(band_power=float(p[slow].sum()*(f[1]-f[0])),f=f.tolist(),psd=p.tolist(),
                candidate_waves=len(durations),
                candidate_wave_median_s=float(np.median(durations)) if durations else None)


def transition_crossings(voltage,peak):
    """Last downward and first upward halfway crossings within a single wave."""
    finite=np.all(np.isfinite(voltage),axis=0)
    v=median_filter(voltage[:,finite],size=(3,1),mode="nearest")
    first,low,last=v[0],v[peak],v[-1]
    down,up=(first+low)/2,(last+low)/2
    down_steps=(v[:peak]>=down)&(v[1:peak+1]<down)
    up_steps=(v[peak:-1]<=up)&(v[peak+1:]>up)
    use=(first>low)&(last>low)&down_steps.any(0)&up_steps.any(0)
    down_index=peak-np.argmax(down_steps[:,use][::-1],axis=0)
    up_index=peak+1+np.argmax(up_steps[:,use],axis=0)
    return down_index,up_index


def local_transition_synchrony(currents,coverage,traces,gids,t0):
    if len(traces)!=2 or not np.any(coverage):
        return {"status":"missing_L4_pyramidal_recordings"}
    ok=coverage>0
    means=[np.interp(np.arange(len(x)),np.flatnonzero(ok),x[ok]/coverage[ok]) for x in currents]
    voltage=np.concatenate(traces,axis=1)[1:]
    waves=positive_waves(lfp_proxy(*means))
    selected=sorted(sorted(waves,key=lambda row:-row[3])[:(len(waves)+1)//2])
    rows=[]
    for a,p,b,amp in selected:
        down,up=transition_crossings(voltage[a:b+1],p-a)
        n=len(down)
        rows.append(dict(peak_time_ms=t0+(p+1)*10,lfp_peak=amp,n_paired_cells=n,
                         up_to_down_sd_ms=float(down.std(ddof=1)*10) if n>1 else None,
                         down_to_up_sd_ms=float(up.std(ddof=1)*10) if n>1 else None))
    both=[x["down_to_up_sd_ms"]-x["up_to_down_sd_ms"] for x in rows if x["n_paired_cells"]>1]
    return dict(status="measured" if both else "no_valid_paired_events",selected_neuron_ids=gids,
                sampling_ms=10,unthresholded_candidates=len(waves),selected_events=rows,
                n_valid_paired_waves=len(both),
                paired_sd_change_ms_p25_p50_p75=np.percentile(both,[25,50,75]).tolist() if both else None)


def cortical_vm(root,t0,t1,up=None,local=False):
    # Mean Vm by recorder plus a global pyramidal-cell current proxy.
    out={}
    nbins=round((t1-t0)/10)
    currents=np.zeros((2,nbins))
    coverage=np.zeros(nbins)
    local_currents=np.zeros((2,nbins))
    local_coverage=np.zeros(nbins)
    local_traces,local_gids=[],[]
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
        local_group=local and is_pc and r["layer_name"].startswith("Vp_L4_")
        if local_group:
            nx,ny=r["layer_shape"]
            gids=[gid for gid,(x,y,_) in sorted(r["locations"].items())
                  if abs((x+.5)/nx-.5)<.05 and abs((y+.5)/ny-.5)<.05]
            lookup=np.full(max(r["gids"])+1,-1,dtype=int)
            lookup[gids]=np.arange(len(gids))
            trace=np.full((nbins,len(gids)),np.nan,dtype=np.float32)
            cols_read+=(cols.index("gid"),)
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
                        ampa,gaba=selected[:,2]*(-vm),selected[:,3]*(-70-vm)
                        currents[0]+=np.bincount(j,weights=ampa,minlength=nbins)
                        currents[1]+=np.bincount(j,weights=gaba,minlength=nbins)
                        coverage+=count
                    if local_group:
                        index=lookup[selected[:,-1].astype(int)]
                        inside=index>=0
                        bins=j[inside]
                        trace[bins,index[inside]]=vm[inside]
                        local_currents[0]+=np.bincount(bins,weights=ampa[inside],minlength=nbins)
                        local_currents[1]+=np.bincount(bins,weights=gaba[inside],minlength=nbins)
                        local_coverage+=np.bincount(bins,minlength=nbins)
        if local_group:
            local_traces.append(trace)
            local_gids.extend(gids)
        ok=counts>0
        if ok.mean()<.98:
            continue
        valid_pc+=int(is_pc)
        v=np.interp(np.arange(nbins),np.flatnonzero(ok),sums[ok]/counts[ok])
        f,p=welch(v,fs=100,nperseg=800)
        band=(f>=.5)&(f<=2)
        delta_mv=None
        if up is not None:
            recorded_up,recorded_down=up&ok,~up&ok
            if recorded_up.any() and recorded_down.any():
                delta_mv=float(v[recorded_up].mean()-v[recorded_down].mean())
        out[r["label"]]=dict(band_power_mv2=float(p[band].sum()*(f[1]-f[0])),
                             sample_fraction=float(ok.mean()),f=f.tolist(),psd=p.tolist(),
                             up_minus_down_mv=delta_mv)
    lfp=dict(pyramidal_recorders=total_pc,covered_recorders=valid_pc,statistics=None)
    if total_pc and valid_pc==total_pc:
        ok=coverage>0
        means=[np.interp(np.arange(nbins),np.flatnonzero(ok),x[ok]/coverage[ok]) for x in currents]
        lfp["statistics"]=lfp_statistics(*means)
    if local:
        lfp["local_transition_proxy"]=local_transition_synchrony(
            local_currents,local_coverage,local_traces,local_gids,t0)
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
    return tree_evidence(runtime,runtime_path)


def tree_evidence(runtime,runtime_path):
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


def paper_run_evidence(rec,condition,runtime_path,sessions):
    if rec.get("condition")!=condition or rec.get("source_commit")!=PAPER_COMMIT:
        raise ValueError(f"Expected a completed paper-linked {condition} simulation")
    runtime=yaml.full_load(runtime_path.read_text())
    want=session_overrides(condition)
    stages=want["simulation"]["params"]["sessions"]
    assert rec["override"]==want and runtime["simulation"]["params"]["sessions"]==stages
    assert rec["runtime_parameter_sha256"]==digest(runtime_path)
    assert list(rec["sessions"])==list(sessions) and all(tuple(rec["sessions"][s])==v for s,v in sessions.items())
    assert rec["nest"]=="NEST nest-2.20.0"
    kernel=runtime["kernel"]["nest_params"]
    assert kernel["resolution"]==.2 and kernel["local_num_threads"]==5
    expected=(2000.,122000.) if condition=="wake" else (8000.,128000.)
    stable=f"{len(stages)-1:02d}_{stages[-1]}"
    if sessions[stable]!=expected:
        raise ValueError(f"Wrong stable interval for {condition}: {sessions[stable]}")
    return stable,tree_evidence(runtime,runtime_path)


def audit(root,rec,authors,condition,thresholds):
    times_path,=root.rglob("session_times.yml")
    sessions=yaml.full_load(times_path.read_text())
    runtime=root/"parameter_tree.yml"
    refs,cv_refs={},{}
    if authors:
        assert rec["model_commit"]==REPLICATE and rec["condition"]==condition and rec["resolution_ms"]==.2
        assert rec["denest_commit"]==DENEST and rec["nest_commit"]==NEST and rec["threads"]==5
        author_dir=authors/"simulations"/CASES[condition]
        param=author_dir/"parameter_tree.yml"
        assert rec["parameter_source"]==str(param.relative_to(authors)) and digest(param)==rec["parameter_sha256"]
        evidence=parameter_evidence(param,runtime)
        stable=STABLE[condition]
        if sessions!=yaml.full_load((author_dir/"session_times.yml").read_text()):
            raise ValueError(f"Session boundaries differ from the {condition} author archive: {sessions}")
        reference=author_dir/"analysis"
        with (reference/"df_rate_stats.csv").open() as f:
            refs={r["metadata_filename"]:r for r in csv.DictReader(f) if r["session"]==stable}
        with (reference/"df_cv_stats.csv").open() as f:
            cv_refs={r["metadata_filename"]:r for r in csv.DictReader(f)
                     if r["session"]==stable and r["min_n_spikes"]=="10"
                     and r["duration"]=="0" and r["t_min"]=="0"}
    else:
        stable,evidence=paper_run_evidence(rec,condition,runtime,sessions)
    recording=evidence["session_recorded"]
    if sessions.keys()!=recording.keys() or not recording[stable]:
        raise ValueError(f"Wrong recording settings for stable {condition}: {recording}")
    start,stop=sessions[stable]
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
        assert by_session[stable]==np.count_nonzero(use)
        rates=np.bincount(which[use],minlength=len(ids))*1000/(stop-start)
        ref=refs[r["label"]] if authors else None
        if ref:
            assert int(float(ref["count"]))==len(ids)
        typ=r["population_name"].rsplit("_",1)[-1]
        stats.append(dict(label=r["label"],layer=r["layer_name"],type=typ,n=len(ids),
            mean_rate=float(rates.mean()),reference_mean=float(ref["mean"]) if ref else None,
            quantiles=np.quantile(rates,[.25,.5,.75]).tolist(),
            reference_quantiles=[float(ref[c]) for c in ("25%","50%","75%")] if ref else None,
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
    ratios=[x["mean_rate"]/x["reference_mean"] for x in stats if x["reference_mean"]]
    rec.update(**evidence,session_intervals_ms=sessions,stable_session=stable,
               duration_seconds=(stop-start)/1000,burn_in_seconds=start/1000,
               neuron_count_recorded=len(id_all),neuron_types=by_type,
               spike_record_sha256=digest(receipt_ref),spike_count_sha256=digest(count_path),
               neuron_dynamics_sha256=digest(dynamics_path))
    voltage,lfp=cortical_vm(root,start,stop,state,local=not authors and condition=="sleep")
    if authors:
        lfp["archive_full_network_candidates"]=archive_wave_stats(reference,start,stop,thresholds)
    return dict(receipt=rec,populations=stats,episode_proxy=phases,dynamics=slim(population_metrics([pop],BIN)),
                decomposition=profile(u,v,t),ei_synchrony=ei,voltage=voltage,lfp_proxy=lfp,
                split_half_timescales=split_timescales(counts,types),
                median_absolute_log_rate_error=float(np.median(np.abs(np.log(ratios)))) if ratios else None,
                rate_ratios=ratios),pop,tau


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


def prepare_condition(root,receipt,authors,condition,thresholds):
    summary,pop,timescales=audit(root,receipt,authors,condition,thresholds)
    with np.load(root/"cortical_spike_record.npz",allow_pickle=False) as spikes:
        info={field:spikes[field] for field in ("neuron_ids","neuron_types","neuron_layers")}
        info["rates_hz"]=np.diff(spikes["spike_times_offsets"])/summary["receipt"]["duration_seconds"]
    counts=np.load(root/"cortical_counts_40ms.npy",mmap_mode="r",allow_pickle=False)
    sampled=[counts[:,cols].sum(1) for _,_,cols in fixed_subsamples(info["neuron_types"])]
    return dict(summary=summary,timescales=timescales,counts=np.stack([pop,*sampled]),**info)


def save_prepared_condition(path,prepared):
    checksum=path.with_suffix(".sha256")
    if path.exists() or checksum.exists():
        raise FileExistsError(path)
    tmp=path.with_suffix(".tmp.npz")
    arrays={k:v for k,v in prepared.items() if k!="summary"}
    summary=json.dumps(shrink(prepared["summary"]),allow_nan=False,separators=(",",":"))
    np.savez_compressed(tmp,summary=summary,**arrays)
    tmp.replace(path)
    checksum.write_text(digest(path)+"\n")


def load_prepared_condition(path,condition):
    checksum=digest(path)
    if checksum!=path.with_suffix(".sha256").read_text().strip():
        raise ValueError(f"Checksum mismatch for {path}")
    with np.load(path,allow_pickle=False) as archive:
        data={k:archive[k] for k in archive.files if k!="summary"}
        data["summary"]=json.loads(archive["summary"].item())
    receipt=data["summary"]["receipt"]
    if receipt.get("source_commit")!=PAPER_COMMIT or receipt.get("condition")!=condition or receipt.get("duration_seconds")!=120:
        raise ValueError(f"Expected a completed 120 s paper-linked {condition} recording in {path}")
    n=receipt["neuron_count_recorded"]
    assert all(data[k].shape==(n,) for k in ("neuron_ids","neuron_types","neuron_layers","rates_hz"))
    assert data["timescales"].shape==(n,2)
    counts=data["counts"]
    assert counts.shape==(7,1500) and np.issubdtype(counts.dtype,np.integer) and counts.min()>=0
    assert (counts[1:]<=counts[0]).all()
    np.testing.assert_allclose(counts[0].var(),data["summary"]["decomposition"]["total"][0],atol=1e-7)
    receipt["transfer_bundle_sha256"]=checksum
    return data


def paper_physiology_gates(rows,rates,types,rate_corr):
    # Directional checks follow the reported regime without a numeric tolerance for Figure 4.
    pc,bc=types=="PC",types=="BC"
    distributions=all(np.median(a[pc])<np.median(a[bc]) and a[pc].mean()>np.median(a[pc]) for a in rates)
    paired=all(rate_corr[k]["pearson"] is not None and rate_corr[k]["pearson"]>0 for k in ("PC","BC"))
    wake_lfp,sleep_lfp=(rows[k]["lfp_proxy"]["statistics"] for k in ("wake","sleep"))
    slow_waves=bool(wake_lfp and sleep_lfp and sleep_lfp["band_power"]>wake_lfp["band_power"])
    sleep=rows["sleep"]
    bouts=sleep["episode_proxy"]
    vm=[z["up_minus_down_mv"] for key,z in sleep["voltage"].items()
        if key.endswith("_PC") and z["up_minus_down_mv"] is not None]
    alternation=min(bouts["up_n"],bouts["down_n"])>0
    silent_down=alternation and all(bouts[k]["up_hz"]>bouts[k]["down_hz"] for k in ("E","I"))
    return dict(layer_rate=bool(distributions and paired),slow_modulation=slow_waves,
                sleep_bouts=bool(silent_down and vm and np.median(vm)>0))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--wake",type=Path)
    p.add_argument("--sleep",type=Path)
    source=p.add_mutually_exclusive_group(required=True)
    source.add_argument("--authors",type=Path,help="Archives for the legacy release")
    source.add_argument("--paper",action="store_true",help="Completed paper-linked mGluR runs")
    source.add_argument("--prepare",choices=("wake","sleep"),help="Audit one completed paper-linked run for transfer")
    source.add_argument("--prepared",action="store_true",help="Combine two checked condition .npz files")
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    def read_receipt(condition,root):
        path=root/"run_receipt.json"
        if not path.is_file():
            p.error(f"Awaiting the completed {condition} NEST receipt: {path}")
        return json.loads(path.read_text())

    if a.prepare:
        condition=a.prepare
        root=getattr(a,condition)
        if root is None:
            p.error(f"--prepare {condition} needs --{condition} with the completed native directory")
        a.out.mkdir(exist_ok=True,parents=True)
        path=a.out/f"bugnon-{condition}-prepared.npz"
        if path.exists() or path.with_suffix(".sha256").exists():
            raise FileExistsError(path)
        prepared=prepare_condition(root,read_receipt(condition,root),None,condition,None)
        save_prepared_condition(path,prepared)
        print("Condition audit archived:",path,flush=True)
        return
    if a.wake is None or a.sleep is None:
        p.error("Both --wake and --sleep are required for a paired analysis")
    roots={"wake":a.wake,"sleep":a.sleep}
    paper=a.paper or a.prepared
    receipts={} if a.prepared else {k:read_receipt(k,r) for k,r in roots.items()}
    a.out.mkdir(exist_ok=True,parents=True)
    bundle=a.out/"bugnon_counts.json.gz"
    if bundle.exists():
        raise FileExistsError(f"Use a fresh output directory to preserve the verified bundle: {bundle}")
    thresholds=None
    archive_info=None
    if a.authors:
        stimulus=a.authors/"simulations"/"fig0_wake_mode_stimuli"
        stimulus_times=yaml.full_load((stimulus/"session_times.yml").read_text())
        stimulus_start=stimulus_times["03_grating"][0]
        stimulus_stop=max(end for _,end in stimulus_times.values())
        stimulus_path=stimulus/"analysis"/"df_sw.csv"
        amplitudes=archive_waves(stimulus_path,stimulus_start,stimulus_stop)[:,0]
        if not len(amplitudes) or not np.all(np.isfinite(amplitudes)):
            raise ValueError(f"Missing finite full-network stimulus wave amplitudes: {stimulus_path}")
        thresholds=dict(fig4_p99=float(np.quantile(amplitudes,.99)),methods_maximum=float(amplitudes.max()))
        archive_info=dict(source_sha256=digest(stimulus_path),stimulus_window_ms=[stimulus_start,stimulus_stop],
                          candidate_waves=len(amplitudes),thresholds=thresholds)
    if a.prepared:
        data={k:load_prepared_condition(path,k) for k,path in roots.items()}
    else:
        data={k:prepare_condition(path,receipts[k],a.authors,k,thresholds) for k,path in roots.items()}
    rows={k:v["summary"] for k,v in data.items()}
    wake_receipt,sleep_receipt=(rows[k]["receipt"] for k in ("wake","sleep"))
    fields=("global_seed","vp_seeds","topology_and_kernel_spec_sha256")
    if paper:
        fields+=("extension_sha256","initialization")
    for field in fields:
        if wake_receipt[field]!=sleep_receipt[field]:
            raise ValueError(f"The two Bugnon conditions differ in {field}")
    wake,sleep=(data[k] for k in ("wake","sleep"))
    ids,types=wake["neuron_ids"],wake["neuron_types"]
    for field in ("neuron_ids","neuron_types","neuron_layers"):
        np.testing.assert_array_equal(wake[field],sleep[field])
    wake_rates,sleep_rates=wake["rates_hz"],sleep["rates_hz"]
    rate_corr={}
    for kind,reference in (("PC",.90),("BC",.80),("MC",None)):
        mask=types==kind
        w,s=wake_rates[mask],sleep_rates[mask]
        pearson=float(np.corrcoef(w,s)[0,1]) if w.std() and s.std() else None
        rate_corr[kind]=dict(n=len(w),pearson=pearson,published_fig4a_r=reference,
                             delta_vs_published=pearson-reference if pearson is not None and reference is not None else None)
    overlap=set(rows["wake"]["voltage"])&set(rows["sleep"]["voltage"])
    ratios={k:rows["sleep"]["voltage"][k]["band_power_mv2"]/rows["wake"]["voltage"][k]["band_power_mv2"]
            for k in sorted(overlap) if rows["wake"]["voltage"][k]["band_power_mv2"]>0}
    if paper:
        gates=paper_physiology_gates(rows,(wake_rates,sleep_rates),types,rate_corr)
    else:
        gates={"layer_rate":all(rows[s]["median_absolute_log_rate_error"] is not None and
                               rows[s]["median_absolute_log_rate_error"]<=np.log(1.5) and
                               np.mean((np.asarray(rows[s]["rate_ratios"])>=.5)&
                                       (np.asarray(rows[s]["rate_ratios"])<=2))>=.70 for s in rows),
               "slow_modulation":bool(ratios and np.median(list(ratios.values()))>=1.2),
               "sleep_bouts":min(rows["sleep"]["episode_proxy"][k] for k in ("up_n","down_n"))>=10}
    case="paper_mglur" if paper else "legacy_archive"
    write(a.out/"bugnon_physiology.json",
          dict(parameter_case=case,conditions=rows,physiology_pass=gates,slow_band_voltage_ratios=ratios,
               paired_cell_rate_correlations=rate_corr,
               paired_timescale_change_sleep_minus_wake=paired_timescales(wake["timescales"],sleep["timescales"],types),
               archive_stimulus_thresholds=archive_info))
    print("Bugnon physiology gate:",gates,flush=True)
    if not all(gates.values()):
        raise SystemExit("Author-model data retained. No Sooter stage until its gate passes.")
    cases=[]
    subsets=list(fixed_subsamples(types))
    for cond,condition_data in data.items():
        m=rows[cond]["receipt"]
        pop,*sample_counts=condition_data["counts"]
        run=dict(model="BugnonEtAl2024",parameter_case=case,connectivity_seed=m["global_seed"],noise_seeds=m["vp_seeds"],
                 provenance_sha256=m["spike_record_sha256"])
        entry(cases,"bugnon_"+cond,cond,cond,"bugnon_original",[pop],
              n_units=m["neuron_count_recorded"],**run)
        for (seed,n,cols),counts in zip(subsets,sample_counts,strict=True):
            name=f"bugnon_{cond}_n{n}_s{seed}"
            selection=dict(sample_seed=seed,sampled_gids=ids[cols].tolist(),
                           population_counts={kind:int(np.sum(types[cols]==kind))
                                              for kind in ("PC","BC","MC")})
            entry(cases,name,cond,name,"bugnon_subsample",[counts],
                  n_units=n,**selection,**run)
    raw=json.dumps(dict(physiology_pass=gates,cases=cases),separators=(",",":")).encode()
    temp=bundle.with_suffix(".tmp")
    temp.write_bytes(gzip.compress(raw,mtime=0))
    temp.replace(bundle)


if __name__=="__main__":
    main()
