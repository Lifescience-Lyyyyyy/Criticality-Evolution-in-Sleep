"""Replay published networks with one weak pulse or one small parameter change."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
from pathlib import Path

import numpy as np

from analyse_snn import BIN, cell_profiles, components, entry, observed_counts, profile, slim
from metrics import population_metrics
from physiology import compare, read_record, write
from run_jercog import CONDITIONS, DT, SEEDS, compile_simulator, simulate

PARAMETERS = {
    "exc_recurrence_up10": (1.10,1.,0.,1.),
    "adaptation_up10": (1.,1.10,0.,1.),
    "drive_up_point2": (1.,1.,.2,1.),
    "inhib_recurrence_up5": (1.,1.,0.,1.05),
}


def stimulus_time(up, requested, m, kicks):
    for k in range(25,min(650,len(up))):
        t=m["burn_in_seconds"]+k*.01
        if np.all(up[k-20:k]==requested) and np.min(np.abs(kicks-t),initial=100.)>.06:
            return round(t,2)


def pulse(out, original_raw, native, exe, info, workers):
    tasks=[]
    for c,n in SEEDS:
        for cond in CONDITIONS:
            name=f"{cond}-c{c}-n{n}"
            k=native/name
            ref=json.loads((k/"manifest.json").read_text())
            with np.load(k/"states.npz",allow_pickle=False) as s:
                for label,value in (("up",True),("down",False)):
                    t=stimulus_time(s["up_mask"],value,ref,s["kick_seconds"])
                    if t is not None:
                        stop=int(np.ceil((t+1.25)/5)*5)
                        tasks.append((name,cond,c,n,label,t,stop))

    def one(case):
        name,cond,c,n,label,t,stop=case
        setting=(1.,1.,0.,1.,t,stop)
        p=simulate(exe,info,out/f"pulse_{label}",cond,(c,n),100,setting=setting)
        m,b=read_record(p)
        _,a=read_record(original_raw/name)
        boundary=round(t/DT)
        for typ in "EI":
            np.testing.assert_array_equal(a[typ][a[typ][:,0]<boundary],b[typ][b[typ][:,0]<boundary])
        k=round(t/.01)
        groups={}
        for group,typ,uid,ncells in (("target_E","E",lambda q:q<=400,401),
                                      ("untargeted_E","E",lambda q:q>400,3599),
                                      ("I","I",lambda q:np.ones_like(q,dtype=bool),1000)):
            series=[]
            for branch in (a,b):
                v=branch[typ]
                ix=uid(v[:,1])&(v[:,0]>=200*(k-10))&(v[:,0]<200*(k+100))
                f=np.bincount(v[ix,0]//200-(k-10),minlength=110)
                assert len(f)==110
                series.append(f/(.01*ncells))
            delta=series[1]-series[0]
            field=np.convolve(delta,np.ones(5)/5,mode="same")
            peak_idx=10+int(np.argmax(np.abs(field[10:35])))
            peak=field[peak_idx]
            threshold=max(.5,abs(peak)*.1)
            recovery=next((j*.01-.10 for j in range(peak_idx+1,105)
                           if np.all(np.abs(field[j:j+5])<=threshold)),None)
            groups[group]={"baseline_hz":series[0].tolist(),"pulse_hz":series[1].tolist(),
                "difference_hz":delta.tolist(),"peak_smoothed_hz":float(peak),
                "integrated_250ms_spikes_per_cell":float(delta[10:35].sum()*.01),
                "recovery_seconds":recovery,"response_threshold_hz":float(threshold)}
        return dict(condition=cond,connectivity_seed=c,noise_seed=n,pre_state=label,
                    pulse_time_s=t,pulse_duration_s=.010,pulse_current=4.,exact_prepulse=True,manifest=m,groups=groups)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        rows=list(pool.map(one,tasks))
    write(out/"pulses.json",rows)
    print("Pulse pairs verified",len(rows),flush=True)


def perturb(out, exe, info, workers):
    tasks=[(label,setting,k,n) for label,setting in PARAMETERS.items() for k,n in SEEDS]

    def one(case):
        label,setting,k,n=case
        raw=simulate(exe,info,out/"raw"/label,"figure7",(k,n),100,setting=(*setting,-1,100))
        phys=compare(raw,out/"native"/label)
        a=out/"native"/label/raw.name
        m=phys["manifest"]
        x=observed_counts(a/"simulation_record.npz",m)
        c,u,v,total=components(x)
        _,_,tau=cell_profiles(x,c)
        cells={}
        for tname,sub in (("E",slice(0,4000)),("I",slice(4000,5000))):
            vals=tau[sub,0]
            ok=np.isfinite(vals)
            cells[tname]={"valid":int(ok.sum()),"median_tau_s":float(np.median(vals[ok])) if ok.any() else None}
        series=x.sum(1)
        case_id=f"{label}-c{k}-n{n}"
        return dict(id=case_id,setting=label,manifest=m,physiology=phys["primary"],cell_timescale=cells,
                    population=slim(population_metrics([series],BIN)),decomposition=profile(u,v,total)),series

    with ThreadPoolExecutor(max_workers=workers) as pool:
        done=list(pool.map(one,tasks))
    cases=[]
    for row,series in done:
        m=row["manifest"]
        entry(cases,row["id"],row["setting"],row["id"],"parameter_perturbation",[series],
              model="JercogEtAl2017",n_units=5000,connectivity_seed=m["connectivity_seed"],
              noise_seed=m["simulation_seed"],parameter_case=m["parameter_config"])
    write(out/"perturbations.json",[row for row,_ in done])
    bundle=dict(physiology_pass={"figure7_baseline":True},cases=cases)
    (out/"perturb_counts.json.gz").write_bytes(gzip.compress(json.dumps(bundle,separators=(",",":")).encode(),mtime=0))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("pulse","parameter"))
    p.add_argument("--upstream",type=Path,required=True)
    p.add_argument("--raw",type=Path,required=True)
    p.add_argument("--native",type=Path,required=True)
    p.add_argument("--matlab-receipt",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    a=p.parse_args()
    gate=json.loads((a.native/"physiology.json").read_text())
    receipt=json.loads(a.matlab_receipt.read_text())
    assert gate["passes"]["figure7"] and receipt["valid_originals"]>=3
    root=a.out.resolve()
    exe,info=compile_simulator(a.upstream.resolve(),root/"build",interventions=True)
    if a.mode=="pulse":
        pulse(root,a.raw,a.native,exe,info,a.workers)
    else:
        perturb(root,exe,info,a.workers)


if __name__=="__main__":
    main()
