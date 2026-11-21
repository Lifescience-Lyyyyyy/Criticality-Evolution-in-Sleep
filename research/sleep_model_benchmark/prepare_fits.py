"""Prepare benchmark-only count fixtures and observation/stationarity diagnostics."""
import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from scipy.signal import correlate

from observation import rate_counts
from validate_original import FORMAL_SEEDS, PARAMETERS, TARGETS, dwells, simulate

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts/data"))
from prepare_sleep_dandi import save_prepared

ANCHOR = "rat_frontal"
BASE = dict(bin=.04, order=5, alternate=10)
SPECS = {
    "primary": BASE,
    "ms20": dict(bin=.02, order=5, alternate=10),
    "ms80_low": dict(bin=.08, order=5, alternate=2),
    "ms80_high": dict(bin=.08, order=3, alternate=3),
    "ms200": dict(bin=.20, order=5, alternate=1),
}
SOURCE = "results/sleep_dandi/000041/cortex.json"


def write(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def acf(x, k=20):
    centered = x - np.mean(x)
    c = correlate(centered, centered, mode="full", method="fft")[len(x)-1:len(x)+k]
    return (c[1:]/c[0]).tolist() if c[0]>0 else [0.]*k


def runs(x, predicate, bin_s):
    mask = predicate(x)
    cuts = np.flatnonzero(np.diff(np.r_[False, mask, False]))
    # Exclude events censored by trace boundaries.
    lengths = (cuts[1::2]-cuts[::2])[(cuts[::2]>0)&(cuts[1::2]<len(x))]
    return dict(events=len(lengths),fraction=float(mask.mean()),
                median_s=float(np.median(lengths)*bin_s) if len(lengths) else None,
                p95_s=float(np.quantile(lengths,.95)*bin_s) if len(lengths) else None)


def stationarity(x, bin_s, seed):
    """Screen first/last blocks for a mean or variance shift using a block permutation."""
    block = round(10/bin_s)
    usable = len(x)//(2*block)*(2*block)
    if usable<12*block:
        return {"pass_":False,"reason":"fewer than 12 complete 10 s blocks"}
    grouped = x[:usable].reshape(-1,block)
    mean = grouped.mean(1)
    logvar = np.log(np.maximum(grouped.var(1), np.finfo(float).tiny))
    values = np.stack((mean, logvar), axis=1)
    n = len(values)//2
    diff = np.abs(values[:n].mean(0)-values[n:].mean(0))
    rng = np.random.default_rng(seed)
    permuted = np.array([rng.permutation(values) for _ in range(499)])
    random_diff = np.abs(permuted[:,:n].mean(1)-permuted[:,n:].mean(1))
    p = (1+(random_diff>=diff).sum(0))/500
    mean_shift = float(diff[0] / max(np.std(x[:usable]),1e-12))
    logvar_shift = float(diff[1])
    drift = (p[0] <= .01 and mean_shift>.25) or (p[1] <= .01 and logvar_shift>np.log(2))
    return dict(pass_=bool(not drift), blocks=len(values), p_mean=float(p[0]),p_variance=float(p[1]),
                mean_shift_in_sd=mean_shift, variance_ratio=float(np.exp(min(logvar_shift,30))))


def empirical_reference():
    path = ROOT/SOURCE
    d = json.loads(path.read_text())
    cohorts = defaultdict(lambda: defaultdict(list))
    for row in d["records"]:
        m,f=row["record"],row["fit"]
        if m["region"]!="mPFC" or m["variant"]!="reference" or m["context"]!="HOME":
            continue
        if not (f["exit_status"]==0 and f["solver_flag"]>0 and f["min_ar_root"]>1):
            continue
        units=len(m["unit_ids"])
        cohorts[m["subject"]][m["state"]].append(dict(session=m["session"],units=units,
                    seconds=f["duration_s"],unit_hz=f["rate_hz"]/units,d2=f["d2"]))
    # Unit count and mean rate anchor observation exposure, avoiding use of empirical d2.
    nrem=[v["NREM"] for v in cohorts.values() if v["NREM"]]
    rate=float(np.median([np.mean([z["unit_hz"] for z in xs]) for xs in nrem]))
    units=round(np.median([np.median([z["units"] for z in xs]) for xs in nrem]))
    return dict(source=SOURCE,sha256=sha256(path.read_bytes()).hexdigest(),dataset="000041 mPFC",
                session_data=dict(cohorts),units=units,calibration_unit_hz=rate,
                source_overlap="DANDI 000041/CRCNS fcx-1 Watson archive; descriptive comparison only")


def physiological_error(trace, sample_s):
    d=dwells(trace)
    if d["status"]!="OK":
        return d
    for label in ("up","down"):
        d[label+"_mean_s"]=sample_s*d[label+"_mean_au"]
    d["z_vs_published"]={k:(d[k]-mean)/sd for k,(mean,sd) in TARGETS[ANCHOR].items()}
    return d


def iaaft(x, seed, iterations=200):
    """Preserve the count marginal exactly and approximate the univariate spectrum."""
    rng=np.random.default_rng(seed)
    sorted_x=np.sort(x)
    target=np.abs(np.fft.rfft(x))
    current=rng.permutation(x)
    for _ in range(iterations):
        angles=np.angle(np.fft.rfft(current))
        y=np.fft.irfft(target*np.exp(1j*angles),n=len(x))
        ranks=np.argsort(y,kind="stable")
        current[ranks]=sorted_x
    return current


def build(cache, dest):
    doc=json.loads((cache/"original_validation.json").read_text())
    original=doc["results"][ANCHOR]
    if original["gate"]!="PASS_SUMMARY":
        raise RuntimeError("Original physiological summary did not pass the benchmark gate")
    empirical=empirical_reference()
    sample_s=original["time_scale_s"]
    u=empirical["units"]
    calibration=[simulate(cache/"levenstein",s,PARAMETERS[ANCHOR]).mean()
                 for s in doc["calibration_seeds"]]
    gain=empirical["calibration_unit_hz"]/float(np.mean(calibration))
    params=PARAMETERS[ANCHOR]
    strata={key:[] for key in SPECS}
    metadata=[]
    latent_controls={
        "scan_recurrence_low":({"w":params["w"]-.40},None),
        "scan_recurrence_high":({"w":params["w"]+.40},None),
        "scan_drive_low":({"I":params["I"]-.15},None),
        "scan_drive_high":({"I":params["I"]+.15},None),
        "scan_adaptation_low":({"b":.75},None),
        "scan_adaptation_high":({"b":1.25},None),
        "scan_noise_low":({"sigma":.20},None),
        "scan_noise_high":({"sigma":.30},None),
        "gate_drive_only":({},[.40*float(np.mean(calibration)),0,500]),
        "gate_recurrence_only":({},[0,.40,500]),
        "recurrence_and_drive_gate":({"w":params["w"]-.40},[.40*float(np.mean(calibration)),0,500]),
    }
    for seed in FORMAL_SEEDS:
        trace=np.load(cache/f"{ANCHOR}-{seed}.npy",allow_pickle=False)
        def counts(source, *, shape=None, rate_multiplier=1):
            return rate_counts([source*gain*rate_multiplier],sample_s,bin_s=.02,units=u,
                               seed=100000+seed,shared_shape=shape)

        def append(label,x,spec,observation,physio,extra=None):
            b=SPECS[spec]["bin"]
            if len(x)*b<120 or np.var(x)<=0 or not (x>=0).all():
                raise ValueError(f"Insufficient variable count data: {seed} {label} {spec}")
            key = f"{label}-{seed}-{spec}"
            burst_cut=int(np.quantile(x,.90))
            m=dict(key=key,asset_id=f"levenstein-{seed}",subject=f"simulation-{seed}",
                   session=f"levenstein-{seed}",region="rat_frontal_model",context=label,
                   variant="reference",state="NREM" if label.startswith("original_") else "unlabeled",
                   seed=seed,bin_seconds=b,spec=spec,unit_count=u)
            strata[spec].append((m,cells(x)))
            metadata.append(dict(record=m,observation=observation,physiology=physio,
                gate_schedule=(extra or {}).get("parameters",{}).get("gating"),
                stationarity=stationarity(x,b,seed),count_acf=acf(x),
                bursts=runs(x,lambda z:z>burst_cut,b),silent=runs(x,lambda z:z==0,b),
                burst_threshold=burst_cut,extra=extra or {}))

        base=counts(trace)
        physiology=physiological_error(trace,sample_s)
        for spec,config in SPECS.items():
            obs=base.rebin(round(config["bin"]/.02))
            append("original_poisson",obs.population()[0],spec,obs.metadata,physiology)
        for label, obs in (("original_shared_gamma",counts(trace,shape=2)),
                           ("original_half_rate",counts(trace,rate_multiplier=.5)),
                           ("original_double_rate",counts(trace,rate_multiplier=2))):
            obs=obs.rebin(2)
            append(label,obs.population()[0],"primary",obs.metadata,physiology)
        for label,(change,gate) in latent_controls.items():
            settings=params|change
            other=simulate(cache/"levenstein",seed,settings,gate=gate)
            obs=counts(other).rebin(2)
            append(label,obs.population()[0],"primary",obs.metadata,physiological_error(other,sample_s),
                   dict(parameters=settings|{"gating":gate}))
        anchor=base.rebin(2)
        x=anchor.population()[0]
        y=iaaft(x,seed)
        error=float(np.sqrt(np.mean((np.array(acf(x))-acf(y))**2)))
        append("matched_count_surrogate",y,"primary",
               {"mode":"IAAFT_count_surrogate","seed":seed,"units":u},None,
               dict(acf_rmse=error,marginal_exact=bool(np.array_equal(np.sort(x),np.sort(y))),
                    second_order_gate=bool(error<=.02),perturbation_defined=False))
    dest.mkdir(parents=True,exist_ok=True)
    plan=[]
    for name,spec in SPECS.items():
        folder=dest/name
        info=save_prepared(strata[name],folder/"counts.mat","states",
                           [f"levenstein-{s}" for s in FORMAL_SEEDS])
        plan.append(dict(name=name,receipt=str(folder/"counts.json"),sha256=info["sha256"],
                         count=len(strata[name]),analysis_spec=list(spec.values())))
    write(dest/"fit_plan.json",plan)
    write(dest/"observation_qc.json",dict(empirical=empirical,sample_s=sample_s,
         rate_scale_hz=gain,rate_calibration_seeds=doc["calibration_seeds"],
         benchmark_burn_AU=doc["burn_au"],records=metadata))
    print(f"Verified synthetic fit inputs: {sum(len(v) for v in strata.values())} rows",flush=True)


def cells(x):
    c=np.empty(1,dtype=object)
    c[0]=x.astype(np.int32)
    return c


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--simulation",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    build(args.simulation,args.out)
