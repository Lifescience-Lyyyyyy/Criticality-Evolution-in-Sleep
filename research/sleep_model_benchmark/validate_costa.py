"""Continuous-signal reproduction using Costa's pinned original C++ equations."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from time import perf_counter

import numpy as np
from scipy.io import loadmat
from scipy.signal import welch

from observation import continuous_signal

PIN="eb9fad3bf9e19e3e7f3018ec471c93cfb0f31ebe"
SEEDS=(1441,1442,1443)
FS=100.


def summarize(trace):
    series=continuous_signal([trace],FS,"mV")
    freq,power=welch(series.segments[0],fs=FS,nperseg=1000,noverlap=500)
    out={}
    for name,(lo,hi) in {"slow":(.2,1.),"spindle":(10,16)}.items():
        ix=np.flatnonzero((freq>=lo)&(freq<=hi))
        highest=ix[np.argmax(power[ix])]
        out[name]=dict(peak_hz=float(freq[highest]), bin_edges_hz=[lo,hi],
                       band_power_mV2=float(np.sum(power[ix])*(freq[1]-freq[0])))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--seconds",type=int,default=180)
    args=p.parse_args()
    root=args.source.resolve()
    got=subprocess.check_output(["git","-C",str(root),"rev-parse","HEAD"],text=True).strip()
    if got!=PIN or subprocess.check_output(["git","-C",str(root),"diff","HEAD","--name-only"]).strip():
        raise ValueError("Use the unchanged pinned ModelDB C++ checkout")
    args.out.mkdir(parents=True,exist_ok=True)
    binary=args.out/"costa_runner"
    wrapper=Path(__file__).with_name("costa_runner.cpp")
    subprocess.run(["g++","-std=c++11","-O3",f"-I{root}",str(wrapper),
                   str(root/"Cortical_Column.cpp"),str(root/"Thalamic_Column.cpp"),
                   "-o",str(binary)],check=True)
    results={"paper":"10.1371/journal.pcbi.1005022","source_commit":PIN,"sampling_hz":FS,
             "burn_s":20,"simulated_seconds_per_seed":args.seconds,"seeds":SEEDS,
             "modality":"continuous_EEG_proxy","spike_counts_created":False,
             "seed_deviation_from_MEX":"explicit std::srand seeds instead of wall time",
             "spectra":"Welch 10 s Hann, 50% overlap, peaks are band-limited",
             "estimator":"continuous spectral validation, no Sooter estimates",
             "reference":{"figure":"Fig. 4–5","slow_oscillation_hz":[.2,1.],
                          "fast_spindle_hz":[12,15],"validation_scope":"band signatures, not full ERP or spindle detector"},
             "results":{}}
    for state in ("N2","N3"):
        parameter_file=root/f"Figures/Data/Parameter_{state}.mat"
        mat=loadmat(parameter_file)
        pieces=[mat[key].ravel(order="F") for key in ("Param_Cortex","Param_Thalamus","Connectivity")]
        params=np.concatenate(pieces)
        if list(map(len,pieces))!=[3,2,4] or not np.isfinite(params).all():
            raise ValueError("Unexpected author parameter array")
        seeds=[]
        for seed in SEEDS:
            now=perf_counter()
            response=subprocess.run([str(binary.resolve()),str(seed),str(args.seconds)],
                        input=" ".join(map(str,params)),capture_output=True,text=True,timeout=200,check=True)
            data=np.fromstring(response.stdout,sep=" ").reshape(-1,2)
            if data.shape!=(args.seconds*int(FS),2) or not np.isfinite(data).all():
                raise ValueError(f"Invalid continuous output {state}/{seed}")
            np.save(args.out/f"{state}-{seed}.npy",data,allow_pickle=False)
            seeds.append(dict(seed=seed,wall_s=round(perf_counter()-now,3),
                              output_sha256=sha256(data.astype("<f8").tobytes()).hexdigest(),
                              cortical=summarize(data[:,0]),thalamic=summarize(data[:,1])))
        results["results"][state]=dict(parameter_sha256=sha256(parameter_file.read_bytes()).hexdigest(),
             parameters=dict(zip(("cortex","thalamus","connections"),[x.tolist() for x in pieces])),
             runs=seeds)
    rows=results["results"]
    comparisons=[]
    for n2,n3 in zip(rows["N2"]["runs"],rows["N3"]["runs"],strict=True):
        slow_ratio=n3["cortical"]["slow"]["band_power_mV2"]/n2["cortical"]["slow"]["band_power_mV2"]
        spindle_hz=n2["thalamic"]["spindle"]["peak_hz"]
        comparisons.append(dict(seed=n2["seed"],n3_to_n2_cortical_slow_power=slow_ratio,
                                n2_thalamic_spindle_candidate_hz=spindle_hz,
                                joint_direction=bool(slow_ratio>1 and 12<=spindle_hz<=15)))
    results["spectral_checks"]=comparisons
    results["gate"]="SMOKE_ONLY" if args.seconds<120 else (
        "PASS_SPECTRAL_SIGNATURE" if all(r["joint_direction"] for r in comparisons) else "FAIL_SPECTRAL_SIGNATURE")
    (args.out/"continuous_validation.json").write_text(json.dumps(results,indent=2,allow_nan=False)+"\n")
    print(results["gate"],comparisons,flush=True)


if __name__=="__main__":
    main()
