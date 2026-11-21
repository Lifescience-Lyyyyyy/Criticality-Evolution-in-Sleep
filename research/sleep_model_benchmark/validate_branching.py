"""Test reported González Fig. 6 ungated RQA and explicit gating interpretations."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from statistics import median
import subprocess

import numpy as np
from scipy.spatial.distance import cdist

from observation import spike_counts

SEEDS=(1441,1442,1443)
# Broad, preregistered pilot acceptance intervals due to Fig. 6 graphical approximations.
REFERENCE=dict(rr=.02,det=.20,lam=.40,div=.30,tt=2.5)
TOLERANCE=dict(rr=.012,det=.12,lam=.15,div=.15,tt=.75)
CONDITIONS={
    "paper_ungated":(1.0,"N",0,500,0),
    "recurrence_only":(.8,"N",0,500,0),
    "drive_gate_250ms":(1.0,"D",250,500,.50),
    "propagation_gate_250ms":(1.0,"P",250,500,.50),
    "recurrence_and_drive_gate":(.8,"D",250,500,.50),
    "drive_gate_4hz_125ms":(1.0,"D",125,250,.50),
}


def line_lengths(x):
    starts=np.flatnonzero(np.diff(np.r_[False,x,False]))
    return starts[1::2]-starts[::2]


def rqa(window, epsilon):
    r=cdist(window,window,"euclidean") <= epsilon
    # Drop the main diagonal from all numerators and denominators.
    np.fill_diagonal(r,False)
    points=int(r.sum())
    diag=[line_lengths(r.diagonal(k)) for k in range(1,len(r))]
    vert=[line_lengths(row) for row in r.T]
    ds=np.concatenate(diag)
    vs=np.concatenate(vert)
    long_diag=ds[ds>=2]
    long_vert=vs[vs>=2]
    return dict(rr=points/(len(r)*(len(r)-1)),det=2*long_diag.sum()/points if points else 0.,
                lam=long_vert.sum()/points if points else 0.,
                div=1/np.max(ds) if len(ds) else 0.,tt=long_vert.mean() if len(long_vert) else 0.)


def native_events(binary,seed,condition):
    cmd=[str(binary.resolve()),str(seed),"1010",*map(str,condition)]
    data=subprocess.run(cmd,capture_output=True,check=True,timeout=55).stdout
    words=np.frombuffer(data,dtype="<u8")[10000:]
    assert len(words)==1000000
    by_unit=[]
    for i in range(50):
        ticks=np.flatnonzero((words>>np.uint64(i))&1)
        by_unit.append(ticks/1000)
    return by_unit,words


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    binary=args.out/"branching"
    src=Path(__file__).with_name("branching.cpp")
    subprocess.run(["g++","-std=c++17","-O3",str(src),"-o",str(binary)],check=True)
    observations={}
    for seed in SEEDS:
        epsilon=None
        for mode,settings in CONDITIONS.items():
            events,words=native_events(binary,seed,settings)
            fifty=spike_counts([(0,1000,events)],.05).segments[0]
            forty=spike_counts([(0,1000,events)],.04).population()[0]
            if forty.sum()!=fifty.sum():
                raise ValueError("Population counts differ on the two exact partitions")
            if epsilon is None:
                epsilon=float(np.std(fifty.sum(1),ddof=1))
            windows=np.linspace(0,99,12,dtype=int)
            stats=[rqa(fifty[i*200:(i+1)*200],epsilon) for i in windows]
            effects={k:float(np.mean([s[k] for s in stats])) for k in REFERENCE}
            silence=line_lengths(forty==0)
            entry=dict(seed=seed,epsilon=epsilon,mean_RQA=effects,
                       rate_hz=float(forty.mean()/.04),
                       off_p95_s=float(np.quantile(silence,.95)*.04) if len(silence) else None,
                       observations_sha256=sha256(words.tobytes()).hexdigest())
            observations.setdefault(mode,[]).append(entry)
    aggregate={mode:{k:median(r["mean_RQA"][k] for r in rs) for k in REFERENCE}
               for mode,rs in observations.items()}
    distances={k:abs(aggregate["paper_ungated"][k]-v) for k,v in REFERENCE.items()}
    passed=all(distances[k]<=v for k,v in TOLERANCE.items())
    results=dict(source="10.1371/journal.pone.0290146 Fig. 6 / Methods",official_code=False,
                 status="PASS_UNGATED_RQA_PILOT" if passed else "FAIL_UNGATED_RQA_PILOT",
                 simulated_seconds=1000,burn_seconds=10,seeds=SEEDS,
                 implementation_assumptions=["undirected ER; no self edges", "independent edge successes conditional on active neighbors",
                 "Euclidean RQA; exclude line of identity; minimum line length 2; 12 prespecified 10 s windows/seeds",
                 "Method reports 250 ms and 4 Hz; attempted 250/500 and 125/250 ms as separate counterfactuals, neither official",
                 "D gates external inputs to 50% of cells, P gates outgoing propagation from the same cells"],
                 original_reference=REFERENCE,tolerance=TOLERANCE,absolute_deviation=distances,
                 conditions=CONDITIONS,seed_results=observations,aggregate=aggregate,
                 d2_status="NOT_RUN_native_branching_not_in_formal_model_table")
    (args.out/"branching_validation.json").write_text(json.dumps(results,indent=2,allow_nan=False)+"\n")
    print(results["status"],aggregate["paper_ungated"],flush=True)


if __name__=="__main__":
    main()
