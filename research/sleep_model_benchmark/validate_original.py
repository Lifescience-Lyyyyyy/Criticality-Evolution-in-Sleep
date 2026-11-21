"""Levenstein Fig. 4 summary reproduction with calibration/validation seeds."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

PARAMETERS = {
    "rat_frontal": dict(I=2.64, w=6.28, b=1.0, sigma=.25, theta=.05),
    "rat_CA1": dict(I=1.90, w=6.00, b=1.0, sigma=.25, theta=.05),
}
TARGETS = {
    "rat_frontal": dict(up_mean_s=[1.70, .92], down_mean_s=[.21, .05],
                       up_cv=[1.10, .27], down_cv=[.38, .06]),
    "rat_CA1": dict(up_mean_s=[.06, .005], down_mean_s=[2.00, .22],
                    up_cv=[.33, .04], down_cv=[1.30, .10]),
}
CALIBRATION_SEEDS = (711, 722, 733)
FORMAL_SEEDS = (1441, 1442, 1443)
BURN_AU = 3000
KEPT_AU = 60000


def simulate(executable, seed, params, substeps=10, kept=KEPT_AU, gate=None):
    args = [seed, kept+BURN_AU, substeps, *(params[k] for k in ("I","w","b","sigma","theta"))]
    if gate is not None:
        args.extend(gate)
    result = subprocess.run([str(executable.resolve()), *map(str,args)],
                            check=True, capture_output=True, text=True, timeout=40)
    rates = np.fromstring(result.stdout, sep="\n")
    if len(rates)!=kept+BURN_AU or not np.isfinite(rates).all() or (rates<0).any() or (rates>1).any():
        raise ValueError("ODE output has incorrect duration or invalid normalized rates")
    return rates[BURN_AU:]


def dwells(rates):
    """Hysteresis between the histogram trough and two modes; first/last bouts censored."""
    histogram, edges = np.histogram(rates, bins=40, range=(0,1))
    centers = (edges[:-1] + edges[1:])/2
    smoothed = gaussian_filter1d(histogram.astype(float), 1)
    peaks, properties = find_peaks(np.r_[0,smoothed,0], prominence=smoothed.max()*.005)
    peaks -= 1
    if len(peaks)<2:
        return {"status":"NO_TWO_MODES"}
    strongest = np.argsort(properties["prominences"])[-2:]
    low, high = sorted(peaks[strongest])
    trough = low + int(np.argmin(smoothed[low:high+1]))
    if trough in (low,high) or centers[high]-centers[low] < .20:
        return {"status":"NO_SEPARATED_MODES"}
    if smoothed[trough] > .75 * min(smoothed[low],smoothed[high]):
        return {"status":"NO_DEEP_TROUGH"}
    down, up = (centers[low]+centers[trough])/2, (centers[high]+centers[trough])/2
    above, below = rates > up, rates < down
    hits = np.flatnonzero(above | below)
    states = above[hits]
    switches = np.flatnonzero(states[1:] != states[:-1]) + 1
    times = hits[switches]
    duration = np.diff(times)
    was_up = states[switches[:-1]]
    durations = {k: duration[was_up == state] for k,state in (("up",True),("down",False))}
    if min(map(len,durations.values())) < 12:
        return {"status":"INSUFFICIENT_COMPLETE_DWELLS"}
    out = {"status":"OK", "up_threshold":float(up), "down_threshold":float(down),
           "up_occupancy_uncensored":float(durations["up"].sum()/duration.sum())}
    for k,v in durations.items():
        out.update({k+"_mean_au":float(v.mean()),k+"_cv":float(v.std(ddof=1)/v.mean()),
                    k+"_events":len(v),k+"_median_au":float(np.median(v))})
    return out


def evaluate(out, executable):
    output = {"paper":"https://doi.org/10.1038/s41467-019-10327-5",
              "source_commit":"989d16ad8fa157d1f0bda7a8b29199c1aa8f9292",
              "implementation":"independent published-equation reproduction (not author MATLAB simulation)",
              "simulator_sha256":sha256(Path(__file__).with_name("levenstein.cpp").read_bytes()).hexdigest(),
              "burn_au":BURN_AU,"kept_au":KEPT_AU,"sample_au":1,"integration_dt_au":.1,
              "calibration_seeds":CALIBRATION_SEEDS,"validation_seeds":FORMAL_SEEDS,
              "threshold_procedure":"histogram two modes + trough midpoint hysteresis; not author's dip test",
              "reference_target":"Fig. 4 and text, mean and SD over recordings",
              "results":{}}
    for case,params in PARAMETERS.items():
        target=TARGETS[case]
        calibration_event="down" if case=="rat_frontal" else "up"
        calibration=[dwells(simulate(executable,s,params)) for s in CALIBRATION_SEEDS]
        if any(r["status"]!="OK" for r in calibration):
            raise RuntimeError(f"Calibration runs did not supply sufficient complete bouts: {case}: {calibration}")
        seconds_per_au = round(target[calibration_event+"_mean_s"][0] /
                               np.mean([r[calibration_event+"_mean_au"] for r in calibration]), 4)
        checks=[]
        for seed in FORMAL_SEEDS:
            trace=simulate(executable,seed,params)
            row=dwells(trace)
            if row["status"]=="OK":
                for k in ("up","down"):
                    row[k+"_mean_s"]=row[k+"_mean_au"]*seconds_per_au
            checks.append(dict(row,seed=seed))
            np.save(out/f"{case}-{seed}.npy",trace,allow_pickle=False)
        means = ({k:float(np.mean([r[k] for r in checks])) for k in target}
                 if all(r["status"]=="OK" for r in checks) else {})
        differences = {k:(v-target[k][0])/target[k][1] for k,v in means.items()}
        tested = ("up_mean_s","up_cv","down_cv") if calibration_event=="down" else ("down_mean_s","up_cv","down_cv")
        passed = (.001 <= seconds_per_au <= .025 and
                  all(k in differences and abs(differences[k])<=2 for k in tested))
        finer=dwells(simulate(executable,FORMAL_SEEDS[0],params,substeps=20))
        output["results"][case]=dict(parameters=params,calibration=calibration,time_scale_s=seconds_per_au,
            calibration_metric=calibration_event+"_mean_s", reference=target,observed=means,
            deviation_in_reference_sd=differences,validation=checks,validated_metrics=tested,
            finer_dt_first_seed=finer, gate="PASS_SUMMARY" if passed else "FAIL_SUMMARY")
    return output


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    src=Path(__file__).with_name("levenstein.cpp")
    executable=args.out/"levenstein"
    subprocess.run(["g++","-O3","-std=c++17",str(src),"-o",str(executable)],check=True)
    result=evaluate(args.out,executable)
    (args.out/"original_validation.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    for case,r in result["results"].items():
        print(case,r["gate"],r["time_scale_s"],r["observed"],flush=True)


if __name__=="__main__":
    main()
