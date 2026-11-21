"""Verify MATLAB receipts and produce small numerical summaries."""
import argparse
from collections import defaultdict
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from statistics import median

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit

from prepare_fits import SPECS, write
from validate_original import FORMAL_SEEDS, PARAMETERS

PIN="f26c782dca3b1c9ab18cc92c6449591c44e0d507"
RUN="44895142233"
COMMIT="18ea304b464eda46137bb955dfe88c47568f9f84"
URL="https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/" + RUN
QC=dict(residual_acf_fail=.15, absolute_residual_excess_fail=3., gaussian_warning_abs_excess=1.)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def model_qc(fit, stationarity, forced):
    p=fit["ar_order"]
    d2, flag, root, innovation=(fit[k] for k in ("d2","solver_flag","min_ar_root","innovation_variance"))
    ok=(all(v is not None and isfinite(v) for v in (d2,root,innovation))
        and d2>=-1e-10 and root>1 and innovation>0 and fit["exit_status"]==0
        and (p==1 or (flag is not None and flag>0)))
    residual=fit["residual_acf"]
    residual_abs=max(map(abs,residual)) if residual and all(v is not None for v in residual) else None
    kurtosis=fit["residual_excess_kurtosis"]
    warnings=[]
    failures=[]
    if forced:
        failures.append("PERIODICALLY_FORCED_CYCLOSTATIONARY")
    if not stationarity["pass_"]:
        failures.append("BLOCK_DRIFT_SCREEN_FAILED")
    if (residual_abs is None or kurtosis is None or
        residual_abs>QC["residual_acf_fail"] or abs(kurtosis)>QC["absolute_residual_excess_fail"]):
        failures.append("RESIDUAL_SCREEN_FAILED")
    elif abs(kurtosis)>QC["gaussian_warning_abs_excess"]:
        warnings.append("RESIDUAL_EXCESS_KURTOSIS")
    if not ok:
        failures.append("NUMERICAL_FIT_FAILED")
    return dict(order=p, d2=d2, exit_status=fit["exit_status"],solver_flag=flag,
                solution_type="p=beta/2 explicit integral" if p==1 else "fmincon SQP",
                root_min=root,innovation_variance=innovation,rate_hz=fit["rate_hz"],
                fano=fit["fano"],zero_fraction=fit["zero_fraction"],
                max_abs_residual_acf=residual_abs,residual_excess_kurtosis=kurtosis,
                residual_abs_tail_over_3sd=fit["residual_abs_tail_over_3sd"],
                acf_rmse=fit["acf_rmse"],numeric_ok=ok,pass_screen=not failures,
                warnings=warnings,failures=failures)


def receipt(path, artifact_id, artifact_sha256):
    if digest(path/"actions.zip")!=artifact_sha256.removeprefix("sha256:"):
        raise ValueError("GitHub artifact SHA-256 mismatch")
    grid=json.loads((path/"inputs/fit_plan.json").read_text())
    assert {r["name"] for r in grid}==set(SPECS) and len(grid)==len(SPECS)
    original=json.loads((path/"original/original_validation.json").read_text())
    assert tuple(original["validation_seeds"])==FORMAL_SEEDS
    assert original["results"]["rat_frontal"]["gate"]=="PASS_SUMMARY"
    obs=json.loads((path/"inputs/observation_qc.json").read_text())
    remaining={row["record"]["key"]:row for row in obs["records"]}
    assert len(remaining)==len(obs["records"])
    summary=[]
    verified={}
    for entry in grid:
        name=entry["name"]
        p=path/"fits"/name/"matlab_fits.json"
        fit=json.loads(p.read_text())
        spec=SPECS[name]
        assert fit["toolkit_commit"]==PIN and str(fit["run_id"])==RUN and "R2024b" in fit["matlab"]
        assert fit["project_commit"]==COMMIT and fit["input"]["sha256"]==entry["sha256"]
        assert [fit["settings"][s] for s in ("bin_seconds","ar_order","long_ar_order")] == list(spec.values())
        assert fit["settings"]["fit_method"]=="YuleWalker" and fit["settings"]["beta"]==2
        assert len(fit["records"])==entry["count"]==fit["input"]["records"]
        verified[name]=dict(fit_sha256=digest(p),input_sha256=entry["sha256"],
                           records=entry["count"],wall_s=fit["elapsed_seconds"])
        for row in fit["records"]:
            m=row["record"]
            checked=remaining.pop(m["key"])
            assert checked["record"]==m and m["spec"]==name and m["seed"] in FORMAL_SEEDS
            assert row["fit"]["duration_s"]>=120 and row["fit"]["ar_order"]==spec["order"]
            assert row["long_reach"]["ar_order"]==spec["alternate"]
            hard_gate=(checked["gate_schedule"] is not None)
            primary=model_qc(row["fit"],checked["stationarity"],hard_gate)
            alt=model_qc(row["long_reach"],checked["stationarity"],hard_gate)
            summary.append(checked | {"fits":[primary,alt]})
    assert not remaining
    rows=sorted(summary,key=lambda r:r["record"]["key"])
    groups=defaultdict(list)
    anchor={r["record"]["seed"]:r for r in rows if r["record"]["context"]=="original_poisson"
            and r["record"]["spec"]=="primary"}
    assert set(anchor)==set(FORMAL_SEEDS)
    for row in rows:
        m=row["record"]
        groups[m["spec"],m["context"]].append(row)
        if m["spec"]=="primary":
            ref=anchor[m["seed"]]
            row["delta_reference_d2"]=[fit["d2"]-base["d2"] if fit["d2"] is not None and base["d2"] is not None else None
                                      for fit,base in zip(row["fits"],ref["fits"],strict=True)]
            row["count_acf_distance_to_reference"]=float(np.sqrt(np.mean(
                (np.array(row["count_acf"])-np.array(ref["count_acf"]))**2)))
    cells=[]
    for (spec,context),xs in sorted(groups.items()):
        assert {x["record"]["seed"] for x in xs}==set(FORMAL_SEEDS)
        phys=[x["physiology"] for x in xs if x["physiology"] and x["physiology"]["status"]=="OK"]
        metrics={"summary_pass_count":sum(all(abs(z)<=2 for z in x["z_vs_published"].values()) for x in phys),
                 "detected":len(phys),
                 "median_z_rms":median(np.sqrt(np.mean(np.array(list(x["z_vs_published"].values()))**2))
                                      for x in phys) if phys else None}
        cells.append(dict(spec=spec,context=context,seed_count=len(xs),physiology=metrics,
                     fit_stats=[{"order":group[0]["order"],
                         "duration_history_ms":1000*SPECS[spec]["bin"]*group[0]["order"],
                         "d2_numerical_median":median(f["d2"] for f in group if f["numeric_ok"])
                                               if any(f["numeric_ok"] for f in group) else None,
                         "d2_screened_median":median(f["d2"] for f in group if f["pass_screen"])
                                               if any(f["pass_screen"] for f in group) else None,
                         "screened_n":sum(f["pass_screen"] for f in group)} for group in
                         zip(*(x["fits"] for x in xs),strict=True)]))
    return dict(provenance=dict(workflow=URL,commit=COMMIT,toolkit_commit=PIN,artifact=artifact_id,
                 artifact_sha256=artifact_sha256,matlab=fit["matlab"],batches=verified,
                 observation_sha256=digest(path/"inputs/observation_qc.json")),
                original=original,empirical=obs["empirical"],rate_scale_hz=obs["rate_scale_hz"],
                fit_qc=QC, rows=rows, groups=cells)


def pulse_response():
    """Compare two intervention channels with identical unperturbed r–a equations."""
    par=PARAMETERS["rat_frontal"]
    I,w,b=(par[k] for k in ("I","w","b"))
    au_s=.004
    def equilibrium(r):
        a=expit(15*(r-.5))
        return expit(w*r-b*a+I-5)-r
    steady=brentq(equilibrium,.65,.99)
    adaptation=expit(15*(steady-.5))
    grid=np.arange(0,150.0+.05,.1)
    replays=[]
    for pulse_au in (5,25):
        for effective in (.25,.75):
            for mode in ("external_drive","recurrent_gain"):
                def field(t,state,forced):
                    r,a=state
                    increase=(effective if mode=="external_drive" else effective*r/steady) if forced else 0
                    return expit(w*r-b*a+I-5+increase)-r,(expit(15*(r-.5))-a)/25
                before=grid[grid<=pulse_au]
                after=grid[grid>=pulse_au]
                first=solve_ivp(lambda t,y:field(t,y,True),(0,pulse_au),(steady,adaptation),
                                t_eval=before,rtol=1e-9,atol=1e-11)
                last=solve_ivp(lambda t,y:field(t,y,False),(pulse_au,grid[-1]),first.y[:,-1],
                               t_eval=after,rtol=1e-9,atol=1e-11)
                r=np.r_[first.y[0,:-1],last.y[0]]-steady
                peak=np.max(np.abs(r))
                final=np.flatnonzero(np.abs(r[len(before)-1:])>.1*peak)
                recovery=float((after[final[-1]]-pulse_au)*au_s) if len(final) else 0.
                replays.append(dict(mode=mode,initial_effective_drive=effective,pulse_s=pulse_au*au_s,
                                   peak_abs_delta=float(peak),last_above_10pct_s=recovery,
                                   absolute_auc_rate_s=float(np.trapezoid(np.abs(r),grid*au_s))))
    return dict(steady_r=steady,steady_adaptation=adaptation,matched_initial_effective_drive=True,
                mode="deterministic fixed-point perturbation, OU disabled only for this causal check",
                perturbations=replays)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--artifact",type=int,required=True)
    p.add_argument("--sha256",required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=receipt(a.receipt,a.artifact,a.sha256)
    a.out.mkdir(parents=True,exist_ok=True)
    original=result.pop("original")
    write(a.out/"original_validation.json",original)
    write(a.out/"matlab_benchmark.json",result)
    write(a.out/"pulse_response.json",pulse_response())
    print(f"Verified {len(result['rows'])} MATLAB count inputs and {len(result['groups'])} strata.")


if __name__=="__main__":
    main()
