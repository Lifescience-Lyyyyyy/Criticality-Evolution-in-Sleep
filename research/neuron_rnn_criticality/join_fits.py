"""Validate the GitHub artifact against its paired count corpus and retain fit warnings."""
import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from physiology import digest, write
from metrics import lag_cov

PIN="f26c782dca3b1c9ab18cc92c6449591c44e0d507"


def valid(f):
    return (f is not None and f["d2"] is not None and f["d2"]>=-1e-10
            and f["exit_status"]==0 and f["solver_flag"] is not None and f["solver_flag"]>0
            and f["min_ar_root"] is not None and f["min_ar_root"]>1
            and f["innovation_variance"] is not None and f["innovation_variance"]>0)


def info(f):
    passed=valid(f)
    if f is None:
        return {"numerical_pass":False,"d2":None}
    acf=[abs(v) for v in f["residual_acf"] if v is not None]
    kurt=f["residual_excess_kurtosis"]
    peak=max(acf) if acf else None
    return dict(numerical_pass=passed,d2=f["d2"] if passed else None,d2_raw=f["d2"],
                ar_order=f["ar_order"],exit_status=f["exit_status"],solver_flag=f["solver_flag"],
                min_ar_root=f["min_ar_root"],innovation_variance=f["innovation_variance"],
                max_residual_acf=peak,residual_excess_kurtosis=kurt,
                residual_abs_tail_over_3sd=f["residual_abs_tail_over_3sd"],acf_rmse=f["acf_rmse"],
                flag_residual_structure=peak>.10 if peak is not None else None,
                flag_excess_kurtosis=abs(kurt)>1 if kurt is not None else None)


def read_artifact(path,receipt):
    assert digest(path)==receipt["digest"].removeprefix("sha256:")
    with ZipFile(path) as z:
        assert not z.testzip()
        def read(name):
            member,=[n for n in z.namelist() if Path(n).name==name]
            return json.loads(z.read(member))
        fit=read("matlab_fits.json")
        plan,=read("plan.json")
    assert str(fit["run_id"])==str(receipt["id"]) and fit["project_commit"]==receipt["commit"]
    assert fit["toolkit_commit"]==PIN and "R2024b" in fit["matlab"]
    assert plan["input_sha256"]==fit["input"]["sha256"]
    settings=fit["settings"]
    assert [settings[k] for k in ("bin_seconds","ar_order","long_ar_order")]==plan["spec"]==[.04,5,10]
    assert settings["beta"]==2 and settings["fit_method"]=="YuleWalker"
    raw=fit["records"]
    records=raw if isinstance(raw,list) else [raw]
    assert len(records)==plan["records"]
    return records,fit,plan


def verify_fit_input(record,blocks):
    population=np.concatenate(blocks)
    expected_lengths=[len(block) for block in blocks]
    covariance=lag_cov(blocks,20,local=False)
    acf=covariance[1:]/covariance[0] if covariance[0]>0 else np.full(20,np.nan)
    for tag,order in (("fit",5),("long_reach",10)):
        fit=record[tag]
        assert fit["ar_order"]==order, (record["record"]["id"],tag)
        np.testing.assert_array_equal(np.atleast_1d(fit["segment_bins"]),expected_lengths)
        assert abs(fit["duration_s"]-len(population)*.04)<1e-8
        assert np.isclose(fit["rate_hz"],population.mean()/.04,rtol=1e-10)
        assert np.isclose(fit["count_variance"],population.var(),rtol=1e-10)
        np.testing.assert_allclose(np.asarray(fit["observed_acf"],dtype=float),acf,atol=1e-10,rtol=1e-10)


def paired(records):
    by=defaultdict(dict)
    bugnon=defaultdict(dict)
    for row in records:
        key=row["id"]
        if row["category"]=="original":
            by[row["connectivity_seed"]][row["context"]]=row["fit"]["d2"]
        elif row["category"]=="within_episode_shift":
            by[(row["connectivity_seed"],row["context"])].setdefault("shift",[]).append(row["fit"]["d2"])
        elif row["category"] in ("bugnon_original","bugnon_subsample"):
            bugnon[(row["n_units"],row.get("sample_seed"))][row["context"]]=row
        elif row["category"]=="empirical":
            by[(row["animal"],row["session"])][row["state"]]=row["fit"]["d2"]
    model=[]
    surrogate=[]
    observed=defaultdict(list)
    for key,v in by.items():
        if isinstance(key,int):
            a,b=v.get("figure7"),v.get("figure7_supp2")
            if a is not None and b is not None:
                model.append(dict(connectivity_seed=key,figure7=a,figure7_supp2=b,delta_supp2_minus_figure7=b-a))
        elif "shift" in v:
            topo,condition=key
            base=by[topo][condition]
            q=[x for x in v["shift"] if x is not None]
            if base is not None and q:
                surrogate.append(dict(connectivity_seed=topo,condition=condition,original=base,
                    shifted_mean=float(np.mean(q)),delta_shift_minus_original=float(np.mean(q))-base,replicates=q))
        else:
            a,b=v.get("Wake"),v.get("NREM")
            if a is not None and b is not None:
                observed[key[0]].append(b-a)
    result=dict(author_conditions=model,within_episode_shift=surrogate,
                empirical_nrem_minus_wake={k:float(np.mean(y)) for k,y in observed.items()})
    if bugnon:
        pairs=[]
        for (n,seed),conditions in bugnon.items():
            wake,sleep=conditions["wake"],conditions["sleep"]
            for field in ("connectivity_seed","noise_seeds","sampled_gids","duration_seconds"):
                if wake.get(field)!=sleep.get(field):
                    raise ValueError(f"Bugnon pair mismatch for {field} at N={n}, subset seed {seed}")
            deltas={}
            for key in ("fit","ar10"):
                w,s=wake[key]["d2"],sleep[key]["d2"]
                deltas[key]=s-w if w is not None and s is not None else None
            pairs.append(dict(n_units=n,sample_seed=seed,delta_sleep_minus_wake=deltas))
        result["bugnon_sleep_minus_wake"]=pairs
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("bundle",type=Path)
    p.add_argument("zip",type=Path)
    p.add_argument("receipt",type=Path)
    p.add_argument("out",type=Path)
    p.add_argument("--prefix",default="baseline")
    a=p.parse_args()
    receipt=json.loads(a.receipt.read_text())
    runs,full,plan=read_artifact(a.zip,receipt)
    with gzip.open(a.bundle,"rt",encoding="utf-8") as f:
        orig={r["id"]:r for r in json.load(f)["cases"]}
    records=[]
    for r in runs:
        m=r["record"]
        z=orig.pop(m["id"])
        assert m=={k:v for k,v in z.items() if k!="counts"}
        verify_fit_input(r,z["counts"])
        records.append(dict(**m,fit=info(r["fit"]),ar10=info(r["long_reach"])))
    assert not orig
    a.out.mkdir(exist_ok=True,parents=True)
    write(a.out/f"{a.prefix}_fits.json",dict(records=records,paired=paired(records)))
    rc=dict(receipt,valid_originals=sum(r["category"] in ("original","bugnon_original")
                                        and r["fit"]["numerical_pass"] for r in records),
            numerical_primary_pass=sum(r["fit"]["numerical_pass"] for r in records),n_cases=len(records),
            numerical_ar10_pass=sum(r["ar10"]["numerical_pass"] for r in records),
            toolkit_commit=PIN,matlab=full["matlab"],matlab_count_sha256=plan["input_sha256"],
            paired_input_sha256=digest(a.bundle))
    write(a.out/f"{a.prefix}_matlab_receipt.json",rc)
    data=json.dumps(full,separators=(",",":")).encode()
    (a.out/f"{a.prefix}_full_matlab.json.gz").write_bytes(gzip.compress(data,mtime=0))
    print(json.dumps(rc),flush=True)


if __name__=="__main__":
    main()
