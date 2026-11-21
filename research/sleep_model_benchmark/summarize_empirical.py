"""Verify the duration-matched empirical run and aggregate windows→sessions→rats."""
import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
from statistics import mean, median

from prepare_empirical import SOURCE_COMMIT, SOURCE_SHA
from prepare_fits import write
from summarize import PIN, model_qc

RUN="44896254896"
COMMIT="29c31e1031726e055d2ea7aeff163dd799fb25b2"
ACTION="https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/"+RUN


def summarize(path,artifact,checksum):
    given=sha256((path/"actions.zip").read_bytes()).hexdigest()
    if given!=checksum.removeprefix("sha256:"):
        raise ValueError("GitHub artifact SHA-256 mismatch")
    plan,=json.loads((path/"inputs/fit_plan.json").read_text())
    select=json.loads((path/"inputs/empirical_selection.json").read_text())
    raw=(path/"fits/empirical/matlab_fits.json").read_bytes()
    all_fits=json.loads(raw)
    assert select["archive_commit"]==SOURCE_COMMIT and select["archive_sha256"]==SOURCE_SHA
    assert all_fits["toolkit_commit"]==PIN and "R2024b" in all_fits["matlab"]
    assert str(all_fits["run_id"])==RUN and all_fits["project_commit"]==COMMIT
    assert all_fits["input"]["sha256"]==plan["sha256"] and plan["analysis_spec"]==[.04,5,10]
    assert all_fits["settings"]["fit_method"]=="YuleWalker" and all_fits["settings"]["beta"]==2
    assert len(all_fits["records"])==plan["count"]==all_fits["input"]["records"]
    stations={(x["session"],x["state"],x["window_index"]):x for x in select["stationarity"]}
    rows=[]
    for item in all_fits["records"]:
        m=item["record"]
        assert (m["dandiset"],m["region"],m["context"])==("000041","mPFC","HOME")
        assert 20<=m["unit_count"]<=36 and item["fit"]["duration_s"]==240
        qc=stations.pop((m["session"],m["state"],m["window_index"]))
        fits=[model_qc(item["fit"],qc,False),model_qc(item["long_reach"],qc,False)]
        rows.append(dict(record=m,stationarity=qc,fits=fits))
    assert not stations
    aggregates=[]
    for order in (5,10):
        k={5:0,10:1}[order]
        for rule,field in (("numerical_only","numeric_ok"),("strict_QC","pass_screen")):
            animals=defaultdict(lambda:defaultdict(dict))
            sessions=defaultdict(list)
            for row in rows:
                m=row["record"]
                if row["fits"][k][field]:
                    sessions[m["subject"],m["session"],m["state"]].append(row["fits"][k]["d2"])
            for (animal,session,state),ds in sessions.items():
                animals[animal][state][session]=mean(ds)
            state_stats={}
            for state in ("Wake","NREM","REM"):
                values={animal:mean(d[state].values()) for animal,d in animals.items() if state in d}
                state_stats[state]=dict(animals=len(values),session_count=sum(s==state for _,_,s in sessions),
                    window_count=sum(len(ds) for (_,_,s),ds in sessions.items() if s==state),
                    d2_median=median(values.values()) if values else None,by_animal=values)
            nw={}
            for animal,d in animals.items():
                w,n=d.get("Wake",{}),d.get("NREM",{})
                pairing=set(n)&set(w)
                if pairing:
                    nw[animal]=mean(n[s]-w[s] for s in pairing)
            aggregates.append(dict(order=order,rule=rule,states=state_stats,animal_paired_NW=nw))
    return dict(provenance=dict(workflow=ACTION,commit=COMMIT,matlab=all_fits["matlab"],toolkit_commit=PIN,
                artifact=artifact,sha256=checksum,fit_sha256=sha256(raw).hexdigest(),input_sha256=plan["sha256"]),
                selection=select, windows=rows, aggregates=aggregates)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--artifact",type=int,required=True)
    p.add_argument("--sha256",required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    o=summarize(a.receipt,a.artifact,a.sha256)
    write(a.out,o)
    print(f"Verified {len(o['windows'])} duration-matched fits and the nested rat summaries.")
