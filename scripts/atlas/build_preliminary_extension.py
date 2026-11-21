"""Prepare the frozen Phase 11.1 unit and temporal spot checks."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZipFile

from build_atlas_data import ROOT, read, write
from build_brainstem_qc import brainstem_tasks, digest
from build_pairwise_fixture import OUT, spans_for
from build_qc_fixture import cells, count_entries, save_prepared
from build_audit_results import covered, spans, sorted_intervals
from recheck_preliminary_sources import selected_reference
from audit_crcns_sleep import RemoteZip
from prepare_phase10_brainstem import recording_dtype, recording_meta


def check_brainstem_clock(case,audit,brain,archive,inventory,rows):
    """Re-read only the 64 KiB prefixes of the independent EEG and EMG members."""
    source,=archive["files"]
    remote=RemoteZip(source["links"]["self"],source["size"])
    obs={}
    for signal in ("EEG","EMG"):
        wanted=f"lozano_hong_jin_2026/{audit['SL_PATH']}/{signal}.mat"
        member=next(x for x in inventory if x["name"]==wanted)
        n,h=recording_dtype(remote,member)
        assert n==audit["clock_partner"][signal.lower()+"_samples"]
        assert h==audit["clock_partner"][signal.lower()+"_header_sha256"]
        obs[signal]={"samples":int(n),"published_member":wanted,"header_sha256":h}
    info,_,score,stages,raw_codes,_,cuts,_,duration=recording_meta({**audit,"root":str(brain)})
    assert digest(score)==audit["scoring_sha256"] and info.get("laser_used","").lower()!="yes"
    for row in rows:
        m=row["record"]
        allowed=sorted_intervals(stages[m["state"]])
        for a,b in spans(m["intervals_s"]):
            assert 0<=a<b<=duration and covered(a+1e-9,b-1e-9,allowed), (case["session"],m["state"],a,b)
    return {**{k:case[k] for k in ("dandiset","subject","session","region","asset_id")},
            "author_scores":str(audit["state_path"]),"scoring_sha256":digest(score),"score_codes":dict(raw_codes),
            "independent_clocks":obs,"physical_duration_s":duration,"shared_accepted_bound_s":audit["clock"]["seconds"],
            "excluded_author_ranges":cuts,"condition":audit["cohort"],"intervention_file_listed":audit["matching_laser_asset"],
            "out_of_bounds_or_crossed_selected_intervals":0,"source_range_bytes":remote.fetched_bytes}


def split_rows(archives,targets):
    expected={(c["asset_id"],c["region"],s) for c in targets for s in ("Wake","NREM")}
    ref=selected_reference()
    result,hashes=[],{}
    for path in archives:
        with ZipFile(path) as z:
            receipt=json.loads(z.read("cohort/phase10_input.json"))
            parts=receipt.get("parts") or [{"name":receipt["fixture"],"sha256":receipt["sha256"]}]
            for p in parts:
                raw=z.read("cohort/"+p["name"])
                assert sha256(raw).hexdigest()==p["sha256"]
                hashes[f"{path.stem}/{p['name']}"]=p["sha256"]
                for m,segments in count_entries(raw):
                    key=(m["asset_id"],m["region"],m["state"])
                    if key not in expected or m["variant"]!="reference" or m["context"]!="REST":
                        continue
                    a=next(x for x in ref[key[:2]] if x["record"]["state"]==key[-1])
                    assert a["record"]==m
                    fit=a["fit"]
                    assert sum(int(s.sum()) for s in segments)==round(fit["duration_s"]*fit["rate_hz"])
                    saved=spans_for(m,segments)
                    i=max(range(len(segments)),key=lambda i: (len(segments[i]),-saved[i][0]))
                    first,last=saved[i]
                    assert last-first>=6000
                    for label,offset in (("first",first),("last",last-3000)):
                        data=segments[i][offset-first:offset-first+3000]
                        record={**m,"variant":"within_state_"+label,
                                "intervals_s":[[offset*.04,(offset+3000)*.04]],"parent_segment_s":[first*.04,last*.04]}
                        result.append((record,cells([data])))
                    expected.remove(key)
    assert not expected,expected
    return result,hashes


def compare_nwb_context(rows):
    old=read("manifests/sleep_sources.json")["assets"]
    for row in rows:
        if row["dandiset"]=="000939":
            acquisition=next(x for x in row["published_clock_observations"] if x["path"]=="/acquisition/ElectricalSeries")
            assert [acquisition["start"],acquisition["end"]]==row["independent_clock"]["seconds"]
            row["acquisition_and_lfp_support_agree"]=True
        if row["dandiset"]=="000552":
            saved=next(x for x in old if x["asset_id"]==row["asset_id"])
            assert "home cage" in saved["description"].lower()
            assert set(saved["digest"].values())<=set(row["digest"].values())
            row["published_session_description"]=saved["description"]
            row["published_intervention_table"]=saved.get("opto_excluded_table")
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("nwb",type=Path)
    p.add_argument("phase10",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    chosen=read(OUT+"sorting_checks/selection.json")["cases"]
    refs=selected_reference()
    brainstem=[{**x,"case":x["session"],"regions":[x["region"]]} for x in chosen if not x["dandiset"].isdecimal()]
    rows={r["asset_id"]:r for r in read("results/phase10/audit/brainstem.json")}
    extra,receipt=brainstem_tasks(a.phase10/"brainstem",brainstem,rows,
                                  target_cases=tuple(c["case"] for c in brainstem),percentiles=(50,))
    inventory=json.loads((a.phase10/"meta/zenodo_zip_inventory.json").read_text())
    archive=json.loads((a.phase10/"meta/zenodo-19462601.json").read_text())
    # The previously untested case which also has an eligible within-state split.
    target=next(x for x in brainstem if x.get("temporal_split"))
    clock=check_brainstem_clock(target,rows[target["asset_id"]],a.phase10/"brainstem",archive,inventory,
                                refs[target["asset_id"],target["region"]])
    temporal=[c for c in chosen if c.get("temporal_split")]
    split,source_sha=split_rows([a.phase10/f"phase10-{name}-cohort.zip" for name in "abc"]+
                                [a.phase10/"phase10-c-rescue.zip"],temporal)
    assert len(split)==4*len(temporal) and len(temporal)==3
    raw=(a.nwb/"nwb_counts.mat").read_bytes()
    receipt_nwb=json.loads((a.nwb/"nwb_counts.json").read_text())
    assert sha256(raw).hexdigest()==receipt_nwb["sha256"]
    nwbrecs=[(m,cells(x)) for m,x in count_entries(raw)]
    quality=json.loads((a.nwb/"nwb_sorting.json").read_text())
    nwbtiming=json.loads((a.nwb/"nwb_timing.json").read_text())
    assert len(quality)==len(nwbtiming)==7
    sorting=ROOT/(OUT+"sorting_checks")
    write(sorting/"source_metrics.json",{"NWB":quality,"brainstem":receipt})
    path=ROOT/(OUT+"temporal_checks/source_clock.json")
    write(path,compare_nwb_context(nwbtiming)+[clock])
    source_sha.update(nwb_fixture=receipt_nwb["sha256"],
                      source_metrics=sha256((sorting/"source_metrics.json").read_bytes()).hexdigest(),
                      source_clock=sha256(path.read_bytes()).hexdigest())
    write(sorting/"colab_receipt.json",dict(runtime="CPU",records=receipt_nwb["records"],
               unit_count_inputs_sha256=receipt_nwb["sha256"],source_cases=len(nwbtiming),
               state_receipt_sha256=source_sha["source_clock"]))
    data=nwbrecs+extra+split
    a.out.mkdir(parents=True,exist_ok=True)
    out=a.out.resolve()/"validation_input.mat"
    info=save_prepared(data,out,"states",sorted({m["asset_id"] for m,_ in data}),source_sha256=source_sha)
    out.with_suffix(".sha256").write_text(f"{info['sha256']}  {out.relative_to(ROOT)}\n")
    print(len(data), "Matlab rows:", len(nwbrecs), "NWB sort,",len(extra),"PAG sort,",len(split),"within-state")


if __name__=="__main__":
    main()
