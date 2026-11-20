"""Cross-check 1 ms brainstem spike counts against the author's companion exports."""
import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
from struct import unpack_from
import zlib

import numpy as np
from scipy.io import whosmat

from audit_crcns_sleep import RemoteZip
from audit_sleep_dandi import subtract
from fetch_phase10_brainstem import download
from prepare_phase10 import finish, run_stage, span

SOURCE = "zenodo:19462601"
STATES = {1:"REM",2:"Wake",3:"NREM"}
NON_REGIONS = {"root","void","fiber tracts","tspc","ml","mcp","scp","icp","ll","sct","poc"}


def parse_config(path):
    rows={}
    for line in path.read_text().splitlines():
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        key, value = (v.strip() for v in line.split(":",1))
        if key=="MOUSE":
            rows[value]=current={"mouse":value,"cohort":path.parent.name}
        else:
            current[key]=value.removeprefix("<BASE_PATH>/")
    return list(rows.values())


def nucleus_counts(path):
    channels=json.loads((path/"channel_locations.json").read_text())
    with (path/"cluster_info.TSV").open(newline="") as fd:
        table=list(csv.DictReader(fd,delimiter="\t"))
    units={}
    for unit in table:
        kind=unit.get("group") or unit.get("KSLabel")
        chan=channels.get("channel_"+unit["ch"])
        if kind=="good" and chan and chan.get("brain_region"):
            key=unit.get("cluster_id",unit.get("id"))+"_good"
            if key in units:
                raise ValueError("Duplicate author cluster ID")
            units[key]=(chan["brain_region"],chan["brain_region_id"],int(unit["n_spikes"]))
    return units


def recording_dtype(remote, item):
    pos=item["header_offset"]
    raw=remote.fetch(pos,min(pos+65536,remote.size))
    name=item["name"]
    n,extra=unpack_from("<HH",raw,26)
    signature,_,_,method=unpack_from("<4sHHH",raw)
    if signature!=b"PK\x03\x04" or raw[30:30+n].decode()!=name:
        raise ValueError("EEG header does not match ZIP inventory")
    prefix=raw[30+n+extra:]
    if method==8:
        prefix=zlib.decompressobj(-15).decompress(prefix,1024**2)
    fields=whosmat(BytesIO(prefix))
    (key,shape,_),=fields
    if key!=Path(name).stem or len(shape)!=2 or min(shape)!=1 or max(shape)<1000:
        raise ValueError("Unexpected companion signal dimension")
    return int(max(shape)),sha256(prefix).hexdigest()


def state_index(path, seconds, cut):
    code, status=np.loadtxt(path,dtype=int,ndmin=2).T
    if abs(len(code)*2.5-seconds)>2.5+1e-9:
        raise ValueError(f"Hypnogram extent differs from recorded EEG by >2.5 s: {len(code)*2.5}, {seconds}")
    boundaries=np.r_[0,np.flatnonzero(code[:-1]!=code[1:])+1,len(code)]
    spans=defaultdict(list)
    for a,b in zip(boundaries[:-1],boundaries[1:],strict=True):
        if code[a] in STATES and a*2.5<seconds and b*2.5<=seconds+1e-9:
            spans[STATES[code[a]]].append((a*2.5,b*2.5))
        elif code[a] in STATES and b==len(code) and (b-1)*2.5>a*2.5:
            spans[STATES[code[a]]].append((a*2.5,(b-1)*2.5))
    cuts=[]
    for entry in cut.split(";") if cut else []:
        a,b=entry.split("-")
        cuts.append((float(a),seconds if b=="$" else float(b)))
    counts=Counter(code.tolist())
    return {s:subtract(spans[s],cuts) for s in STATES.values()},counts,Counter(status.tolist()),cuts


def recording_meta(asset):
    root=Path(asset["root"])
    np_path=root/asset["NP_PATH"]
    sl_path=root/asset["SL_PATH"]
    info={k.strip():v.strip() for k,v in (line.split(":",1) for line in (sl_path/"info.txt").read_text().splitlines() if ":" in line)}
    units=nucleus_counts(np_path)
    score=sl_path/("remidx_"+sl_path.name+".txt")
    partner=asset.get("clock_partner")
    sr=float(info["SR"])
    observed=info.get("actual_duration")
    raw=np.loadtxt(score,dtype=int,ndmin=2)
    duration=partner["eeg_samples"]/sr if partner else len(raw)*2.5
    states,counts,notes,cuts=state_index(score,duration,asset.get("KCUT",""))
    return info,units,score,states,counts,notes,cuts,observed,duration


def inspect(asset,cfg):
    root=Path(asset["root"])
    row={k:v for k,v in asset.items() if k!="root"}
    row.update(dandiset=SOURCE,session=Path(asset["SL_PATH"]).name,subject=asset["subject"],
               state_path=asset["SL_PATH"]+"/remidx_"+Path(asset["SL_PATH"]).name+".txt",
               version="Ver2")
    info,units,score,states,code,notes,cuts,observed,duration=recording_meta(asset)
    regions=Counter(r for r,_,_ in units.values())
    row.update(units_total=len(units),source_regions=dict(regions),score_codes=dict(code),score_notes=dict(notes),
               annotated_spans={s:len(a) for s,a in states.items()},
               annotated_seconds={s:span(a) for s,a in states.items()},
               quality_cut=cuts,info_recorded_duration=observed,laser_used_info=info.get("laser_used"),
               source_package_md5="b4283da46ef2989eab44ccd20a2d77cc")
    if asset["cohort"]!="WT":
        reason="Intervention cohort. Primary analysis uses the separate unperturbed WT mice."
        return {**row,"gate":"PARTIAL","reason":reason},[]
    partner=asset["clock_partner"]
    if float(info["SR"])!=1000 or not partner or partner["eeg_samples"]!=partner["emg_samples"]:
        return {**row,"gate":"BLOCKED_METADATA","reason":"No concordant 1000 Hz EEG/EMG recording clock"},[]
    size=partner["eeg_samples"]
    row["clock"]={"source":"1k_train verified against EEG, EMG and 2.5 s author firing-rate export",
                  "seconds":[0.,size/1000],"rate_hz":1000,"clock_partner":partner}
    highres=root/asset["NP_PATH"]/"1k_train.npz"
    coarse_path=root/asset["NP_PATH"]/"traind.csv"
    good=defaultdict(dict)
    labels=[]
    with coarse_path.open(newline="") as fd:
        reader=csv.reader(fd)
        columns=next(reader)
        coarse=np.array([list(map(float,r)) for r in reader])
    if len(coarse)!=size//2500:
        return {**row,"gate":"BLOCKED_METADATA","reason":"2.5 s export and 1 ms physical clock differ"},[]
    index={u:i for i,u in enumerate(columns)}
    excluded={r.strip(" '\"") for r in asset.get("EXCLUDE","").split(",")}
    with np.load(highres,allow_pickle=False) as arrays:
        for u,(region,annotation,spike_count) in units.items():
            qc={"id":u,"region":region,"ontology_id":annotation,"sorter_spikes":spike_count,"eligible":False}
            labels.append(qc)
            if u in excluded or region.lower() in NON_REGIONS or u not in arrays.files or u not in index:
                qc["exclusion"]="Author exclusion, non-neural mapping, or missing export column"
                continue
            signal=arrays[u]
            if signal.shape!=(size,) or signal.dtype.kind not in "bui" or np.any(signal>1) or np.any(signal<0):
                qc["exclusion"]="Unexpected 1 ms vector shape or representation"
                continue
            counts=signal[:len(coarse)*2500].reshape(-1,2500).sum(axis=1)
            matches=float(np.mean(np.isclose(counts,coarse[:,index[u]]*2.5,atol=1e-7,rtol=0)))
            t=np.flatnonzero(signal)*.001
            qc.update(spikes=len(t),exact_coarse_bin_fraction=matches)
            if matches<1 or len(t)<2:
                qc["exclusion"]="Fine counts do not exactly reproduce companion export"
                continue
            good[region][u]=t
            qc.update(eligible=True,first=float(t[0]),last=float(t[-1]),
                      refractory=float(np.mean(np.diff(t)<cfg["refractory_seconds"])))
    def checksum(path):
        with path.open("rb") as stream:
            from hashlib import file_digest
            return file_digest(stream,"sha256").hexdigest()
    row.update(digest={"extracted_sha256":checksum(highres)},companion_sha256=checksum(coarse_path),
               scoring_sha256=checksum(score),units=labels,units_by_region={r:len(u) for r,u in good.items()})
    return finish(row,{"REST":states},good,cfg)


def catalog(metadata, inventory, root, scope):
    source,=metadata["files"]
    entries={r["name"]:r for r in inventory}
    remote=RemoteZip(source["links"]["self"],source["size"])
    base="lozano_hong_jin_2026/"
    candidates=[]
    cohorts=("WT",) if scope!="all" else ("WT","Opto","eYFP","vGlut","mPFC-LH")
    configs=[base+f"{c}/{c}_config.txt" for c in cohorts]
    download(metadata,inventory,root,configs)
    for c in cohorts:
        for asset in sorted(parse_config(root/c/(c+"_config.txt")),key=lambda r:r["mouse"]):
            if scope=="first" and asset["mouse"]!="DL108":
                continue
            np_path=asset["NP_PATH"]
            sl_path=asset["SL_PATH"]
            mouse=asset["mouse"]
            identity=(mouse.removesuffix("Red").removesuffix("Green") if c=="mPFC-LH"
                      else mouse.removesuffix("b") if c=="vGlut" else mouse)
            names=[base+f"{np_path}/{file}" for file in ("cluster_info.TSV","channel_locations.json")]
            names.extend([base+f"{sl_path}/info.txt",base+f"{sl_path}/remidx_{Path(sl_path).name}.txt"])
            highres=base+f"{np_path}/1k_train.npz"
            asset.update(root=str(root),asset_id=f"{SOURCE}:{c}:{mouse}",path=highres,subject=identity,
                         highres_crc32=f'{entries[highres]["crc"]:08x}', highres_bytes=entries[highres]["file_size"],
                         matching_laser_asset=base+f"{sl_path}/laser_{Path(sl_path).name}.mat" in entries)
            if c=="WT":
                names.extend([highres,base+f"{np_path}/traind.csv"])
            download(metadata,inventory,root,names)
            if c=="WT":
                clock={}
                clock_file=root/np_path/"header_receipt.json"
                if clock_file.exists():
                    clock=json.loads(clock_file.read_text())
                else:
                    for signal in ("EEG","EMG"):
                        key=signal.lower()
                        count,digest=recording_dtype(remote,entries[base+f"{sl_path}/{signal}.mat"])
                        clock.update({f"{key}_samples":count,f"{key}_header_sha256":digest})
                    clock_file.write_text(json.dumps(clock,indent=2)+"\n")
                asset["clock_partner"]=clock
            candidates.append(asset)
    return candidates


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata",type=Path,required=True)
    p.add_argument("--inventory",type=Path,required=True)
    p.add_argument("--root",type=Path,required=True)
    p.add_argument("--scope",choices=("first","WT","all"),required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--config",type=Path,default=Path("config/sleep_dandi.json"))
    p.add_argument("--mode",choices=("pilot","cohort","audit"),required=True)
    args=p.parse_args()
    data=json.loads(args.metadata.read_text())
    inventory=json.loads(args.inventory.read_text())
    cfg=json.loads(args.config.read_text())
    candidates=catalog(data,inventory,args.root,args.scope)
    run_stage(candidates,inspect,cfg,args.mode,args.output,2)


if __name__=="__main__":
    main()
