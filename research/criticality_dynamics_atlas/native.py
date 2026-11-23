"""Plan audited three-state replays, and summarize cached 20 ms unit observations."""
import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from analyse import DT, seed_for, fold, wrap
from inputs import ROOT, key, reference, source, write
from metrics import ar_modes, population_metrics, unit_metrics


def plan(out):
    files=[
      ("results/validation/phase11_1/pairwise_equal_duration/matlab-fits.json","equal_duration_WNR"),
      ("results/sleep_dandi/rem-duration.json","rem_equal_duration")]
    result=defaultdict(list)
    for file,v in files:
        p=ROOT/file
        raw=p.read_bytes()
        archive=json.loads(raw)
        receipt=sha256(raw).hexdigest()
        groups=defaultdict(list)
        for r in archive["records"]:
            m=r["record"]
            if m["variant"]!=v or m["region"] not in {"PoSub","CA1","mPFC","BLA","OFC","ACC","MotorCtx"}:
                continue
            groups[key(m)[:-1]].append(r)
        # Metadata ordering, independent of the fits. Preserve simultaneous CA1/BLA when chosen.
        winners={}
        for k,g in sorted(groups.items()):
            states={r["record"]["state"] for r in g}
            if states!={"Wake","NREM","REM"}:
                continue
            m=g[0]["record"]
            if min(len(np.atleast_2d(r["record"]["intervals_s"])) for r in g)<3:
                continue
            w=(m["dandiset"],m["subject"],m["context"],m["region"])
            winners.setdefault(w,g)
        for g in winners.values():
            for r in g:
                m,f=r["record"],r["fit"]
                fixed= r["long_reach"] or {}
                assert f["exit_status"]==0 and f["solver_flag"]>0 and f["min_ar_root"]>1
                result[m["asset_id"]].append({
                    "meta":wrap(m),"variant":v,"index":sum(map(len,result.values())),
                    "original_spike_total":round(f["duration_s"]*f["rate_hz"]),
                    "d2":f["d2"],"d2_AR10":fixed.get("d2"),
                    "coeff":f["ar_coefficients"],"fit_sha256":receipt,
                    "run":archive["run_id"]})
    write(out,list(result.values()))
    print("Planned",len(result),"assets",sum(map(len,result.values())),"rows")


def raw_pags(cache):
    """Verify checksummed previously extracted 1 ms vectors, then align time bins."""
    batches=defaultdict(list)
    for part in ("brainstem1","brainstem2"):
        refs=reference(part)
        for m,b,receipt in source(part):
            if key(m) not in refs or m["variant"]!="reference" or m["region"]!="PAG":
                continue
            batches[m["subject"]].append((m,b,refs[key(m)]))
    for subject,records in sorted(batches.items()):
        path=cache/"WT"/subject/"1k_train.npz"
        if not path.is_file():
            continue
        check=sha256(path.read_bytes()).hexdigest()
        ids=records[0][0]["unit_ids"]
        assert all(m["unit_ids"]==ids and check==m["digest"] for m,_,_ in records)
        spans=[]
        for m,local,r in records:
            starts=[a+60*i for a,b in np.atleast_2d(m["intervals_s"]) for i in range(int((b-a+1e-7)//60))]
            selected=sorted({0,(len(starts)-1)//2,len(starts)-1}) if starts else []
            spans.append([(starts[i],starts[i]+60) for i in selected])
        separate=[[np.empty((3000,len(ids)),dtype=np.uint16) for _ in ss] for ss in spans]
        with np.load(path,allow_pickle=False) as z:
            for j,uid in enumerate(ids):
                data=z[uid]
                for ss,blocks in zip(spans,separate,strict=True):
                    for (a,b),x in zip(ss,blocks,strict=True):
                        i,k=round(a*1000),round(b*1000)
                        x[:,j]=data[i:k].reshape(-1,20).sum(1)
        for (m,b,r),ss,blocks in zip(records,spans,separate,strict=True):
            if len(ss)==0:
                continue
            # The published population fixture contains the exact same 40 ms bins.
            for (start,stop),value in zip(ss,blocks,strict=True):
                for (a,end),archived in zip(np.atleast_2d(m["intervals_s"]),b,strict=True):
                    if a<=start and stop<=end+1e-8:
                        offset=round((start-a)/DT)
                        np.testing.assert_array_equal(fold([value.sum(1)],2)[0],archived[offset:offset+1500])
                        break
                else:
                    raise ValueError((m["session"],m["state"],ss))
            yield dict(meta=wrap(m)|{"intervals_s":ss},index=None,d2=None,d2_AR10=None,
                       fixture_sha256=check),blocks


def calc(row,blocks,kind):
    m=row["meta"]
    base=fold(blocks,2)
    counts=[x.sum(1) for x in base]
    k=key(m)
    s=seed_for(k[:-1])
    modes=ar_modes(row["coeff"],DT) if row.get("coeff") else []
    return dict(meta=m,index=row["index"],layer=kind,duration=sum(map(len,base))*DT,
                units=len(m["unit_ids"]),d2=row["d2"],d2_AR10=row["d2_AR10"],
                tau_AR_longest=modes[0]["tau"] if modes else None,modes=modes,
                raw_receipt={k:v for k,v in row.items() if k not in {"meta","coeff"}},
                **population_metrics(counts,DT),**unit_metrics(base,DT,s,True),
                bins=[dict(bin_s=dt,**population_metrics([u.sum(1) for u in fold(blocks,step)],dt))
                      for step,dt in ((1,.02),(4,.08),(10,.20))])


def verified_blocks(archive):
    """Read the already-audited NWB count blocks with their original receipt."""
    from io import BytesIO
    with ZipFile(archive) as z:
        rec=json.loads(z.read("receipt.json"))
        for r in rec["verified"]:
            raw=z.read(f"{r['index']}.npz")
            assert sha256(raw).hexdigest()==r["matrix_sha256"]
            with np.load(BytesIO(raw),allow_pickle=False) as b:
                yield r,[b[n].astype(float) for n in b.files]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=["plan","analysis"])
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--archive",type=Path,action="append",default=[])
    p.add_argument("--brainstem",type=Path)
    a=p.parse_args()
    if a.mode=="plan":
        plan(a.out)
        return
    rows=[]
    for archive in a.archive:
        for r,blocks in verified_blocks(archive):
            rows.append(calc(r,blocks,"NWB_matched"))
        with ZipFile(archive) as z:
            write(a.out.with_name(f"{archive.stem}-receipt.json"),json.loads(z.read("receipt.json")))
    if a.brainstem:
        for row,blocks in raw_pags(a.brainstem):
            rows.append(calc(row,[x.astype(float) for x in blocks],"PAG_1ms"))
    write(a.out,rows)
    print("Native metrics",len(rows))


if __name__=="__main__":
    main()
