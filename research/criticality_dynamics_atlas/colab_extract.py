"""Fetch only sorted spike columns for the exact archived 60 s windows in plan.json."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
import json
from pathlib import Path
from time import sleep
from urllib.request import urlopen
from zipfile import ZipFile, ZIP_DEFLATED

import h5py
import numpy as np
import remfile

VERSIONS={"000041":"0.261104.0300","000552":"0.261113.1035",
          "000061":"0.250624.0426","000939":"0.261016.2310"}


def fetch(rows):
    first=rows[0]["meta"]
    ds=first["dandiset"]
    url=f"https://api.dandiarchive.org/api/dandisets/{ds}/versions/{VERSIONS[ds]}/assets/{first['asset_id']}/"
    for n in range(3):
        try:
            with urlopen(url,timeout=35) as f:
                meta=json.load(f)
            known=first["digest"]
            assert set(known.values() if isinstance(known,dict) else [known]) <= set(meta["digest"].values()),url
            storage,=[u for u in meta["contentUrl"] if u.startswith("https://dandiarchive.s3.amazonaws.com/")]
            with remfile.File(storage) as r, h5py.File(r,"r") as f:
                u=f["units"]
                ids={str(v):i for i,v in enumerate(u["id"][:])}
                edges=np.r_[0,u["spike_times_index"][:]]
                raw=u["spike_times"][:]
                trials={str(i):raw[edges[ids[str(i)]]:edges[ids[str(i)]+1]]
                        for i in {i for x in rows for i in x["meta"]["unit_ids"]}}
            output=[]
            for row in rows:
                m=row["meta"]
                result=[]
                for a,b in np.atleast_2d(m["intervals_s"]):
                    start,stop=round(a/.02),round(b/.02)
                    matrix=[]
                    for i in m["unit_ids"]:
                        t=trials[str(i)]
                        ia,ib=np.searchsorted(t,(start*.02-1e-10,stop*.02-1e-10))
                        bins=np.floor(t[ia:ib]/.02+1e-9).astype(int)-start
                        assert not len(bins) or (bins[0]>=0 and bins[-1]<stop-start)
                        matrix.append(np.bincount(bins,minlength=stop-start))
                    v=np.column_stack(matrix)
                    assert v.max()<65536
                    result.append(v.astype(np.uint16))
                assert sum(v.sum() for v in result)==row["original_spike_total"],(m["session"],m["state"])
                p=Path(f"{row['index']}.npz")
                np.savez_compressed(p,*result)
                output.append(row | {"matrix_sha256":sha256(p.read_bytes()).hexdigest(),
                                     "source_version":VERSIONS[ds],"source_receipt":url})
            print("VERIFIED",ds,first["session"],len(output),flush=True)
            return output
        except (OSError,TimeoutError):
            if n==2:
                raise
            sleep(3)


def main():
    cases=json.loads(Path("plan.json").read_text())
    yes,no=[],[]
    with ThreadPoolExecutor(max_workers=6) as pool:
        fs={pool.submit(fetch,x):x[0]["meta"]["session"] for x in cases}
        for run in as_completed(fs):
            try:
                yes.extend(run.result())
            except Exception as e:
                no.append(dict(session=fs[run],type=type(e).__name__,details=str(e)[:350]))
                print("UNAVAILABLE",no[-1],flush=True)
    Path("receipt.json").write_text(json.dumps(dict(verified=sorted(yes,key=lambda x:x["index"]),
        unavailable=no,read_scope="NWB units/id and ragged spike_times only; no neural signal download")))
    with ZipFile("unit-windows.zip","w",ZIP_DEFLATED) as z:
        z.write("receipt.json")
        for row in yes:
            z.write(f"{row['index']}.npz")
    print("RESULT",len(yes),len(no),Path("unit-windows.zip").stat().st_size)


if __name__=="__main__":
    main()
