"""Read existing SHA-pinned count fixtures without copying them into Git."""
from collections import defaultdict
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[2]
DIR=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts/data"))
from prepare_sleep_dandi import save_prepared

A="data/fixtures/"
FILES={
 "unit_phase8":("ec3b391a0dbdcb554fb61de6c15e01a7d75269f1",A+"sleep_dandi_input.sha256"),
 "unit_phase10":("56952e69678cbc732b438c463cee25eaf40f570d",A+"phase10/phase10_input.sha256"),
 "000041":("5fe4ab5efbb4888dd63f7b2972d4195a628c83b0",A+"sleep_dandi_input.sha256"),
 "000061":("a67a750a00737117528dad9d68a8dce3267f5530",A+"sleep_dandi_input.sha256"),
 "000552":("8b4212e2f0fc4abee4f6995d16810a2ff1a2d321",A+"sleep_dandi_input.sha256"),
 "000056":("2e0fb612e478f51f4dde4136189d3dc5619527fa",A+"phase10/phase10_input.sha256"),
 "PoSub_ADN":("51d57c648998940305a0de7b77e691496731c7f9",A+"phase10/phase10_input.sha256"),
 "brainstem1":("89664ff4042d14e1d9a665f8a6592fc038d3583d",A+"phase10/phase10_input.sha256"),
 "brainstem2":("e6ae689e77d37afd46d7cbb97fa1fed6c601d334",A+"phase10/phase10_input.sha256"),
 "quiet":("9a1e1c0613ab6acbe64ccab31ddc448ec7c8fb06",A+"phase10/phase10_input.sha256"),
 "hc23":("ff32d190d1c3ec385288088f30886ab81b53d609",A+"phase9/crcns_input.sha256"),
 "hc24rest":("75ae1b2104c637260059f72e47ca269233b8cbe9",A+"phase9/crcns_input.sha256"),
 "hc24run":("9177741a4055da019a28429ea61543ff56936310",A+"phase9/crcns_input.sha256")
}
FITS={
 "000041":"results/sleep_dandi/000041/cortex.json",
 "000061":"results/sleep_dandi/000061/cross.json",
 "000552":"results/sleep_dandi/000552/mouse-ca1.json",
 "000056":"results/phase10/fits/a_cohort.json",
 "PoSub_ADN":"results/phase10/fits/posub-adn-supplement.json",
 "brainstem1":"results/phase10/fits/c_cohort.json",
 "brainstem2":"results/phase10/fits/c-rescue_cohort.json",
 "quiet":"results/phase10/fits/quiet-wake.json",
 "hc23":"results/sleep_crcns/hc-23/cohort.json",
 "hc24rest":"results/sleep_crcns/hc-24/pilot.json",
 "hc24run":"results/sleep_crcns/hc-24/post-run.json",
}


def key(m):
    return tuple(m[k] for k in ("dandiset","subject","session","region","context","state"))


def git(*args):
    return subprocess.check_output(["git","-C",str(ROOT),*args])


def source(label,fetch=False):
    commit, checksum=FILES[label]
    if fetch:
        try:
            git("cat-file","-e",commit+":"+checksum)
        except subprocess.CalledProcessError:
            git("fetch","--depth=1","--no-tags","origin",commit)
    for l in git("show",commit+":"+checksum).decode().splitlines():
        expected, name=l.split()
        raw=git("show",commit+":"+name)
        assert sha256(raw).hexdigest()==expected,(label,name)
        for x in loadmat(BytesIO(raw),simplify_cells=True)["cohort"].values():
            m=json.loads(x["metadata"])
            if m["variant"] not in {"reference","matched_controls"}:
                continue
            a=np.atleast_1d(x["counts"])
            blocks=[v.astype(np.float64) for v in a] if a.dtype==object else [a.astype(np.float64)]
            if m["variant"]=="matched_controls":
                blocks=[v.T for v in blocks]
                assert all(v.shape==(1500,len(m["unit_ids"])) for v in blocks)
            yield m,blocks,dict(commit=commit,fixture=name,sha256=expected)


def reference(label):
    p=ROOT/FITS[label]
    d=json.loads(p.read_bytes())
    rows=d["records"]
    return {key(r["record"]):r for r in (rows if isinstance(rows,list) else [rows])
            if r["record"]["variant"]=="reference"}


def shrink(x):
    if isinstance(x,dict):
        return {k:shrink(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):
        return [shrink(v) for v in x]
    if isinstance(x,np.integer):
        return int(x)
    if isinstance(x,(float,np.floating)):
        return round(float(x),8) if np.isfinite(x) else None
    return x


def write(p,doc):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(shrink(doc),ensure_ascii=False,allow_nan=False,separators=(",",":"))+"\n")


def export(out,groups):
    plan=[]
    for name,dt,p,alt,rows in groups:
        if not rows:
            continue
        output=out/name/"counts.mat"
        series=[]
        for m,blocks in rows:
            entries=np.empty((len(blocks),1),dtype=object)
            for i,x in enumerate(blocks):
                entries[i,0]=x
            series.append((m,entries))
        doc=save_prepared(series,output,"states",sorted({m["session"] for m,_ in rows}))
        assert doc["fixture"]==output.name
        plan.append(dict(name=name,spec=[dt,p,alt],receipt=str(output.with_suffix(".json")),
                         input_sha256=doc["sha256"],records=len(rows)))
    write(out/"plan.json",plan)
