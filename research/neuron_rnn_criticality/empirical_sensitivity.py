"""Separate unit-count and duration sensitivity in already audited cortical recordings."""
import argparse
from hashlib import sha256
import gzip
import json
from pathlib import Path

import numpy as np

from analyse_snn import entry,source


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    cases=[]
    selected={}
    states={}
    for m,blocks,src in source("unit_phase8"):
        if m.get("dandiset")!="000041" or m["region"]!="mPFC" or m["variant"]!="matched_controls":
            continue
        ids=m["unit_ids"]
        if len(ids)<10:
            continue
        key=(m["subject"],m["session"],m["context"])
        if key not in selected:
            seed=int.from_bytes(sha256("/".join(key).encode()).digest()[:4],"little")
            selected[key]=(tuple(ids),np.sort(np.random.default_rng(seed).choice(len(ids),10,False)))
        original,pick=selected[key]
        assert tuple(ids)==original and len(blocks)==3 and all(len(x)==1500 for x in blocks)
        states.setdefault(key,set()).add(m["state"])
        for selection,cols in (("all",np.arange(len(ids))),("N10",pick)):
            series=[x[:,cols].sum(1) for x in blocks]
            choices={"180s":series}
            choices.update({f"w{i+1}_60s":[v] for i,v in enumerate(series)})
            for window,counts in choices.items():
                case_id=f"emp_{m['session']}_{m['state']}_{selection}_{window}"
                entry(cases,case_id,m["context"],m["session"],"empirical_sensitivity",counts,
                      region="mPFC",state=m["state"],animal=m["subject"],model="DANDI-000041",
                      n_units=len(cols),provenance_sha256=src["sha256"],sample=selection,
                      window=window,unit_ids=[ids[i] for i in cols])
    assert len(states)==4 and all(x=={"Wake","NREM"} for x in states.values())
    doc={"physiology_pass":{"prior_recording_qc":True},"cases":cases}
    a.out.write_bytes(gzip.compress(json.dumps(doc,separators=(",",":")).encode(),mtime=0))
    print(f"{len(cases)} fits, {len(states)} independent recordings",flush=True)


if __name__=="__main__":
    main()
