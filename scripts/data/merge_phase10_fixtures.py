"""Combine independently audited MATLAB inputs without changing counts or intervals."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from scipy.io import loadmat

from prepare_sleep_dandi import save_prepared


def load(folder):
    receipt=json.loads((folder/"phase10_input.json").read_text())
    parts=receipt.get("parts") or [{"name":receipt["fixture"],"sha256":receipt["sha256"]}]
    rows=[]
    for part in parts:
        path=folder/part["name"]
        if sha256(path.read_bytes()).hexdigest()!=part["sha256"]:
            raise ValueError(f"Input digest failure: {path}")
        mat=loadmat(path,struct_as_record=False,squeeze_me=False)["cohort"][0,0]
        for field in mat._fieldnames:
            item=getattr(mat,field)[0,0]
            rows.append((json.loads(item.metadata[0]),item.counts))
    if len(rows)!=receipt["records"]:
        raise ValueError(f"Source receipt does not match MAT: {folder}")
    return receipt,rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs",type=Path,nargs="+")
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    inputs=[load(path) for path in a.inputs]
    modes={receipt["mode"] for receipt,_ in inputs}
    if len(modes)!=1:
        raise ValueError("Fixture inputs must have the same mode")
    mode=modes.pop()
    rows=[row for _,source in inputs for row in source]
    if mode=="controls":
        strongest={}
        for row,_ in rows:
            key=row["dandiset"],row["region"],row["subject"]
            score=(len(row["unit_ids"]),row["asset_id"])
            strongest[key]=max(score,strongest.get(key,score))
        rows=[(r,c) for r,c in rows if r["asset_id"]==strongest[r["dandiset"],r["region"],r["subject"]][1]]
    keys=[tuple(r[k] for k in ("dandiset","asset_id","context","region","state","variant")) for r,_ in rows]
    if len(keys)!=len(set(keys)):
        raise ValueError("Repeated input record")
    received=[{"receipt_sha256":sha256((path/"phase10_input.json").read_bytes()).hexdigest(),
               "receipt":receipt} for path,(receipt,_) in zip(a.inputs,inputs,strict=True)]
    extra={"constituent_inputs_json":json.dumps(received,sort_keys=True)}
    if mode=="controls":
        extra.update(detailed_controls=True,repetitions=40)
    print(save_prepared(rows,a.out,mode,sorted({r["asset_id"] for r,_ in rows}),**extra))


if __name__=="__main__":
    main()
