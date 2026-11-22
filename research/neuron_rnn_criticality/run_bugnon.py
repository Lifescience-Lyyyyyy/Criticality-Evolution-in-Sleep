"""Run the two published Bugnon reference parameter trees using pinned deNEST."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from time import monotonic

REPLICATE="028ce0aad1b9a9fa002341a43728c2be79bbfc48"
DENEST="9bbec4c9b45c028c651429a08c140921b2d44e16"
NEST="37722d5e35a19bcd2103025b7a581cc76f6dd5fc"
CASES={"wake":"fig0_wake_mode_spontaneous","sleep":"fig0_sleep_mode"}


def revision(path,expected):
    current=subprocess.check_output(["git","-C",str(path),"rev-parse","HEAD"],text=True).strip()
    if current!=expected:
        raise ValueError(f"{path} is at {current}, expected {expected}")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("condition",choices=tuple(CASES))
    p.add_argument("--upstream",required=True,type=Path)
    p.add_argument("--denest",required=True,type=Path)
    p.add_argument("--nest-source",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    for path,expected in ((a.upstream,REPLICATE),(a.denest,DENEST),(a.nest_source,NEST)):
        revision(path,expected)
    import nest
    import denest
    assert "2.20.0" in nest.version()
    source=a.upstream.resolve()/"simulations"/CASES[a.condition]/"parameter_tree.yml"
    raw=source.read_bytes()
    params=denest.ParamsTree.read(source)
    out=a.out.resolve()
    if out.exists() and list(out.iterdir()):
        raise FileExistsError(out)
    # The two author reference cases use no grating-array inputs.
    sim=denest.Simulation(params,output_dir=out,
                          input_dir=a.upstream.resolve()/"network_default_params"/"input")
    t=monotonic()
    sim.run()
    out.mkdir(parents=True,exist_ok=True)
    receipt={"model_name":"BugnonEtAl2024","condition":a.condition,
        "parameter_source":str(source.relative_to(a.upstream.resolve())),
        "parameter_sha256":sha256(raw).hexdigest(),"model_commit":REPLICATE,
        "denest_commit":DENEST,"nest_commit":NEST,"nest_version":nest.version(),
        "resolution_ms":nest.GetKernelStatus("resolution"),
        "threads":nest.GetKernelStatus("local_num_threads"),
        "global_seed":nest.GetKernelStatus("grng_seed"),
        "vp_seeds":list(nest.GetKernelStatus("rng_seeds")),
        "wall_seconds":round(monotonic()-t,2)}
    (out/"run_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt),flush=True)


if __name__=="__main__":
    main()
