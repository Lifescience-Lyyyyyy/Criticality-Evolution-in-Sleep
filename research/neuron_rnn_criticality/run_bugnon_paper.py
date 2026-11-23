"""Run the paper-linked Bugnon network under spontaneous retinal input."""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path

from bugnon_paper_smoke import COMMIT, EXAMPLE
from run_bugnon import revision


def session_overrides(condition):
    phases = ["warmup", "noise"]
    if condition != "wake":
        phases += ["sleep_transition", "sleep"]
    # The publication used spontaneous retinal input and 120 s stable recordings.
    # The network check only verifies model construction and the parameter switch.
    durations = ({phase: 1.0 for phase in phases} if condition == "network-check"
                 else {phases[-1]: 120000.0})
    return {
        "simulation": {"params": {"sessions": phases}},
        "session_models": {phase: {"params": {"simulation_time": ms}}
                           for phase, ms in durations.items()},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("condition", choices=("wake", "sleep", "network-check"))
    parser.add_argument("source", type=Path)
    parser.add_argument("module", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    revision(args.source, COMMIT)
    module = args.module.absolute()
    if module.name != "htneuronmodule.so" or not module.is_file():
        raise ValueError(f"Missing paper-linked NEST module: {module}")
    out = args.out.absolute()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(out)
    os.environ["NEST_MODULE_PATH"] = str(module.parent)
    import denest
    import nest
    import numpy as np

    assert Path(denest.__file__).resolve().parent.parent == args.source.resolve()
    assert "2.20.0" in nest.version()
    overrides = session_overrides(args.condition)
    phases = overrides["simulation"]["params"]["sessions"]
    example = args.source / EXAMPLE
    parameters = denest.load_trees(example / "params/tree_paths.yml", overrides)
    sim = denest.Simulation(parameters, input_dir=example / "input", output_dir=out)
    sizes = {key: int(nest.GetKernelStatus(key)) for key in ("network_size", "num_connections")}
    assert int(nest.GetKernelStatus("local_num_threads")) == 5
    transition = sim.session_models["sleep_transition"].params["unit_changes"]
    checks = {}
    if "sleep_transition" in phases:
        for change in transition:
            population = change["population_name"]
            layer = sim.network.layers[change["layers"][0]]
            gid = int(layer.gids(population=population)[0])
            keys = change["nest_params"]
            baseline = nest.GetStatus([gid])[0]
            before = {key: baseline[key] for key in keys}
            checks[population] = dict(gid=gid, before=before, mode=change["change_type"], values=keys)

    sim.run()
    for check in checks.values():
        end = nest.GetStatus([check["gid"]])[0]
        after = {key: end[key] for key in check["values"]}
        for key, value in check["values"].items():
            expected = check["before"][key] * value if check["mode"] == "multiplicative" else value
            np.testing.assert_allclose(after[key], expected, atol=1e-12, rtol=1e-10)
        check["after"] = after

    receipt = dict(
        condition=args.condition, source_commit=COMMIT, nest=nest.version(),
        extension_sha256=sha256(module.read_bytes()).hexdigest(),
        initialization=sizes, sessions=sim.session_times,
        runtime_parameter_sha256=sha256((out / "parameter_tree.yml").read_bytes()).hexdigest(),
        override=overrides, checked_transition=checks,
        global_seed=int(nest.GetKernelStatus("grng_seed")),
        vp_seeds=list(nest.GetKernelStatus("rng_seeds")),
    )
    text = json.dumps(receipt, allow_nan=False, indent=2) + "\n"
    temporary = out / "run_receipt.json.tmp"
    temporary.write_text(text)
    temporary.replace(out / "run_receipt.json")
    print(text, flush=True)


if __name__ == "__main__":
    main()
