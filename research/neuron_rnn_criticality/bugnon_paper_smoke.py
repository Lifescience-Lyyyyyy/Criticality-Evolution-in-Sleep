"""Verify the mGluR extension from the paper-linked deNEST example."""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path

from run_bugnon import revision

COMMIT = "d8facf6323aa95151398ee4945623147e33491a8"
EXAMPLE = Path("example_networks/MC_sleep")
MODEL = "ht_neuron_minis_GABA_B1a_soft_gap_noprelimspike_mGluR"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Checkout of the paper-linked deNEST source")
    parser.add_argument("module", type=Path, help="Isolated build of htneuronmodule.so via an ASCII path")
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    revision(args.source, COMMIT)
    module = args.module.absolute()
    if module.name != "htneuronmodule.so" or not module.is_file():
        raise ValueError(f"Missing paper-linked NEST module: {module}")
    os.environ["NEST_MODULE_PATH"] = str(module.parent)
    import denest
    import nest
    import numpy as np
    from denest.network import Network

    assert Path(denest.__file__).resolve().parent.parent == args.source.resolve()
    assert "2.20.0" in nest.version()
    parameters = denest.load_trees(args.source / EXAMPLE / "params/tree_paths.yml")
    tree = parameters.asdict()
    stages = tree["simulation"]["params"]["sessions"]
    sessions = tree["session_models"]
    assert stages == ["warmup", "noise", "grating", "sleep_transition", "sleep"]
    changes = sessions["sleep_transition"]["params"]["unit_changes"]
    transition = {change["population_name"]: change["nest_params"] for change in changes}
    assert transition["Vp_PC"] == {"g_KL": 1.5, "g_peak_AMPA": 1.5, "g_peak_mGluR": 2.0}
    assert transition["Vp_BC"] == {"g_KL": 1.4, "g_peak_AMPA": 1.9, "g_peak_mGluR": 2.0}
    for kind in ("PC", "BC"):
        assert transition["poisson_nrn_to_cortex_mGluR_to_" + kind]["rate"] == 200.0

    nest.SetKernelStatus({"resolution": 0.2, "local_num_threads": 1})
    nest.Install("htneuronmodule")
    network = Network(parameters.children["network"])
    receptor = nest.GetDefaults(MODEL)["receptor_types"]["mGluR"]
    cells = nest.Create(MODEL, 3, params=[{"g_peak_mGluR": p} for p in (-.005, -.010, -.005)])
    spikes = nest.Create("spike_generator", params={"spike_times": [2., 4.]})
    nest.Connect(spikes, cells[:2], syn_spec={"weight": 1., "delay": .2, "receptor_type": receptor})
    meter = nest.Create("multimeter", params={"interval": .2, "record_from": ["V_m", "g_mGluR"]})
    nest.Connect(meter, cells)
    nest.Simulate(12.)
    events = nest.GetStatus(meter, "events")[0]
    standard, doubled, control = [events["g_mGluR"][events["senders"] == gid] for gid in cells]
    assert len(standard) == len(doubled) == len(control) > 50 and np.all(np.isfinite(events["V_m"]))
    max_standard = float(np.max(np.abs(standard)))
    assert max_standard > 0
    np.testing.assert_allclose(doubled, 2 * standard, rtol=3e-6, atol=1e-10)
    np.testing.assert_array_equal(control, np.zeros_like(control))

    receipt = dict(
        status="single_cell_compatibility_pass", source_commit=COMMIT, nest=nest.version(),
        extension_sha256=sha256(module.read_bytes()).hexdigest(), model=MODEL, mglur_receptor=receptor,
        simulation_ms=12., samples_per_cell=len(standard),
        parsed_network=dict(layers=len(network.layers), projections=len(network.projections),
                            population_recorders=len(network.population_recorders)),
        peak_absolute_mglur=[max_standard, float(np.max(np.abs(doubled))), float(np.max(np.abs(control)))],
        example_session_ms={stage: sessions[stage]["params"]["simulation_time"] for stage in stages},
        transition_parameters=transition,
    )
    args.out.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
