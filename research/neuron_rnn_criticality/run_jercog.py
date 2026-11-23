"""Compile the pinned author simulator and capture native spikes with explicit seeds."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import shlex
import subprocess
from time import monotonic

COMMIT = "11603e5a32e71664aee8b0d2945ef7ec40d61b37"
SOURCE = "https://github.com/djercog/JercogEtAl2017.git"
SEEDS = ((1729, 3401), (2718, 3402), (3141, 3403))
CONDITIONS = {"figure7": "1", "figure7_supp2": "0"}
DT = .00005
SEDGE = 10.0


def digest(path):
    with path.open("rb") as source:
        from hashlib import file_digest
        return file_digest(source, "sha256").hexdigest()


def replace_once(source, before, after):
    if source.count(before) != 1:
        raise ValueError(f"Pinned source structure changed: {before}")
    return source.replace(before, after)


def verify_source(upstream):
    revision = subprocess.check_output(["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
    if revision != COMMIT:
        raise ValueError(f"Expected upstream {COMMIT}, found {revision}")
    contents = {}
    for name in ("main.cpp", "functions.cpp", "functions.h"):
        contents[name] = subprocess.check_output(
            ["git", "-C", str(upstream), "show", f"{COMMIT}:spikingModel/{name}"]
        ).decode().replace("\r\n", "\n")
    return contents


def compile_simulator(upstream, build, interventions=False):
    source = verify_source(upstream)
    build.mkdir(parents=True, exist_ok=True)
    edits = [
        ("#include <vector>", "#include <vector>\n#include <array>"),
        ("double iniSeed=time(NULL);", "double iniSeed=atoi(argv[3]);"),
        ("double timeTot=100000;", "double timeTot=atof(argv[5])*1000.0;"),
        ("int matDelaySize=1000;", "int delaySeed=atoi(argv[4]);\n\tint matDelaySize=1000;"),
        ("\t\tmatDelayE[i]=(int)", "\t\tr8_uniform_01(&seed);\n\t\tr8_uniform_01(&seed);\n\t\tmatDelayE[i]=(int)"),
        ("(1.0+nDelay) *r8_uniform_01 (&seed)", "(1.0+nDelay) *r8_uniform_01 (&delaySeed)"),
        ("(1.0+(nDelay/2.0)) *r8_uniform_01 (&seed)", "(1.0+(nDelay/2.0)) *r8_uniform_01 (&delaySeed)"),
        ("double commNoiseTime[50+(int)(timeTot/lambdaKicks)][2];",
         "vector<array<double,2>> commNoiseTime(1);"),
        ("\t\tcommNoiseTime[nrKicks][0]= tmKick;\n\t\tcommNoiseTime[nrKicks][1]= commonNoise;",
         "\t\tcommNoiseTime.push_back({double(tmKick), commonNoise});"),
        ("int nrCellsE=800, nrCellsI=200;", "int nrCellsE=Ne, nrCellsI=Ni;"),
        ('\t\tpFile = fopen(fileRasterI, "w+b");', '\t\tfclose(pFile);\n\t\tpFile = fopen(fileRasterI, "w+b");'),
    ]
    if interventions:
        edits += [
            ("\t// integration times", """    // Protocol interventions, supplied only after the baseline gate.
    Jee *= atof(argv[6]); Jie *= atof(argv[6]); Ja *= atof(argv[7]);
    Ie += atof(argv[8]); Ii += atof(argv[8]);
    Jei *= atof(argv[9]); Jii *= atof(argv[9]);
    const double pulseSec=atof(argv[10]);
    const double stopSec=atof(argv[11]);
    // integration times"""),
            ("timestep < nTime; timestep++", "timestep < round(stopSec*1000/dt); timestep++"),
            ("+ Ie + commNoise)/tau_e;",
             "+ Ie + commNoise + ((pulseSec>=0 && timestep>=round(pulseSec*1000/dt) && timestep<round((pulseSec+.010)*1000/dt)) ? 4.0 : 0.0))/tau_e;"),
        ]
    for before, after in edits:
        source["main.cpp"] = replace_once(source["main.cpp"], before, after)
    # The upstream signed int overflows after 21.47 s at 5000 cells and dt=0.05 ms.
    source["functions.cpp"] = replace_once(source["functions.cpp"], "static int used = -1;",
                                          "static long long used = -1;")
    for name, code in source.items():
        (build / name).write_text(code, encoding="utf-8")
    command = ["g++", "-std=gnu++17", "-O3", "-flto", "main.cpp", "functions.cpp", "-lm", "-o", "jercog"]
    subprocess.run(command, cwd=build, check=True, timeout=120)
    return build / "jercog", {
        "command": shlex.join(command),
        "compiler": subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
        "source_sha256": {name: digest(build / name) for name in source},
        "binary_sha256": digest(build / "jercog"),
    }


def simulate(executable, build_info, output, condition, seeds, seconds, setting=None):
    connection, noise = seeds
    path = output / f"{condition}-c{connection}-n{noise}"
    path.mkdir(exist_ok=False, parents=True)
    command = [str(executable), CONDITIONS[condition], "1", str(noise), str(connection), str(seconds)]
    stop=seconds
    params={"author_switch": CONDITIONS[condition]}
    if setting is not None:
        command.extend(map(str,setting))
        stop=setting[-1]
        params.update(dict(zip(("exc_recurrence_multiplier","adaptation_multiplier",
                                "mean_current_delta","inhib_recurrence_multiplier",
                                "pulse_time_seconds","stop_seconds"),setting)))
    started = monotonic()
    with (path / "console.log").open("wb") as log:
        subprocess.run(command, cwd=path, check=True, stdout=log, stderr=subprocess.STDOUT,
                       timeout=max(600, 24 * seconds))
    files = {p.name: digest(p) for p in sorted(path.glob("*.bin"))}
    if len(files) != 7:
        raise ValueError(f"{path.name}: expected seven native outputs, got {list(files)}")
    manifest = dict(
        model_name="JercogEtAl2017", model_commit=COMMIT, condition=condition,
        parameter_config=params,
        simulation_seed=noise, connectivity_seed=connection, duration_seconds=stop,
        burn_in_seconds=SEDGE if stop >= 30 else 0, integration_step_seconds=DT,
        neuron_counts={"E": 4000, "I": 1000}, files=files,
        network_metadata={"topology": "all-to-all with population-specific delayed input",
                          "neurons": "E adapted LIF and I LIF"},
        provenance={"source": SOURCE, "build": build_info,
                    "simulation_argv": command[1:], "wall_seconds": round(monotonic() - started, 2)},
    )
    (path / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(path.name, f"{manifest['provenance']['wall_seconds']:.1f} s wall", flush=True)
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--upstream", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--seconds", type=int, default=100)
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--smoke", action="store_true")
    a = p.parse_args()
    if a.seconds < 5 or a.seconds % 5:
        p.error("duration must be a positive multiple of five seconds (author's output buffer)")
    root = a.out.resolve()
    executable, info = compile_simulator(a.upstream.resolve(), root / "build")
    jobs = [("figure7", (1729, 1729)), ("figure7_supp2", (1729, 1729))] if a.smoke else [
        (condition, seed) for seed in SEEDS for condition in CONDITIONS
    ]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        list(pool.map(lambda job: simulate(executable, info, root, *job, a.seconds), jobs))


if __name__ == "__main__":
    main()
