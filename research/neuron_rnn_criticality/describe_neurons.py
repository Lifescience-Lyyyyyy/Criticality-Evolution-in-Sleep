"""Summarize rate dependence and event statistics per independent network."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from physiology import digest, write
from metrics import events as count_events


def describe_cells(rate, tau, tau_integrated):
    mask = np.isfinite(tau)
    valid_rates, valid_tau = rate[mask], tau[mask]
    rho = (spearmanr(valid_rates, valid_tau).statistic
           if len(valid_tau) > 2 and np.ptp(valid_rates) and np.ptp(valid_tau) else None)
    return dict(
        n=len(rate), n_fit=int(mask.sum()),
        rate_hz_p10_p50_p90=np.percentile(rate, [10, 50, 90]).tolist(),
        fit_rate_hz_p10_p50_p90=np.percentile(valid_rates, [10, 50, 90]).tolist() if mask.any() else None,
        tau_s_p10_p50_p90=np.percentile(valid_tau, [10, 50, 90]).tolist() if mask.any() else None,
        integrated_tau_s_median=float(np.nanmedian(tau_integrated)) if np.isfinite(tau_integrated).any() else None,
        spearman_rate_tau=rho,
    )


def describe(run, counts):
    cells = np.asarray(run["individual"], dtype=float)
    metrics = {name:describe_cells(*rows.T)
               for name, rows in (("E", cells[:4000]), ("I", cells[4000:]))}
    events = []
    dt = counts[run["id"]]["bin_seconds"]
    original = np.asarray(counts[run["id"]]["counts"][0])
    threshold = float(np.quantile(original,.90))
    for record in [run, *(r for r in run["controls"] if r["kind"] == "shift")]:
        q = record["population"]
        blocks=[np.asarray(a) for a in counts[record["id"]]["counts"]]
        bins = dict(silent_fraction=q["silent_fraction"], silence=q["off"],
                    burst_own_p90=q["burst"], burst_original_p90=count_events(blocks,dt,"high",threshold))
        events.append(dict(id=record["id"], **bins))
    return dict(id=run["id"], condition=run["manifest"]["condition"],
                connectivity_seed=run["manifest"]["connectivity_seed"], cell_types=metrics, burst_threshold_count_per_bin=threshold, events=events)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("analysis", type=Path)
    p.add_argument("counts", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    receipt=json.loads((a.counts.parent/"baseline_matlab_receipt.json").read_text())
    assert digest(a.counts)==receipt["paired_input_sha256"]
    with gzip.open(a.counts,"rt") as f:
        cases={r["id"]:r for r in json.load(f)["cases"]}
    rows = [describe(r,cases) for r in json.loads(a.analysis.read_text())["runs"]]
    write(a.out, rows)
    for r in rows:
        details = [(k,q["n_fit"],q["spearman_rate_tau"]) for k,q in r["cell_types"].items()]
        events = r["events"]
        print(r["id"], details)
        print([(round(x["silent_fraction"],4),x["silence"]["p95_s"],x["burst_own_p90"]["p95_s"],x["burst_original_p90"]["p95_s"]) for x in events])


if __name__ == "__main__":
    main()
