"""Fit the same native windows at four widths, exporting only Yule–Walker moments."""
import argparse
from collections import defaultdict
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from analyse import fold, key_id
from inputs import write
from native import raw_pags, verified_blocks
from summarize import PAIRS, PIN, avg, ci, fit_ok

# Each tuple gives width in ms, number of native 20 ms bins, and alternate AR order.
SPECS = ((20, 1, 10), (40, 2, 10), (80, 4, 2), (200, 10, 1))


def moments(blocks, order):
    """Mirror pinned calc_db global sample SD and yw_ar valid-lag exposure."""
    joined = np.concatenate(blocks, dtype=float)
    variance = joined.var(ddof=1)
    if variance <= 0:
        raise ValueError("Constant population counts")
    centered = [x - joined.mean() for x in blocks]
    return [float(sum(np.dot(x[:-k], x[k:]) for x in centered) /
                  (variance * sum(len(x) - k for x in centered)))
            for k in range(1, order + 1)]


def prepare(cache, brainstem, out):
    rows = []
    studies = [("NWB_matched", verified_blocks(cache / f"{batch}-unit-windows.zip"))
               for batch in ("phase8", "posub")]
    studies.append(("PAG_1ms", raw_pags(brainstem)))
    for layer, cases in studies:
        for r, arrays in cases:
            m = r["meta"]
            native = [x.sum(1) for x in arrays]
            # These exposures are equal across states in the verified inputs.
            assert all(len(x) % 20 == 0 for x in native), key_id(m)
            for width, step, alternate in SPECS:
                counts = fold(native, step)
                rows.append(dict(
                    id=f"{key_id(m)}_{width}ms", key=key_id(m), layer=layer,
                    meta={k: m[k] for k in ("dandiset", "subject", "session", "region",
                                            "context", "state", "intervals_s", "asset_id")},
                    width_ms=width, order=5, alternate_order=alternate,
                    units=len(m["unit_ids"]), bins=sum(map(len, counts)),
                    gamma=moments(counts, max(5, alternate)),
                    reference_d2=r["d2"] if width == 40 else None,
                    reference_gamma_coeff=r.get("coeff") if width == 40 else None,
                    source_sha256=r.get("matrix_sha256") or r["fixture_sha256"],
                ))
    assert rows and len({r["id"] for r in rows}) == len(rows)
    assert any(r["layer"] == "PAG_1ms" for r in rows)
    assert len({(r["key"], r["width_ms"]) for r in rows}) == len(rows)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Preserve full floating point precision to compare against archived MATLAB fits.
    out.write_text(json.dumps(rows, allow_nan=False, separators=(",", ":")) + "\n")
    print(f"{len(rows)} fits from {len(rows) // len(SPECS)} native observations")


def summarize(input_path, artifact, receipt_path, out):
    receipt = json.loads(receipt_path.read_bytes())
    raw = artifact.read_bytes()
    assert sha256(raw).hexdigest() == receipt["digest"].removeprefix("sha256:")
    with ZipFile(BytesIO(raw)) as z:
        result = json.loads(z.read("native-d2-fits.json"))
    reference_rows = json.loads(input_path.read_bytes())
    assert result["input_sha256"] == sha256(input_path.read_bytes()).hexdigest()
    assert result["toolkit_commit"] == PIN and "R2024b" in result["matlab"]
    assert (str(result["run_id"]), result["commit"]) == (str(receipt["id"]), receipt["commit"])
    fits = {r["id"]: r for r in result["rows"]}
    assert len(fits) == len(result["rows"]) == len(reference_rows)
    rows, errors, kernels = [], [], []
    groups = defaultdict(dict)
    for r in reference_rows:
        f = fits.pop(r["id"])
        primary, alternate = f["fit"], f["alternate"]
        value = primary["d2"] if fit_ok(primary, r["order"]) else None
        alt = alternate["d2"] if fit_ok(alternate, r["alternate_order"]) else None
        if r["reference_d2"] is not None:
            assert value is not None, r["id"]
            errors.append(abs(value - r["reference_d2"]))
            kernels.append(float(np.max(np.abs(np.array(primary["ar_coefficients"]) -
                                               r["reference_gamma_coeff"]))))
        rows.append({k: v for k, v in r.items() if k not in ("gamma", "reference_gamma_coeff")} |
                    {"d2": value, "d2_alt": alt, "fit": primary, "alternate": alternate})
        m = r["meta"]
        g = (r["layer"], m["dandiset"], m["region"], m["context"], r["width_ms"])
        groups[g, m["subject"], m["session"]][m["state"]] = rows[-1]
    assert errors and max(errors) < 1e-5 and max(kernels) < 1e-9, (max(errors), max(kernels))
    animal = defaultdict(list)
    for (g, subject, session), states in sorted(groups.items()):
        assert len({(s["bins"], s["units"]) for s in states.values()}) == 1, (g, subject, session)
        for a, b in PAIRS:
            if a in states and b in states:
                for metric in ("d2", "d2_alt"):
                    x, y = states[a][metric], states[b][metric]
                    if x is not None and y is not None:
                        animal[g, a + "_" + b, metric, subject].append(x - y)
    comparisons = defaultdict(dict)
    for (g, contrast, metric, subject), values in animal.items():
        comparisons[g, contrast, metric][subject] = avg(values)
    paired = [dict(group=g, contrast=contrast, metric=metric, n=len(values), mean=avg(list(values.values())),
                   ci95=ci(list(values.values())), positive=sum(v > 0 for v in values.values()), animals=values)
              for (g, contrast, metric), values in sorted(comparisons.items())]
    receipt.update(matlab=result["matlab"], input_sha256=result["input_sha256"], n_fits=len(rows),
                   reference_replay_n=len(errors), max_reference_d2_error=max(errors),
                   max_reference_kernel_error=max(kernels))
    write(out, dict(rows=rows, paired=paired, receipt=receipt))
    print(f"Replayed {len(errors)} original d2 fits. Max |error|={max(errors):.3g}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="mode", required=True)
    prep = sub.add_parser("prepare")
    for path in ("cache", "brainstem", "out"):
        prep.add_argument(path, type=Path)
    join = sub.add_parser("summarize")
    for path in ("input", "artifact", "receipt", "out"):
        join.add_argument(path, type=Path)
    a = p.parse_args()
    if a.mode == "prepare":
        prepare(a.cache, a.brainstem, a.out)
    else:
        summarize(a.input, a.artifact, a.receipt, a.out)
