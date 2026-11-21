"""Replay pair-specific 180 s rows from the verified 000939 Phase 10 count archive."""
import argparse
from collections import defaultdict
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile
import json

from build_qc_fixture import cells, count_entries, save_prepared
from build_atlas_data import ROOT, read, write
from prepare_crcns_sleep import EQUAL_DURATION_STATES, equal_duration_controls
from prepare_sleep_dandi import whole_windows

OUT = "results/validation/phase11_1/"
REFERENCE = "results/phase10/fits/posub-adn-supplement.json"


def archived_counts(archive):
    with ZipFile(archive) as z:
        receipt = json.loads(z.read("cohort/phase10_input.json"))
        raw = z.read("cohort/" + receipt["fixture"])
    assert sha256(raw).hexdigest() == receipt["sha256"]
    return raw, receipt["sha256"]


def spans_for(m, segments):
    intervals = m["intervals_s"]
    if intervals and isinstance(intervals[0], (int, float)):
        intervals = [intervals]
    out = [(round(a / .04), round(b / .04)) for a, b in intervals]
    assert all(b - a == len(s) for (a, b), s in zip(out, segments, strict=True))
    assert all(a >= previous[1] for previous, (a, _) in zip(out, out[1:]))
    return out


def build(archive, out):
    raw, original_sha = archived_counts(archive)
    published = read(REFERENCE)
    rows = {(r["record"]["asset_id"], r["record"]["state"]): r for r in published["records"]
            if r["record"]["dandiset"] == "000939" and r["record"]["variant"] == "reference"}
    cohort = defaultdict(dict)
    for m, segments in count_entries(raw):
        if m["variant"] != "reference":
            continue
        asset, state = m["asset_id"], m["state"]
        source = rows[asset, state]
        assert source["record"] == m and m["region"] == "PoSub" and m["context"] == "REST"
        fit = source["fit"]
        lengths = fit["segment_bins"]
        assert [len(s) for s in segments] == (lengths if isinstance(lengths, list) else [lengths])
        assert sum(int(s.sum()) for s in segments) == round(fit["rate_hz"] * fit["duration_s"])
        spans = spans_for(m, segments)
        cohort[asset][state] = m, dict(zip(spans, segments, strict=True))
    assert len(cohort) == 29 and sum(map(len, cohort.values())) == len(rows)
    selected, eligibility = [], []
    for asset, states in sorted(cohort.items()):
        m = next(iter(states.values()))[0]
        assert all(r["unit_ids"] == m["unit_ids"] for r, _ in states.values())
        windows = {s: whole_windows(list(segments),1500) for s, (_, segments) in states.items()}
        audit = {k:m[k] for k in ("dandiset","region","subject","session","asset_id")}
        audit["candidates"] = {s:len(v) for s,v in windows.items()}
        audit["exclusions"] = {
            v: {s: ("缺少有效参考状态" if s not in windows else f"只有 {len(windows[s])} 个完整 60 秒窗")
                for s in sc if len(windows.get(s, ())) < 3}
            for v, sc in EQUAL_DURATION_STATES.items()
        }
        eligibility.append(audit)
        for variant, state, chosen in equal_duration_controls(windows):
            record, segments = states[state]
            counts = []
            for a,b in chosen:
                for (x,y), s in segments.items():
                    if x <= a < b <= y:
                        counts.append(s[a-x:b-x])
                        break
                else:
                    raise AssertionError((asset, state, a, b))
            assert [len(s) for s in counts] == [1500] * 3
            selected.append(({**record,"variant":variant,"intervals_s":[[a*.04,b*.04] for a,b in chosen]},
                             cells(counts)))
    ids = {m["asset_id"] for m,_ in selected}
    info = save_prepared(selected,out,"states",sorted(ids),
                         source_sha256={"000939_phase10_counts":original_sha,
                                        "phase10_primary_fits":sha256((ROOT/REFERENCE).read_bytes()).hexdigest()})
    out.with_suffix(".sha256").write_text(f"{info['sha256']}  {out.relative_to(ROOT)}\n")
    write(ROOT/(OUT+"pairwise_equal_duration/eligibility.json"),eligibility)
    print(f"{len(selected)} verified count rows, eligible: " +
          ", ".join(f"{v}={sum(not x['exclusions'][v] for x in eligibility)}" for v in EQUAL_DURATION_STATES))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    build(a.archive,a.out.resolve()/"validation_input.mat")
