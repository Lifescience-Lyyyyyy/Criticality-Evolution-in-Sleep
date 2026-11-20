#!/usr/bin/env python3
"""Recount original 000939 complete windows without changing eligibility rules."""
from statistics import mean, median

from build_atlas_data import ROOT, read, write
from catalog import STATES

FIT_PATH = "results/phase10/fits/posub-adn-supplement.json"


def window_counts(fit, seconds):
    segments = fit["segment_bins"]
    if isinstance(segments, int):
        segments = [segments]
    width = round(seconds / .04)
    return sum(n // width for n in segments)


def audit_windows():
    groups = {}
    unit_counts = {}
    for row in read(FIT_PATH)["records"]:
        m, f = row["record"], row["fit"]
        if m["dandiset"] != "000939" or m["context"] != "REST" or m["region"] != "PoSub":
            continue
        key = m["subject"], m["session"], m["asset_id"]
        groups.setdefault(key, {}).setdefault(m["variant"], {})[m["state"]] = f
        if m["variant"] == "reference":
            assert unit_counts.setdefault(key, len(m["unit_ids"])) == len(m["unit_ids"])
    animals = []
    for (subject, session, asset), methods in sorted(groups.items()):
        ref = methods["reference"]
        counts = {str(width): {s: window_counts(f, width) for s, f in ref.items()} for width in (60, 45)}
        all_ok = min(counts["60"].values()) >= 3
        needed = [s for s in STATES if s in ref and counts["60"][s] < 3]
        aux = methods.get("equal_duration", {})
        assert all_ok == (set(aux) == set(ref)), (subject, counts, list(aux))
        assert {"NREM", "Wake"} <= ref.keys()
        arm = ref["NREM"]["d2"] - ref["Wake"]["d2"]
        animals.append(dict(subject=subject, session=session, asset_id=asset,
                            seconds={s: ref[s]["duration_s"] if s in ref else None for s in STATES},
                            candidate_windows={w: {s: n.get(s) for s in STATES} for w, n in counts.items()},
                            retained=all_ok, required_states=[s for s in STATES if s in ref], blocking_states=needed,
                            nrem_minus_wake=arm, units=unit_counts[subject, session, asset], matched_60_effect=aux["NREM"]["d2"] - aux["Wake"]["d2"] if all_ok else None,
                            feasible_equal_exposure_45=min(counts["45"].values()) >= 4))
    stats = {}
    for name, accepted in (("included", True), ("excluded", False)):
        selected = [r for r in animals if r["retained"] == accepted]
        stats[name] = dict(
            n=len(selected), mean_reference=mean(r["nrem_minus_wake"] for r in selected),
            positive=sum(r["nrem_minus_wake"] > 0 for r in selected),
            median_units=median(r["units"] for r in selected),
            state_median_seconds={s: median(r["seconds"][s] for r in selected if r["seconds"][s] is not None)
                                  for s in STATES},
            median_windows_60={s: median(r["candidate_windows"]["60"][s] for r in selected
                                         if r["candidate_windows"]["60"][s] is not None) for s in STATES},
        )
    assert len(animals) == 29 and stats["included"]["n"] == 15
    stats["matched_60_mean"] = mean(r["matched_60_effect"] for r in animals if r["retained"])
    stats["eligible_at_45_same_180s"] = sum(r["feasible_equal_exposure_45"] for r in animals)
    return dict(source=FIT_PATH, window_gate="at least 3 complete 60 s windows in every available reference state",
                exploratory_45_gate="at least 4 complete 45 s windows in every available reference state",
                animals=animals, summary=stats)


if __name__ == "__main__":
    evidence = audit_windows()
    write(ROOT / "results/validation/000939_equal_duration.json", evidence)
    print(evidence["summary"])
