"""Small read-only export of Phase 11.1 assessment and verified sensitivities."""
from catalog import BASE

P111 = "results/validation/phase11_1/"
CONTRASTS = ("NW", "RN", "RW")


def export_preliminary(read):
    evidence=read(P111+"evidence_assessment.json")
    pairs=read(P111+"pairwise_equal_duration/summary.json")
    sorting=[]
    for row in read("results/validation/audit_comparisons.json")["cases"]:
        if row["case"] not in ("B","G"):
            continue
        area="PoSub" if row["case"]=="B" else "LDT"
        ds=row["dataset"]
        models={f["variant"]:f for f in row["fits"] if f["region"]==area and f["context"]=="REST"}
        alt=[v for v in models if v.startswith(("lratio_", "contamination_"))]
        assert alt
        controls=[dict(id=v,units=next(iter(models[v]["states"].values()))["units"],effects=models[v]["effects"])
                  for v in sorted(alt,reverse=True)]
        ref=models["reference"]
        provenance=next(iter(models[alt[0]]["states"].values()))["provenance"]
        sorting.append(dict(cohort=f"mouse|{area}|{ds}|REST",session=row["session"],variants=controls,
                            reference=dict(units=next(iter(ref["states"].values()))["units"],effects=ref["effects"]),
                            run=f"{BASE}/actions/runs/{provenance['run_id']}",
                            source=f"{BASE}/blob/main/results/validation/audit_comparisons.json"))
    extra=read(P111+"sorting_checks/summary.json")
    for case in extra["cases"]:
        if not case.get("alternate"):
            continue
        sorting.append(dict(cohort="mouse|"+case["region"]+"|"+case["dandiset"]+"|"+case["context"],
                            session=case["session"],
                            variants=[dict(id=case["method"],units=case["alternate_units"],effects=case["alternate"])],
                            reference=dict(units=case["original_units"],effects=case["reference"]),
                            run=extra["provenance"]["action"],source=f"{BASE}/blob/main/{P111}sorting_checks/summary.json"))
    display={}
    for c in CONTRASTS:
        row=pairs["contrasts"][c]
        display[c]={k:row[k] for k in ("legacy","corrected","three_state_fixed","matched_animals",
                                      "corrected_minus_legacy_same_animals")}
        for mode in ("legacy","corrected","three_state_fixed"):
            display[c][mode]={k:v for k,v in row[mode].items() if k!="animal_effects"}
    return dict(evidence={r["id"]:dict(category=r["category"],reason=r["reason"]) for r in evidence},
                pairwise=dict(cohort="mouse|PoSub|000939|REST",contrasts=display,provenance=pairs["provenance"],
                              source=f"{BASE}/blob/main/{P111}pairwise_equal_duration/summary.json"),
                sorting=sorting)
