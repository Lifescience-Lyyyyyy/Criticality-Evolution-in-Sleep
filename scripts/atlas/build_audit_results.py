"""Evaluate the fixed Phase 11 sample against source tables and actual fits."""
from collections import Counter, defaultdict

from build_atlas_data import ROOT, read, write
from catalog import BASE, FIT_FILES, PAIRS, V1_FILE

OUT = "results/validation/"
TABLES = "results/phase10/audit/"
ATLAS = BASE + "/blob/main/"
STAGING = {"awake": "Wake", "wake": "Wake", "wakestate": "Wake", "nrem": "NREM",
           "sws": "NREM", "non-rem": "NREM", "nremstate": "NREM", "rem": "REM", "remstate": "REM"}


def sorted_intervals(raw):
    return [[a, b] for a, b in sorted(raw)]


def spans(value):
    return [value] if value and isinstance(value[0], (float, int)) else value


def intersection(a, b, x, y):
    return max(0., min(b, y) - max(a, x))


def covered(a, b, ranges):
    end = a
    for x, y in ranges:
        if x > end + 1e-8:
            break
        if y > end:
            end = y
        if end >= b - 1e-8:
            return True
    return False


def entries(path):
    doc = read(path)
    rows = doc["records"]
    for r in rows if isinstance(rows, list) else [rows]:
        yield r, {"run_id": str(doc["run_id"]), "commit": doc["project_commit"], "file": path}


def state_summary(states):
    out = {}
    for s, row in states.items():
        r, f, p = row
        assert f["exit_status"] == 0 and f["solver_flag"] > 0 and f["min_ar_root"] > 1
        out[s] = dict(d2=f["d2"], units=len(r["unit_ids"]) if "unit_ids" in r else None,
                      unit_ids=r.get("unit_ids"), cutoff=r.get("lratio_cutoff", r.get("contamination_cutoff_pct")),
                      seconds=f["duration_s"], segment_bins=f["segment_bins"], bins=r.get("intervals_s"),
                      min_root=f["min_ar_root"], solver=f["solver_flag"], kurtosis=f["residual_excess_kurtosis"],
                      residual_acf=f["residual_acf"], observed=f["observed_acf"], modeled=f["model_acf"],
                      acf_rmse=f["acf_rmse"], max_silence_s=f["max_silence_s"], provenance=p)
    effects = {key: out[b]["d2"] - out[a]["d2"] for key, (a, b) in PAIRS.items() if a in out and b in out}
    return dict(states=out, effects=effects)


def reject_reason(u):
    if not u["anatomy_agrees"]:
        return "单位与电极解剖不一致"
    if u.get("outside_spikes"):
        return "尖峰越出采集时钟"
    if u.get("duplicates") or u.get("reversals"):
        return "时间戳重复或逆序"
    if u.get("on_source_grid") is False:
        return "偏离采样格"
    if u["spikes"] < 2:
        return "尖峰过少"
    if u["cell_type"].lower() in {"mua", "mu", "noise", "bad", "unsorted"}:
        return "作者标注为多单位或噪声"
    return "其他，见原逐单元审计"


def scores(case, selected, raw, previous):
    table = {x["path"]: x for x in raw["clocks"]}
    sample = next(x for x in table.values() if x["path"].startswith("/processing/ecephys/LFP/"))
    clock = [sample["start"], sample["end"]]
    if case["dandiset"] == "000939":
        assert table["/acquisition/ElectricalSeries"]["end"] == clock[1]
    if previous and previous.get("clock"):
        assert previous["clock"]["seconds"] == clock
    if case["case"] == "D":
        stages = raw["stages"]["intervals/sleep_states"]
        failed = [r for r in stages if r[1] > clock[1] or r[0] < clock[0]]
        assert failed and raw["spikes_last_s"] > clock[1] + 5000
        return dict(clock=clock, stages=len(stages), out_of_range=failed, latest_spike=raw["spikes_last_s"],
                    retained=0, physical_clock_conflict=True)
    result = {}
    for (context, area, variant), states in selected.items():
        if variant != "reference":
            continue
        for state, (m, f, _) in states.items():
            labels = raw["stages"].get(m["state_path"].lstrip("/"))
            if labels is None:
                continue
            timeline = sorted_intervals([(a, b) for a, b, name in labels if b > a and STAGING.get(name.lower()) == state])
            annotations = spans(m["intervals_s"])
            for a, b in annotations:
                assert clock[0] <= a < b <= clock[1] + 1e-8, (m["session"], a, b)
                assert covered(a, b, timeline), (m["session"], state, a, b)
                if previous and case["dandiset"] == "000939":
                    assert all(intersection(a, b, x, y) < 1e-8 for x, y in previous["opto_exclusion"])
                    assert covered(a, b, sorted_intervals(previous["context_support"]))
                elif previous and case["dandiset"] == "000056":
                    assert covered(a, b, sorted_intervals(previous["source_intervals"]))
                    assert covered(a, b, sorted_intervals(previous["context_support"]))
            result[context, area, state] = f["duration_s"]
    if case["case"] == "B":
        adn = selected["REST", "ADN", "reference"]
        posub = selected["REST", "PoSub", "reference"]
        for state in adn:
            assert adn[state][0]["intervals_s"] == posub[state][0]["intervals_s"]
    return dict(clock=clock, fields={k: len(x) for k, x in raw["stages"].items()},
                retained=[dict(context=c, region=r, state=s, seconds=v) for (c, r, s), v in sorted(result.items())],
                behavior_exclusions=(previous or {}).get("behavior_conflicting_scores", []),
                excluded_endpoints=(previous or {}).get("invalid_score_rows", []),
                source_scores_overlap=(previous or {}).get("paper_sleep_score", {}).get("overlap_seconds"))


def gaps():
    early = read("manifests/sleep_sources.json")["assets"]
    phase10 = []
    for name in ("000056", "000939", "brainstem"):
        phase10.extend(read(TABLES + name + ".json"))
    early_checked = [r for r in early if "global_spike_gaps_over_10s" in r]
    late_checked = [r for r in phase10 if r.get("eligible", {}).get("REST")]
    assert not any(r["global_spike_gaps_over_10s"] for r in early_checked)
    assert not any(r["gaps"] for r in late_checked)
    h23 = read("manifests/phase9_sources.json")["hc-23"]["per_session_qc"]
    eligible = [r for r in h23 if r["status"] == "PARTIAL"]
    more = [dict(subject=r["subject"], session=r["session"], retained_bouts=r["quality"]["retained_bouts"],
                 support_segments=r["total_support_segments"])
            for r in eligible if r["total_support_segments"] > r["quality"]["retained_bouts"]]
    h24 = [(r["record"]["session"], r["record"]["region"], r["fit"]["max_silence_s"])
           for p in ("pilot", "post-run") for r, _ in entries(f"results/sleep_crcns/hc-24/{p}.json")
           if r["record"]["variant"] == "reference"]
    return dict(phase8_no_trigger_n=len(early_checked), phase10_no_trigger_n=len(late_checked),
                hc23_eligible=len(eligible), hc23_split=more, hc24_preserved=[
                    dict(session=s, region=r, silence_seconds=g) for s, r, g in h24],
                limitation="独立采集时钟可验证 recorded bounds。仅靠 .t 尖峰导出无法判定这些 hc-23 切点的采集是否真实中断。")


def build():
    chosen = read(OUT + "selected_sessions.json")["cases"]
    nwb = {r["case"]: r for r in read(OUT + "nwb_source_recheck.json")}
    previous = {r["asset_id"]: r for d in ("000056", "000939", "brainstem") for r in read(TABLES + d + ".json")}
    docs = [*FIT_FILES, OUT + "matlab-v1.json", OUT + "matlab-qc.json"]
    if (ROOT / (OUT + "matlab-shorter.json")).exists():
        docs.append(OUT + "matlab-shorter.json")
    indexed = defaultdict(lambda: defaultdict(dict))
    for file in docs:
        for row, source in entries(file):
            m = row["record"]
            if m["context"] == "REST_QUIET_EXPLORATORY":
                continue
            ctx, r, v = (m[k] for k in ("context", "region", "variant"))
            s = m["state"]
            key = m["asset_id"]
            if v == "reference" and row.get("long_reach"):
                indexed[key][ctx, r, "AR10"][s] = m, row["long_reach"], source
            indexed[key][ctx, r, v][s] = m, row["fit"], source
    records = []
    legacy = next(x for x in read("manifests/mouse-preparation.json")["sessions"] if x["session"] == "YMV02")
    earlier_nwb = {r["asset_id"]: r for r in read("manifests/sleep_sources.json")["assets"]}
    brainstem_path = ROOT / "data/fixtures/phase11/shorter/brainstem_source_recheck.json"
    brainstem_sorting = {x["case"]: x for x in read(brainstem_path)} if brainstem_path.exists() else {}
    for case in chosen:
        k, asset = case["case"], case["asset_id"]
        meta = previous.get(asset)
        if meta and k in nwb:
            assert meta["digest"] == nwb[k]["digest"]
        original = indexed[asset]
        subsets = {p: fits for p, fits in original.items() if p[1] in case["regions"]}
        summary = [dict(context=c, region=r, variant=v, **state_summary(fits)) for (c, r, v), fits in sorted(subsets.items())]
        if k in nwb:
            clock = scores(case, subsets, nwb[k], meta)
        else:
            clock = (meta or {}).get("clock")
        quality = []
        for region in case["regions"]:
            original_units = next((m for m in summary if m["region"] == region and m["variant"] == "reference"
                                   and m["context"] != "REST_PAPER"), None)
            fit_units = next(iter(original_units["states"].values()))["units"] if original_units else 0
            if k == "I":
                fit_units = legacy["units"]
                for methods in summary:
                    for f in methods["states"].values():
                        f["units"] = fit_units
            units = [u for u in (meta or {}).get("units", []) if u["region"] == region]
            discarded = Counter(reject_reason(u) if "anatomy_agrees" in u else "精确时钟／排序校验未通过"
                                for u in units if not u["eligible"])
            sorting = brainstem_sorting.get(k)
            if sorting:
                assert sorting["accepted_original"] == fit_units
                counts = sorting["unit_count_by_sort_label"]
                discarded.update({f"作者排序 {label}": count for label,count in counts.items() if label != "good"})
                regional_n = sum(counts.values())
                if extra := counts["good"] - fit_units - sum(not u["eligible"] for u in units):
                    assert extra > 0
                    discarded["原参考的共同稳定单元门槛"] = extra
            elif k == "H":
                regional_n = earlier_nwb[asset]["unit_location_labels"].get(region, nwb[k]["unitCount"])
                if fit_units < regional_n:
                    discarded["原区域／单元质控；逐单元原因未另存"] = regional_n - fit_units
            else:
                regional_n = len(units) if units else fit_units
            quality.append(dict(region=region, reference_units=fit_units, regional_audited_units=regional_n,
                                regional_rejected=dict(discarded), raw_file_units=nwb.get(k, {}).get("unitCount"),
                                raw_metrics={f: v["description"] for f, v in nwb.get(k, {}).get("metrics", {}).items()}))
        records.append(dict(case=k, session=case["session"], dataset=case["dandiset"], asset=asset,
                            quality=quality, source_clock=clock, fits=summary))
    old = read(V1_FILE)
    original = {(r["session"], r["state"]): r["d2_bits_per_second"] for r in old["primary"]}
    for row, _ in entries(OUT + "matlab-v1.json"):
        r, f = row["record"], row["fit"]
        assert abs(original[r["session"], r["state"].upper()] - f["d2"]) < 1e-10
    return dict(cases=records, gap_rule=gaps(), window_audit=read(OUT + "000939_equal_duration.json")["summary"])


def thresholds(record):
    by_region = defaultdict(dict)
    for f in record["fits"]:
        if f["context"] == "REST_PAPER":
            continue
        by_region[f["region"]][f["variant"]] = f
    return by_region


def contrast_changes(a, b):
    return [k for k in a.keys() & b.keys() if (a[k] > 0) != (b[k] > 0) or (a[k] < 0) != (b[k] < 0)]


def impact(doc):
    assessed = []

    def add(concern, ds, session, orig, audit, inference, status, file, region=None, category=None):
        assessed.append(dict(concern=concern, dataset=ds, session=session, original=orig, audited=audit,
                             inference=inference, status=status, region=region, evidence=ATLAS + file, category=category))

    for case in doc["cases"]:
        name, ds = case["session"], case["dataset"]
        for region, methods in thresholds(case).items():
            refs = methods["reference"]["effects"]
            strict = {k: m["effects"] for k, m in methods.items() if k in
                      {"strict_units", "lratio_q50", "lratio_q75", "contamination_q50", "contamination_q75"}}
            differences = sorted({d for eff in strict.values() for d in contrast_changes(refs, eff)})
            if strict:
                add("单元质量敏感性", ds, name, refs, strict,
                    "所试子集改变 " + "、".join(differences) + " 的方向。子集是敏感性对照，保留原参考。"
                    if differences else "可比较的差值方向保持。更严格子集对幅度的影响见逐态结果。",
                    "CONFIRMED_IMPACT" if differences else "NO_OBSERVED_IMPACT",
                    OUT + "audit_comparisons.json", region, "sensitivity")
            elif case["case"] in ("F", "G", "H", "I", "B"):
                add("独立单元质量子集", ds, name, "已核验参考群体",
                    "无独立可复算子集或既有阈值保留全组", "无法据此排除未测试的排序敏感性。",
                    "INCONCLUSIVE", OUT + "audit_comparisons.json", region, "insufficient")
            if "AR10" in methods:
                alternative = methods["AR10"]["effects"]
                changes = contrast_changes(refs, alternative)
                add("AR 阶数与残差", ds, name, refs, alternative,
                    "模型阶数改变 " + "、".join(sorted(changes)) + " 的符号，模型拟合仍稳定。"
                    if changes else "在两种阶数下方向保持。残差偏离是模型约束提示。",
                    "CONFIRMED_IMPACT" if changes else "NO_OBSERVED_IMPACT",
                    OUT + "audit_comparisons.json", region, "model")
    d = next(c for c in doc["cases"] if c["case"] == "D")
    add("采集时钟与分期冲突", "000939", d["session"], "Phase 10 整场排除",
        d["source_clock"], "宽带和 LFP 末端一致，尖峰超过两者；保持整场排除。",
        "CONFIRMED_IMPACT", OUT + "nwb_source_recheck.json", "PoSub", "exclusion")
    for case in (x for x in doc["cases"] if x["case"] in "ABCE"):
        add("状态标注及行为边界", case["dataset"], case["session"], "原主要状态区间",
            case["source_clock"], "抽查的已纳入区间落在已核验的 LFP、作者相应分期和行为支持内。未更改原区间。",
            "NO_OBSERVED_IMPACT", OUT + "nwb_source_recheck.json", category="metadata")
    add("整窗等时纳入", "000939", "29 只动物", "参考 Wake／NREM 29 只", doc["window_audit"],
        "60 秒整窗门槛只保留 15 只，缺口均在已有 REM 的连续窗。两组参考方向仍为正。",
        "CONFIRMED_IMPACT", OUT + "000939_equal_duration.json", "PoSub", "selection")
    c = next(x for x in doc["cases"] if x["case"] == "C")
    c_models = thresholds(c)["PoSub"]
    wide, short = [c_models[f"window_{w}s_180s"]["effects"] for w in (60, 45)]
    add("固定 180 秒的窗长敏感性", "000939", c["session"], wide, short,
        "两种窗长给出相同方向。NREM−Wake 的平均差对窗长与选窗位置敏感。",
        "CONFIRMED_IMPACT", OUT + "audit_comparisons.json", "PoSub", "sensitivity")
    e = next(x for x in doc["cases"] if x["case"] == "E")
    e_models = thresholds(e)["PoSub"]
    if "window_45s_180s" in e_models:
        add("较短整窗是否恢复样本", "000939", e["session"], "60 秒整窗未纳入",
            e_models["window_45s_180s"]["effects"],
            "45 秒窗在这场记录构成三态 180 秒，可作探索估计。原 60 秒队列仍独立。",
            "CONFIRMED_IMPACT", OUT + "audit_comparisons.json", "PoSub", "selection")
    gap = doc["gap_rule"]
    add("静默切点", "hc-23", f"{gap['hc23_eligible']} 场合格记录", "原处理在 >10 s 群体尖峰静默处切段",
        gap["hc23_split"], "确认进入中间支持区间的切点。是否为采集断档及对 d₂ 的实际影响未证实。",
        "INCONCLUSIVE", OUT + "audit_comparisons.json", "MEC II/III", "insufficient")
    add("长零计数与连续记录", "hc-24", "跑后 CA1", "按作者 NREM 与 LFP 边界保留",
        gap["hc24_preserved"], "最大 55.28 秒的连续零计数仍进入拟合。长静默本身不足以判定断档。",
        "NO_OBSERVED_IMPACT", OUT + "audit_comparisons.json", "CA1", "metadata")
    add("缺失 REM", "PNAS-V1", "原公开 13 场", "只报告 Wake／NREM", "Phase 11 V1 回放仍为两态",
        "该数据无法回答此处的 REM−NREM。", "NOT_APPLICABLE", V1_FILE, "V1", "coverage")
    return assessed


if __name__ == "__main__":
    results = build()
    write(ROOT / (OUT + "audit_comparisons.json"), results)
    statuses = impact(results)
    write(ROOT / (OUT + "impact_assessment.json"), statuses)
    print(Counter(r["status"] for r in statuses))
    print(results["gap_rule"]["hc23_split"])
