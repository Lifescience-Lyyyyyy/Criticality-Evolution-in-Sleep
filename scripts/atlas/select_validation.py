#!/usr/bin/env python3
"""Freeze targeted audit sessions using source coverage only."""
from collections import defaultdict

from build_atlas_data import ROOT, read, write
from catalog import FIT_FILES, SOURCE_COMMIT, V1_FILE


def select():
    series = defaultdict(dict)
    for path in FIT_FILES:
        rows = read(path)["records"]
        for row in rows if isinstance(rows, list) else [rows]:
            r = row["record"]
            if r["context"] in {"REST", "HOME"}:
                key = (r["dandiset"], r["region"], r["asset_id"], r["variant"])
                series[key][r["state"]] = r

    def candidates(dataset, region, states):
        out = []
        for (ds, area, _, variant), rows in series.items():
            if (ds, area, variant) == (dataset, region, "reference") and set(states) <= rows.keys():
                row = rows["NREM"]
                out.append((len(row["unit_ids"]), row["path"], row))
        return sorted(out, key=lambda a: a[:2])

    def best(xs):
        return min(xs, key=lambda p: (-p[0], p[1]))[-1]

    def item(code, row, rule, regions=None):
        return {**{k: row[k] for k in ("dandiset", "subject", "session", "asset_id", "path")},
                "case": code, "regions": regions or [row["region"]], "selection": rule}

    all_states = ("Wake", "NREM", "REM")
    adn = candidates("000056", "ADN", all_states)
    post = candidates("000056", "PoSub", all_states)
    posub = candidates("000939", "PoSub", all_states)
    audited = read("results/phase10/audit/000939.json")
    rejected = min((r for r in audited if r["gate"] == "BLOCKED_QUALITY" and r["clock"]["seconds"]),
                   key=lambda r: r["path"])
    adn_units = {r["asset_id"]: n for n, _, r in adn}
    common = [(min(n, adn_units[r["asset_id"]]), path, r) for n, path, r in post if r["asset_id"] in adn_units]
    excluded = [r for _, _, r in candidates("000939", "PoSub", ("Wake", "NREM"))
                if not {"Wake", "NREM"} <= series["000939", "PoSub", r["asset_id"], "equal_duration"].keys()]
    legacy = candidates("000552", "CA1", all_states)
    ldt = candidates("zenodo:19462601", "LDT", all_states)
    assert len(ldt) == 1
    v1 = min({r["session"] for r in read(V1_FILE)["primary"]})
    cases = [
        item("A", best(adn), "取三态 ADN 单元数最多的会话，平局取归档路径字典序。"),
        item("B", best(common), "只取同步三态 ADN／PoSub。最大化两区较小的单元数，平局取归档路径字典序。", ["ADN", "PoSub"]),
        item("C", posub[len(posub)//2][-1], "000939 三态合格会话按单元数和归档路径排序，取中位。"),
        item("D", rejected, "有可核验采集时钟且因质量整场排除的 000939 会话，取归档路径字典序第一场。", ["PoSub"]),
        item("E", min(excluded, key=lambda r: r["path"]), "参考有 Wake／NREM、60 秒整窗等时长未纳入，取归档路径字典序第一场。"),
        item("F", best(candidates("zenodo:19462601", "PAG", all_states)), "合格 WT PAG 的三态会话中取单元数最多的一场。"),
        item("G", ldt[0][-1], "复核 Phase 10 预先已报告的 LDT 例外。只存在一个合格动物，无按新结果选择。"),
        item("H", legacy[len(legacy)//2][-1], "小鼠 CA1 三态会话按单元数和归档路径排序，取中位。"),
        dict(case="I", dandiset="PNAS-V1", subject=None, session=v1, asset_id=v1 + ".mat", path=v1 + ".mat",
             regions=["V1"], selection="原公开 V1 会话按名称取首场，补足 AR 诊断覆盖。"),
    ]
    assert len({(r["dandiset"], r["asset_id"]) for r in cases}) == len(cases) == 9
    return {"baseline": SOURCE_COMMIT, "rule": "仅使用已发布的会话覆盖、单元数与既有例外标记。冻结样本后才计算本阶段对照。",
            "cases": cases}


if __name__ == "__main__":
    sample = select()
    write(ROOT / "results/validation/selected_sessions.json", sample)
    for row in sample["cases"]:
        print(row["case"], row["dandiset"], row["session"], "/".join(row["regions"]))
