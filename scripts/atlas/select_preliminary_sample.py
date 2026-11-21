"""Freeze additional quality/timing cases from unit counts and coverage alone."""
from collections import defaultdict

from build_atlas_data import ROOT, read, write
from build_pairwise_fixture import OUT
from catalog import FIT_FILES

BASELINE = "7f100ec04da3ae0725af76da1901e9188aa15f38"
GROUPS = [
    ("000056","ADN",2,"lratio_q50"),
    ("000056","PoSub",1,"lratio_q50"),
    ("000939","PoSub",3,"isi_q50"),
    ("zenodo:19462601","PAG",2,"contamination_q50"),
    ("000552","CA1",1,"isi_q50"),
]


def select():
    excluded = {(r["dandiset"],r["subject"]) for r in read("results/validation/selected_sessions.json")["cases"]}
    series = defaultdict(dict)
    for file in FIT_FILES:
        rows = read(file)["records"]
        for row in rows if isinstance(rows,list) else [rows]:
            r, fit = row["record"], row["fit"]
            if r["variant"] != "reference" or r["context"] not in ("REST","HOME"):
                continue
            if (r["dandiset"],r["subject"]) in excluded:
                continue
            spans = fit["segment_bins"]
            series[r["dandiset"],r["region"],r["asset_id"]][r["state"]] = (
                {k:r[k] for k in ("dandiset","subject","session","region","asset_id","path","context")},
                len(r["unit_ids"]), fit["duration_s"], max(spans) if isinstance(spans,list) else spans,
            )
    cases=[]
    for ds,region,n,method in GROUPS:
        animals={}
        for (source,area,_),states in series.items():
            if (source,area)!=(ds,region) or not {"Wake","NREM","REM"} <= states.keys():
                continue
            m, units, _, _ = states["NREM"]
            if units < 10:
                continue
            lengths = {s: a[3] for s,a in states.items()}
            durations = {s: a[2] for s,a in states.items()}
            continuous = min(lengths[s] for s in ("Wake","NREM"))>=6000
            row = dict(m, original_units=units, seconds=durations,max_continuous_bins=lengths,
                       supports_120s_halves=continuous, method=method)
            score = (-sum(durations.values()), m["path"])
            key = m["subject"]
            if key not in animals or score < animals[key][0]:
                animals[key] = score,row
        ordered = sorted((v[1] for v in animals.values()),key=lambda r:(r["original_units"],r["path"]))
        assert len(ordered)>=n,(ds,region,len(ordered))
        indices = [len(ordered)//2] if n==1 else [0,len(ordered)-1] if n==2 else [0,len(ordered)//2,len(ordered)-1]
        picked = [ordered[i] for i in indices]
        if region in ("ADN","PAG") or ds=="000939":
            qualifying = [r for r in picked if r["supports_120s_halves"]]
            if not qualifying:
                others = [r for r in ordered if r["supports_120s_halves"] and r not in picked]
                if others:
                    picked[-1] = others[-1]
                    qualifying = [picked[-1]]
            if qualifying:
                target = max(qualifying,key=lambda r:(min(r["max_continuous_bins"][s] for s in ("Wake","NREM")),r["path"]))
                target["temporal_split"] = True
        for row in picked:
            row["source_recheck"] = row is picked[0] or row.get("temporal_split",False)
            cases.append(row)
    assert len(cases)==9
    return dict(baseline=BASELINE,selection=[
        "排除 Phase 11 原抽样的同数据集动物。每只动物优先三态有效时长总和最大的会话，同长按归档路径。",
        "至少 10 个共同单元。按单元数和路径取低、高或中位，范围覆盖两个或三个不同动物。",
        "每个目标组中，若已选会话均缺合格的 Wake、NREM 240 秒连续段，且另有合格动物，以最高单元候选替换高端样本。",
        "每场只试一个排序子集。原始群体内取低侧中位数并保留并列。L-ratio 与 ContamPct 来自发布字段，其余使用已发表的 1.5 ms ISI 检查定义计算相对排序。该规则为抽样敏感性规格。",
        "同会话分别选用 Wake 与 NREM 各自最长的连续段，等长取该段前 120 秒和后 120 秒；平局取最早段。"
    ],cases=cases)


if __name__=="__main__":
    result=select()
    write(ROOT/(OUT+"sorting_checks/selection.json"),result)
    for c in result["cases"]:
        print(c["dandiset"],c["region"],c["subject"],c["session"],c["original_units"],
              c.get("temporal_split",False),c["source_recheck"])
