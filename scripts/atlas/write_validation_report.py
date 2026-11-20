"""Render the verified Phase 11 comparisons as a Chinese audit report."""
from build_audit_results import OUT, thresholds
from build_atlas_data import ROOT, read


def fmt(n, digits=3):
    return "—" if n is None else f"{n:+.{digits}f}"


def pair(effects):
    return "／".join(fmt(effects.get(x)) for x in ("NW", "RN"))


SUBSETS = {
    "reference": "参考 AR(5)", "strict_units": "原严格 ISI",
    "lratio_q75": "L-ratio q75", "lratio_q50": "L-ratio q50",
    "contamination_q75": "污染率 q75", "contamination_q50": "污染率 q50",
}


def make_report():
    data = read(OUT + "audit_comparisons.json")
    cases = {c["case"]: c for c in data["cases"]}
    windows = data["window_audit"]
    d = cases["D"]["source_clock"]
    b = thresholds(cases["B"])["PoSub"]
    wide, shorter = [thresholds(cases["C"])["PoSub"][f"window_{w}s_180s"]["effects"] for w in (60, 45)]
    e = thresholds(cases["E"])["PoSub"].get("window_45s_180s")
    v1 = thresholds(cases["I"])["V1"]
    ldt = thresholds(cases["G"])["LDT"]
    runs = [read(OUT + f"matlab-{x}.receipt.json") for x in ("v1", "qc", "shorter")
            if (ROOT / (OUT + f"matlab-{x}.receipt.json")).exists()]
    results_file = "../results/validation/audit_comparisons.json"
    lines = [f"""# Phase 11 针对性验证

按归档路径、已通过质量门槛的单元数、状态覆盖和已报告的 LDT 例外冻结九场样本，再求补充估计。所有新增拟合均用 Sooter 工具箱 f26c782dca3b1c9ab18cc92c6449591c44e0d507、MATLAB R2024b、40 ms、Yule–Walker 和固定 β=2。共有 {sum(r['count'] for r in runs)} 条输入序列，V1 参考回放另拟合 AR(10)。历史分期、计数、拟合及报告保留。[抽样清单](../results/validation/selected_sessions.json) · [逐态数据与诊断]({results_file})

## 实际影响

| 关注点 | 数据与会话 | 原结果 | 核验结果 | 对推断的影响 | 状态 |
|---|---|---|---|---|---|
| 单元选择 | 000056，Mouse28-140313 PoSub | NW {fmt(b['reference']['effects']['NW'], 5)} | q75 {fmt(b['lratio_q75']['effects']['NW'], 5)}、q50 {fmt(b['lratio_q50']['effects']['NW'], 5)} | 单会话的微弱正向变为近零或小幅负向 | CONFIRMED_IMPACT，[曲线](../figures/phase11/sorting_sensitivity.pdf) |
| 脑干排序 | DL108_122320n1 LDT | NW {fmt(ldt['reference']['effects']['NW'])} | 污染率 q50 NW {fmt(ldt['contamination_q50']['effects']['NW'])} | 唯一一只 LDT 的 NW 对单元选择敏感，RN 仍为正 | CONFIRMED_IMPACT，[曲线](../figures/phase11/sorting_sensitivity.pdf) |
| 阶数敏感 | 000056，Mouse28-140313 PoSub | RW {fmt(b['reference']['effects']['RW'])} | AR(10) RW {fmt(b['AR10']['effects']['RW'])} | REM−Wake 变号，两阶数值仍稳定 | CONFIRMED_IMPACT，[逐拟合]({results_file}) |
| 时钟冲突 | 000939，A3705／200306 | 整场排除 | 宽带和 LFP 同在 {d['clock'][1]:.3f} 秒截止，尖峰到 {d['latest_spike']:.5f} 秒 | {len(d['out_of_range'])}/{d['stages']} 个评分行越界。维持排除 | CONFIRMED_IMPACT，[复核](../results/validation/nwb_source_recheck.json) |
| 整窗样本 | 000939，29 只 | 参考 NW 为 29/29 正向 | 60 秒等时分析保留 15 只，另 14 只均卡在 REM | 子集参考均值偏高，两个组的方向一致 | CONFIRMED_IMPACT，[逐动物窗数](../results/validation/000939_equal_duration.json) |
| 窗长敏感 | 000939，A6216／210623 | 3 × 60 秒 NW {fmt(wide['NW'])} | 4 × 45 秒 NW {fmt(shorter['NW'])} | 同时长下方向保持，NW 幅度约缩小 30% | CONFIRMED_IMPACT，[逐拟合](../results/validation/matlab-qc.json) |
| 区间一致性 | 000056 的 A、B 和 000939 的 C、E | 原参考区间 | 落在独立采集时钟、作者评分和核验过的行为边界内 | 原区间继续使用 | NO_OBSERVED_IMPACT，[逐区间]({results_file}) |
| 原估计回放 | 小鼠 V1，YMV02 | NW {fmt(v1['reference']['effects']['NW'])} | AR(5) 逐态精确重放。AR(10) NW {fmt(v1['AR10']['effects']['NW'])} | 方向保持，幅度对阶数敏感 | NO_OBSERVED_IMPACT（回放及方向），[实测](../results/validation/matlab-v1.json) |
| 群体长静默 | hc-23，21 场 | 按尖峰支持切段 | 8 场支持段数超过保留的作者 bout 数 | 预处理切段生效。采集连续性与 d₂ 变化仍缺独立证据 | INCONCLUSIVE，[切段证据]({results_file}) |
| REM 缺失 | 小鼠 V1，公开 13 场 | 仅 Wake 和 NREM | 本轮回放仍为两态 | 无法评估该公开复现的 REM 差值 | NOT_APPLICABLE，[既有结果](02_mouse_replication.md) |

NW 为 NREM−Wake，RN 为 REM−NREM，RW 为 REM−Wake。估计和差值单位为 bits/s。状态判断对应表中具体问题。抽样结论的适用范围为所列会话和所试规格。完整状态索引保存在 [影响表](../results/validation/impact_assessment.json)。

## A · 排序与解剖

000056 的 pE、pI 是由波形和放电特征推定的细胞类型。L-ratio 和 Mahalanobis isolation distance 表征聚类分离。作者 isi_index 为不应期违反率。in_lab_cell_list 标记进入实验室的分析名单。源文件说明原论文的逐单元名单无法完整恢复。000939 的三场没有这些隔离指标，原严格 ISI 使用从排序时钟重新计算的不应期违反率。000552 的 putative_cell_type 是 CellExplorer 分类。脑干检查作者 good、mua、noise 标签及通道核团，检查原单位 ID 与 1 ms 尖峰哈希。ContamPct 是排序表中的污染估计。[原 NWB 字段与时钟](../results/validation/nwb_source_recheck.json) · [脑干原排序表收据](../data/fixtures/phase11/shorter/brainstem_source_recheck.json)

公开资料缺少可推广的质量阈值。探索对照在原群体内按有限 L-ratio 或 ContamPct 的第 75 或 50 百分位取较低值，跨状态保持同一组单元、原时间箱和原分期。图和数据列出每场单元数量及跨阈值变化。原严格 ISI 和新增百分位作为各自独立对照。[排序曲线](../figures/phase11/sorting_sensitivity.pdf)

| 会话／脑区 | 原区域单元数 | 参考保留 | 归档排除理由与数量 |
|---|---:|---:|---|"""]
    for c in cases.values():
        for q in c["quality"]:
            reasons = "、".join(f"{k} {v}" for k, v in q["regional_rejected"].items()) or "原审计未列区域内排除"
            lines.append(f"| {c['session']} {q['region']} | {q['regional_audited_units'] or '—'} | {q['reference_units']} | {reasons} |")
    lines.append("""
### 单元对照逐态结果

| 会话／脑区 | 群体 | 单元 | Wake | NREM | REM | NW／RN |
|---|---|---:|---:|---:|---:|---:|""")
    for c in cases.values():
        for region, models in thresholds(c).items():
            for key, title in SUBSETS.items():
                if key not in models:
                    continue
                f = models[key]
                states = f["states"]
                vals = " | ".join(fmt(states.get(s, {}).get("d2")).lstrip("+") for s in ("Wake", "NREM", "REM"))
                units = next(iter(states.values()))["units"]
                lines.append(f"| {c['session']} {region} | {title} | {units} | {vals} | {pair(f['effects'])} |")
    included, excluded = (windows[k] for k in ("included", "excluded"))
    e_line = f"A1821／190910 用 4 × 45 秒得到 NW／RN {pair(e['effects'])} bits/s。" if e else ""
    lines.append(f"""
## B · 分期、行为和独立时钟

逐场重新读出公开 NWB 的时钟长度、评分与行为表。000056 先区分原论文评分和后续 Buzcode，同名状态重叠及参考有效暴露在逐态数据中给出。两者为同一批记录的不同评分。A、B 所用区间在休息和原 sub-recording 支持内。B 同步 ADN／PoSub 的每态时间箱一致。C、E 在独立 LFP、作者状态标注和家笼支持内，避开光刺激及其前后 60 秒。D 的原时钟不相容并已排除。全程保留原始录制时间单位。

## C · 等时选择

[逐动物审计](../results/validation/000939_equal_duration.json) 保存 29 只的三态时长、候选窗数、参考单元数和效应。原规则要求所有可用参考状态各至少有三个连续 60 秒窗。14 只未纳入均由于已有 REM 短或碎，Wake 与 NREM 窗均达门槛。

| 60 秒整窗组 | 动物 | 参考 NW | 中位 Wake 秒 | 中位 NREM 秒 | 中位 REM 秒 | 中位单元数 |
|---|---:|---:|---:|---:|---:|---:|
| 纳入 | {included['n']} | {fmt(included['mean_reference'])} | {included['state_median_seconds']['Wake']:.2f} | {included['state_median_seconds']['NREM']:.2f} | {included['state_median_seconds']['REM']:.2f} | {included['median_units']} |
| 排除 | {excluded['n']} | {fmt(excluded['mean_reference'])} | {excluded['state_median_seconds']['Wake']:.2f} | {excluded['state_median_seconds']['NREM']:.2f} | {excluded['state_median_seconds']['REM']:.2f} | {excluded['median_units']} |

两组参考效应各自全为正。原纳入组的等时平均 NW 为 {fmt(windows['matched_60_mean'])}。固定 180 秒和同组单元，先用典型 C 比较 3 × 60 秒和 4 × 45 秒。两规格均收敛并保留 NW 和 RN 方向，RN 为 {fmt(wide['RN'])} 与 {fmt(shorter['RN'])}。通过这场方向门槛后才用较短窗核验排除个体 E。{e_line} 固定 180 秒的 45 秒候选窗规则理论上有 {windows['eligible_at_45_same_180s']} 只达标。较短窗的幅度敏感，保持原 60 秒规格和独立样本口径。

## D · AR 诊断

| 会话／脑区 | AR(5) NW／RN | AR(10) NW／RN | AR(5) 最大绝对残差 ACF | AR(10) 最大绝对残差 ACF | 两阶最小 AR 根 |
|---|---:|---:|---:|---:|---:|""")
    for c in cases.values():
        for region, models in thresholds(c).items():
            if "AR10" not in models:
                continue
            vals = [models[x]["states"] for x in ("reference", "AR10")]
            residual = [max(abs(x) for f in group.values() for x in f["residual_acf"] if x is not None) for group in vals]
            root = min(f["min_root"] for group in vals for f in group.values())
            lines.append(f"| {c['session']} {region} | {pair(models['reference']['effects'])} | {pair(models['AR10']['effects'])} | {residual[0]:.3f} | {residual[1]:.3f} | {root:.4f} |")
    dgap = data["gap_rule"]
    lines.append(f"""
补充估计的求解正常且 AR 根在单位圆外。残差峰度、残差 ACF、观测和模型 ACF 均保存在 [逐态诊断]({results_file})。YMV02 Wake 的 ACF RMSE 从约 0.137 降至 0.047。残差峰度与相关性阈值描述模型偏离，数值收敛和 AR 稳定性是另一组判断。LDT 的 RN 在两个阶数下均为正，单只证据保留为探索性。β=2 固定参考模型流形，AR 系数逐会话和状态估计，生理相变仍需独立实证。

## E · 连续性与长静默

Phase 8 逐会话有完整静默记录的 {dgap['phase8_no_trigger_n']} 场，以及 Phase 10 的 {dgap['phase10_no_trigger_n']} 个主要合格记录资产，原群体审计均无大于 10 秒的峰间隙触发。hc-23 的 8/21 场出现额外的支持切段。例如 Rat-26 Session2 从 54 个保留作者 bout 形成 65 个支持段。尖峰时钟经原 Neuralynx .t 核验，但这些切点缺少能够证明采集是否实际中断的独立连续记录证据。恢复切点后 d₂ 会改变多少仍不确定。hc-24 的跑后 CA1 最大连续零计数 **55.28 秒**保留在拟合中。继续用独立时钟界定采集有效范围，并保留原预处理。

## 输入与执行

[抽样脚本](../scripts/atlas/select_validation.py) · [公开源复核](../scripts/atlas/recheck_sources.py) · [固定时长重选窗](../scripts/atlas/build_qc_fixture.py) · [Colab 与本机计数逐数组核验](../results/validation/colab_preparation_receipt.json) · [完整影响表](../results/validation/impact_assessment.json)
""")
    for r in runs:
        lines.append(f"- [{r['action'].split('/')[-1]}]({r['action']}) · {r['count']} 条输入 · 提交 {r['commit']} · 输入 SHA-256 {r['input_sha256']}。")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    (ROOT / "reports/26_targeted_validation.md").write_text(make_report(), encoding="utf-8")
