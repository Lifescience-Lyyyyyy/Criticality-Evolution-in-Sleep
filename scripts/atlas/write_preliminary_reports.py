"""Render the Phase 11.1 scientific reports from measured and audited outputs."""
from build_atlas_data import ROOT, read
from build_pairwise_fixture import OUT

LINK="../"+OUT
PAIRS={"NW":"NREM−Wake","RN":"REM−NREM","RW":"REM−Wake"}


def fmt(x,d=3):
    return "—" if x is None else f"{x:+.{d}f}"


def seq(it):
    return "／".join(fmt(it.get(k)) for k in ("NW","RN"))


def assessment(stats,unit):
    if not stats or not stats["n"]:
        return "—"
    return f"{fmt(stats['mean'])}（{stats['positive']}正、{stats['negative']}负／{stats['n']}{unit}）"


def write_validation():
    pairing=read(OUT+"pairwise_equal_duration/summary.json")
    independence=read(OUT+"independence.json")
    sorting=read(OUT+"sorting_checks/summary.json")
    temporal=read(OUT+"temporal_checks/within_state.json")["cases"]
    sources=read(OUT+"temporal_checks/source_clock.json")
    hc23=read(OUT+"temporal_checks/hc23_gap.json")
    model=read(OUT+"model_checks/adequacy.json")["cohorts"]
    receipt=sorting["provenance"]
    text=[f"""# Phase 11.1 初步结果验证

以 Phase 10 参考序列和 Phase 11 复核作为基线，仓库冻结在 {independence['baseline']}。新增两批均使用原 Sooter 提交 f26c782dca3b1c9ab18cc92c6449591c44e0d507、MATLAB R2024b、40 ms 分箱、AR(5)、固定 β=2 和 Yule–Walker。AR(10) 只作既有模型敏感性比较。原数值、三态控制和已确认排除均保留原名称及来源。

## A · 000939 PoSub 成对等时

同场按指定两态各取三段完整 60 秒窗。每段来自同一连续有效作者分期。按开始时间取第一、靠前的中位和最后一个完整窗，每态共 180 秒，各态使用同一组已排序单元。三态固定等时另列，要求三态同时符合门槛。旧等时估计保留原变量时长。均差单位为 bits/s。区间以动物有放回抽样 10,000 次，种子 11101。

| 对比 | 原等时动物 | 新 动物／会话 | 均值／中位数 | 正／负 | 动物 bootstrap 95% | 新减旧全样本均差 | 共同动物均差变化 |
|---|---:|---:|---:|---:|---|---:|---:|"""]
    for k,r in pairing["contrasts"].items():
        n=r["corrected"]
        ci="、".join(fmt(x) for x in n["bootstrap95"])
        text.append(f"| {PAIRS[k]} | {r['legacy']['n']} | {n['n']}／{n['sessions']} | {fmt(n['mean'])}／{fmt(n['median'])} | {n['positive']}／{n['negative']} | [{ci}] | {fmt(r['corrected_minus_legacy_all_animals'])} | {fmt(r['corrected_minus_legacy_same_animals'])} |")
    text.append(f"""
原 NW 等时受到其他状态门槛限制，14 只仅因 REM 完整窗不足被排除。新 NW 允许这 14 只纳入。新 RN／RW 各排除 18 只，其中 4 只无有效 REM，14 只 REM 够参考分析但少于三窗。固定三态对照保留 11 只，其三种配对统计另存。NW 在新规则下仍以正向为主，同时出现两只负向。RN 维持负向，RW 有混合方向。跨样本均值变化可能同时来自样本组成与配平规格，共同动物差异单独给出。[逐动物排除与窗数]({LINK}pairwise_equal_duration/eligibility.json) · [结果与保留会话]({LINK}pairwise_equal_duration/summary.json) · [MATLAB {pairing['provenance']['action'].split('/')[-1]}]({pairing['provenance']['action']})

## B · 统计独立性

对 {len(independence['cohorts'])} 个区域／场景记录表逐会话求差，再重组为动物内平均。覆盖 {len(independence['source_metadata'])} 个已分析的哺乳动物睡眠数据源，逐对核验 {len(independence['inter_dataset_checks'])} 组资产 ID、归档摘要、动物／会话编码和源论文。未找到同一原始记录跨来源计入的证据。未披露的原动物重用仍无法排除。

000061 的 PRE／POST 与同步 CA1／BLA 共享 Rat08、Rat10、Rat11。000056 的 ADN／PoSub 共享 Mouse24、Mouse25、Mouse28，且两份评分属于同一记录。脑干同一次探针记录的核团共享动物。原小鼠 V1 仅有 13 个可核实会话，独立动物数未知，未将另一个归档映射 000166 重复计入。不同物种与两态缺失分别保留。历史动物区间和本轮 bootstrap 使用动物层级；技术重复百分位不代表动物区间。[各数据源身份与复算结果]({LINK}independence.json)

## C · 最近一轮排序敏感性样本

选择按归档路径、单元数、动物身份、有效时长和连续段覆盖冻结，排除 Phase 11 同来源已抽样动物。本轮每场固定一个较低质量指标的中位数子集，保留并列。L-ratio 与 ContamPct 来自作者排序字段。000939、CA1 及无有限 L-ratio 的场次用短于 1.5 ms 的 ISI 比例给原群体内单元排序，分数来自有效作者状态区间内的同段尖峰。细胞类型不参与质量排序。表内 NW／RN 均为 bits/s。方向改变或绝对效应变化达到 0.25 bits/s 记为样本内影响，0.25 仅是本表的描述标记。

| 数据／脑区 | 会话 | 质量指标 | 原／子集单元 | 原 NW／RN | 新 NW／RN | 新减原 NW／RN | 状态 |
|---|---|---|---:|---:|---:|---:|---|""")
    for c in sorting["cases"]:
        status=c["assessment"]
        if c.get("direction_changes"):
            status+="（"+"、".join(c["direction_changes"])+"变号）"
        text.append(f"| {c['dandiset']} {c['region']} | {c['session']} | {c['method']} | {c['original_units']}／{c['alternate_units']} | {seq(c['reference'])} | {seq(c['alternate'] or {})} | {seq(c.get('delta',{}))} | {status} |")
    rw_changes=[f"{c['session']}（{fmt(c['reference']['RW'])} → {fmt(c['alternate']['RW'])} bits/s）"
                for c in sorting["cases"] if "RW" in c.get("direction_changes", [])]
    if rw_changes:
        text.append("\nREM−Wake 变号的会话为" + "、".join(rw_changes) + "。9 个新增会话的 NW 与 RN 都保持原方向。")
    text.append("\nPhase 11 原抽样保留两个已知的 NW 变号，000056 PoSub 的 Mouse28 由 +0.03872 变为 −0.02888 bits/s，LDT 的 DL108 由 +0.23273 变为 −0.14774 bits/s。两个子集均按各自质量指标取 q50。[原审计](26_targeted_validation.md)")
    fallbacks=[c["session"] for c in sorting["cases"] if c.get("fallback_reason")]
    detail="、".join(fallbacks)
    text.append(f"""
{detail} 的源文件没有可用 L-ratio，按预定 ISI 替代口径保留同一已冻结会话。完整数据含每态 d₂、AR 稳定根、求解标志、有效时长、观测／拟合 ACF、残差 ACF 和峰度。该抽样可发现会话级依赖，不能估计整个队列的排序影响率。[冻结样本]({LINK}sorting_checks/selection.json) · [逐单元分数]({LINK}sorting_checks/source_metrics.json) · [逐态数值及诊断]({LINK}sorting_checks/summary.json) · [MATLAB {receipt['action'].split('/')[-1]}]({receipt['action']})

## D · 分期、行为与时间边界

复核 {len(sources)-1} 个未参与 Phase 11 原抽样的 NWB 会话及一个 PAG 会话。NWB 回读公开分期、独立 LFP 或录制偏移与 LED 时钟、有效休息时段。小鼠 CA1 依据该归档会话的家笼描述，光标记阶段前后各留 60 秒。脑干从公开档案仅读取 EEG、EMG 文件前缀，独立复核信号长度，同时核对已下载的作者评分和 1 ms 尖峰。使用的分箱均须处于有效单状态段。允许相应历史分箱的端点舍入。复查在库参考序列的尖峰总数和单元身份，没有改写作者分期。000056 保留原论文与 Buzcode 两个同源场景，000939 已确认两个不兼容时钟会话仍保持排除。[独立时钟逐项记录]({LINK}temporal_checks/source_clock.json)

hc-23 原静默分割切点仍未决。Rat-26 Session2 的该段归档目录共有 {hc23['matching_members']} 个排序尖峰文件，未找到可独立证实该切点的连续采集时钟或对应 LFP。已有 8/21 场在输入处理中发生静默处切段，现有证据无法把切点判为采集停顿，也无法估计恢复这些输入段的数值影响。[归档检查]({LINK}temporal_checks/hc23_gap.json)

## E · 模型适配和同状态时间敏感性

| 新时段对照 | Wake、NREM 各 120 秒前段 NW | 后段 NW | 后减前 | 方向 |
|---|---:|---:|---:|---|""")
    for c in temporal:
        text.append(f"| {c['dandiset']} {c['region']} {c['session']} | {fmt(c['first_effect']['NW'])} | {fmt(c['last_effect']['NW'])} | {fmt(c['delta']['NW'])} | {'改变' if c['direction_changes'] else '保持'} |")
    text.append(f"""
每态先选同场最长的同状态连续段，再分别取不重叠的前后 120 秒。两次估计均沿用同一组单元和原有 AR(5) 规格。两段之间的变化体现该会话对有效时段选取的敏感性，单次保持方向无法证明严格平稳性。[两个时间段的原始拟合]({LINK}temporal_checks/within_state.json)

| 来源与脑区 | AR(5) 通过／拟合 | AR(10) 通过／拟合 | AR(5) 残差 ACF／峰度警示次数 | AR(10) 残差 ACF／峰度警示次数 |
|---|---:|---:|---:|---:|""")
    for c in model:
        if "diagnostic_sample_session" in c:
            text.append(f"| V1 公开会话 | 26／26 | YMV02 两态 | YMV02 单场有完整 ACF | YMV02 单场有完整 ACF |")
            continue
        counts=[]
        warnings=[]
        for s in ("reference","AR10"):
            rows=list(c["state_fit_summaries"][s].values())
            counts.append(f"{sum(r['converged'] for r in rows)}／{sum(r['fits'] for r in rows)}")
            warnings.append(f"{sum(r['residual_acf_warnings'] for r in rows)}／{sum(r['kurtosis_warnings'] for r in rows)}")
        text.append(f"| {c['id'].replace('|','／')} | {' | '.join(counts+warnings)} |")
    text.append(f"""
**数值有效性。** 核对退出状态、求解与 AR 根。**模型充分性。** 分别检查观测／拟合 ACF、残差 ACF、残差峰度、时长和有效单元数。警示次数按拟合计，不表示独立动物，也不用于按方向排除。**生物学解释。** β=2 定义固定参考流形，谱斜率未参与该参数的估计。现有对照无法分别辨识放电率、方差与离散程度、时间相关、同步及状态非平稳性。时段或排序敏感均可作为机制研究的线索，因果贡献仍需独立设计。[完整模型与 AR 阶数结果]({LINK}model_checks/adequacy.json)""")
    (ROOT/"reports/29_preliminary_validation.md").write_text("\n".join(text)+"\n",encoding="utf-8")


def write_summary():
    assessments=read(OUT+"evidence_assessment.json")
    text=["""# Sleep-State-Dependent Proximity to Criticality Across Mammalian Neural Circuits: Preliminary Cross-Regional Evidence

## 可作为初步结果展示的范围

在固定 β=2 的高斯 AR 估计下，几个独立小鼠和大鼠记录来源的 NREM d₂ 往往高于 Wake 和 REM。可报告的是模型依赖、区内配对的重复方向。物理相变、普遍发生比例和具体生理机制需要独立证据。

**复现的数值观察。** 小鼠 000056 ADN、000939 PoSub、000552 CA1 和 WT PAG 的参考组值重复出现 NREM 相对较高的方向。大鼠额叶的方向相似。小鼠 V1 的结论以 13 个公开会话为单位，对 REM 与独立动物数保留缺失标记。

**较稳健的对照。** 独立 PoSub 的新成对等时中，NW 保留 29 只，27 正 2 负。RN 在 11 只中全部为负。不同脑区的组方向在已检验的 AR 阶数、单元数或放电率控制中多有保持。已查阅的有限排序和分期样本不能代替整个队列复核。

**敏感结果。** 个别会话与小样本脑干核团的效应幅度或方向随排序、时段与 AR 阶数变化。REM−Wake 的方向不统一。原等时分析对 REM 的要求限制了 NW 分母，新配对规格与原结果并列。

**未决问题。** hc-23 静默切点缺独立采集连续性验证。小鼠 V1 的生物学重复单位仍未知。部分数据只有 NREM 或缺有效 REM。低移动清醒、两套同源分期和稀疏核团资料各自限制外推。

**机制假设。** 群体率、方差与过离散、时间相关、跨单元同步及状态内非平稳都可能贡献于估计。本项目当前控制尚不足以鉴别这些机制独立的因果作用。

## 各队列的证据评估

下列标签描述证据边界。数值为参考 AR(5) 的均差，单位 bits/s。每列的分母由同一会话中可以配对的状态决定。独立动物数显示为未知时，保留原会话单位。共同动物在同来源分区和 PRE／POST 间重复出现，表行不可累加为生物样本。

| 物种／区域／来源／场景 | 动物／会话 | NREM−Wake | REM−NREM | REM−Wake | 证据评估 |
|---|---:|---:|---:|---:|---:|"""]
    for r in assessments:
        n=r["animals"]
        unit="场" if n is None else "只"
        label="／".join(r["id"].split("|"))
        vals=[assessment(r["paired"].get(k),unit) for k in PAIRS]
        text.append(f"| {label} | {n if n is not None else '未知'}／{r['sessions']} | {' | '.join(vals)} | {r['category']} |")
    text.append(f"""
各行的具体限制、AR(10) 比较和抽样影响连到 [证据评估记录]({LINK}evidence_assessment.json)。较稳健的标签表示已测主要方向保持，规格敏感、探索和未决标签对应不同层级的证据。已证实的时钟错误记录保留排除。主要参考分析可以作为带有条件与局限的初步证据公开展示。[验证报告](29_preliminary_validation.md) · [原 Phase 10 跨脑区结果](24_phase10_multiregion_results.md)""")
    (ROOT/"reports/30_preliminary_research_summary.md").write_text("\n".join(text)+"\n",encoding="utf-8")


if __name__=="__main__":
    write_validation()
    write_summary()
