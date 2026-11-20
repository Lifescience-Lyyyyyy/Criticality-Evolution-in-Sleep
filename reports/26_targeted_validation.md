# Phase 11 针对性验证

按归档路径、已通过质量门槛的单元数、状态覆盖和已报告的 LDT 例外冻结九场样本，再求补充估计。所有新增拟合均用 Sooter 工具箱 f26c782dca3b1c9ab18cc92c6449591c44e0d507、MATLAB R2024b、40 ms、Yule–Walker 和固定 β=2。共有 41 条输入序列，V1 参考回放另拟合 AR(10)。历史分期、计数、拟合及报告保留。[抽样清单](../results/validation/selected_sessions.json) · [逐态数据与诊断](../results/validation/audit_comparisons.json)

## 实际影响

| 关注点 | 数据与会话 | 原结果 | 核验结果 | 对推断的影响 | 状态 |
|---|---|---|---|---|---|
| 单元选择 | 000056，Mouse28-140313 PoSub | NW +0.03872 | q75 -0.00010、q50 -0.02888 | 单会话的微弱正向变为近零或小幅负向 | CONFIRMED_IMPACT，[曲线](../figures/phase11/sorting_sensitivity.pdf) |
| 脑干排序 | DL108_122320n1 LDT | NW +0.233 | 污染率 q50 NW -0.148 | 唯一一只 LDT 的 NW 对单元选择敏感，RN 仍为正 | CONFIRMED_IMPACT，[曲线](../figures/phase11/sorting_sensitivity.pdf) |
| 阶数敏感 | 000056，Mouse28-140313 PoSub | RW -0.135 | AR(10) RW +0.061 | REM−Wake 变号，两阶数值仍稳定 | CONFIRMED_IMPACT，[逐拟合](../results/validation/audit_comparisons.json) |
| 时钟冲突 | 000939，A3705／200306 | 整场排除 | 宽带和 LFP 同在 11168.136 秒截止，尖峰到 16198.55965 秒 | 25/92 个评分行越界。维持排除 | CONFIRMED_IMPACT，[复核](../results/validation/nwb_source_recheck.json) |
| 整窗样本 | 000939，29 只 | 参考 NW 为 29/29 正向 | 60 秒等时分析保留 15 只，另 14 只均卡在 REM | 子集参考均值偏高，两个组的方向一致 | CONFIRMED_IMPACT，[逐动物窗数](../results/validation/000939_equal_duration.json) |
| 窗长敏感 | 000939，A6216／210623 | 3 × 60 秒 NW +1.501 | 4 × 45 秒 NW +1.043 | 同时长下方向保持，NW 幅度约缩小 30% | CONFIRMED_IMPACT，[逐拟合](../results/validation/matlab-qc.json) |
| 区间一致性 | 000056 的 A、B 和 000939 的 C、E | 原参考区间 | 落在独立采集时钟、作者评分和核验过的行为边界内 | 原区间继续使用 | NO_OBSERVED_IMPACT，[逐区间](../results/validation/audit_comparisons.json) |
| 原估计回放 | 小鼠 V1，YMV02 | NW +1.462 | AR(5) 逐态精确重放。AR(10) NW +0.939 | 方向保持，幅度对阶数敏感 | NO_OBSERVED_IMPACT（回放及方向），[实测](../results/validation/matlab-v1.json) |
| 群体长静默 | hc-23，21 场 | 按尖峰支持切段 | 8 场支持段数超过保留的作者 bout 数 | 预处理切段生效。采集连续性与 d₂ 变化仍缺独立证据 | INCONCLUSIVE，[切段证据](../results/validation/audit_comparisons.json) |
| REM 缺失 | 小鼠 V1，公开 13 场 | 仅 Wake 和 NREM | 本轮回放仍为两态 | 无法评估该公开复现的 REM 差值 | NOT_APPLICABLE，[既有结果](02_mouse_replication.md) |

NW 为 NREM−Wake，RN 为 REM−NREM，RW 为 REM−Wake。估计和差值单位为 bits/s。状态判断对应表中具体问题。抽样结论的适用范围为所列会话和所试规格。完整状态索引保存在 [影响表](../results/validation/impact_assessment.json)。

## A · 排序与解剖

000056 的 pE、pI 是由波形和放电特征推定的细胞类型。L-ratio 和 Mahalanobis isolation distance 表征聚类分离。作者 isi_index 为不应期违反率。in_lab_cell_list 标记进入实验室的分析名单。源文件说明原论文的逐单元名单无法完整恢复。000939 的三场没有这些隔离指标，原严格 ISI 使用从排序时钟重新计算的不应期违反率。000552 的 putative_cell_type 是 CellExplorer 分类。脑干检查作者 good、mua、noise 标签及通道核团，检查原单位 ID 与 1 ms 尖峰哈希。ContamPct 是排序表中的污染估计。[原 NWB 字段与时钟](../results/validation/nwb_source_recheck.json) · [脑干原排序表收据](../data/fixtures/phase11/shorter/brainstem_source_recheck.json)

公开资料缺少可推广的质量阈值。探索对照在原群体内按有限 L-ratio 或 ContamPct 的第 75 或 50 百分位取较低值，跨状态保持同一组单元、原时间箱和原分期。图和数据列出每场单元数量及跨阈值变化。原严格 ISI 和新增百分位作为各自独立对照。[排序曲线](../figures/phase11/sorting_sensitivity.pdf)

| 会话／脑区 | 原区域单元数 | 参考保留 | 归档排除理由与数量 |
|---|---:|---:|---|
| Mouse12-120809 ADN | 50 | 50 | 原审计未列区域内排除 |
| Mouse28-140313 ADN | 22 | 22 | 原审计未列区域内排除 |
| Mouse28-140313 PoSub | 38 | 38 | 原审计未列区域内排除 |
| 210623 PoSub | 75 | 74 | 时间戳重复或逆序 1 |
| 200306 PoSub | 117 | 0 | 尖峰越出采集时钟 117 |
| 190910 PoSub | 74 | 71 | 时间戳重复或逆序 3 |
| DL159_7721n1 PAG | 234 | 115 | 作者排序 noise 44、作者排序 mua 75 |
| DL108_122320n1 LDT | 36 | 18 | 作者排序 mua 12、作者排序 noise 6 |
| e13_1m1_200228 CA1 | 83 | 83 | 原审计未列区域内排除 |
| YMV02 V1 | 85 | 85 | 原审计未列区域内排除 |

### 单元对照逐态结果

| 会话／脑区 | 群体 | 单元 | Wake | NREM | REM | NW／RN |
|---|---|---:|---:|---:|---:|---:|
| Mouse12-120809 ADN | 参考 AR(5) | 50 | 0.933 | 2.291 | 0.766 | +1.358／-1.526 |
| Mouse12-120809 ADN | 原严格 ISI | 47 | 0.798 | 2.191 | 0.636 | +1.393／-1.555 |
| Mouse12-120809 ADN | L-ratio q75 | 37 | 0.838 | 2.265 | 0.650 | +1.427／-1.615 |
| Mouse12-120809 ADN | L-ratio q50 | 25 | 0.773 | 2.103 | 0.597 | +1.330／-1.506 |
| Mouse28-140313 ADN | 参考 AR(5) | 22 | 0.623 | 1.420 | 0.372 | +0.797／-1.048 |
| Mouse28-140313 ADN | L-ratio q75 | 16 | 0.608 | 1.451 | 0.426 | +0.843／-1.025 |
| Mouse28-140313 ADN | L-ratio q50 | 11 | 0.546 | 1.421 | 0.380 | +0.875／-1.041 |
| Mouse28-140313 PoSub | 参考 AR(5) | 38 | 1.165 | 1.204 | 1.031 | +0.039／-0.174 |
| Mouse28-140313 PoSub | L-ratio q75 | 28 | 1.228 | 1.228 | 1.061 | -0.000／-0.166 |
| Mouse28-140313 PoSub | L-ratio q50 | 19 | 1.276 | 1.247 | 1.084 | -0.029／-0.163 |
| 210623 PoSub | 参考 AR(5) | 74 | 0.507 | 1.548 | 0.971 | +1.041／-0.577 |
| 210623 PoSub | 原严格 ISI | 68 | 0.506 | 1.558 | 0.938 | +1.053／-0.620 |
| 190910 PoSub | 参考 AR(5) | 71 | 0.583 | 1.742 | 0.724 | +1.159／-1.018 |
| 190910 PoSub | 原严格 ISI | 68 | 0.586 | 1.737 | 0.732 | +1.151／-1.005 |
| DL159_7721n1 PAG | 参考 AR(5) | 115 | 0.975 | 1.312 | 0.837 | +0.337／-0.475 |
| DL159_7721n1 PAG | 原严格 ISI | 114 | 0.984 | 1.322 | 0.845 | +0.337／-0.476 |
| DL159_7721n1 PAG | 污染率 q75 | 87 | 1.021 | 1.327 | 0.879 | +0.306／-0.448 |
| DL159_7721n1 PAG | 污染率 q50 | 62 | 1.073 | 1.462 | 1.016 | +0.390／-0.446 |
| DL108_122320n1 LDT | 参考 AR(5) | 18 | 1.167 | 1.399 | 1.532 | +0.233／+0.133 |
| DL108_122320n1 LDT | 污染率 q75 | 13 | 1.371 | 1.460 | 2.034 | +0.089／+0.574 |
| DL108_122320n1 LDT | 污染率 q50 | 9 | 1.702 | 1.554 | 2.183 | -0.148／+0.629 |
| e13_1m1_200228 CA1 | 参考 AR(5) | 83 | 1.591 | 2.210 | 1.264 | +0.619／-0.946 |
| YMV02 V1 | 参考 AR(5) | 85 | 0.561 | 2.023 | — | +1.462／— |

## B · 分期、行为和独立时钟

逐场重新读出公开 NWB 的时钟长度、评分与行为表。000056 先区分原论文评分和后续 Buzcode，同名状态重叠及参考有效暴露在逐态数据中给出。两者为同一批记录的不同评分。A、B 所用区间在休息和原 sub-recording 支持内。B 同步 ADN／PoSub 的每态时间箱一致。C、E 在独立 LFP、作者状态标注和家笼支持内，避开光刺激及其前后 60 秒。D 的原时钟不相容并已排除。全程保留原始录制时间单位。

## C · 等时选择

[逐动物审计](../results/validation/000939_equal_duration.json) 保存 29 只的三态时长、候选窗数、参考单元数和效应。原规则要求所有可用参考状态各至少有三个连续 60 秒窗。14 只未纳入均由于已有 REM 短或碎，Wake 与 NREM 窗均达门槛。

| 60 秒整窗组 | 动物 | 参考 NW | 中位 Wake 秒 | 中位 NREM 秒 | 中位 REM 秒 | 中位单元数 |
|---|---:|---:|---:|---:|---:|---:|
| 纳入 | 15 | +1.433 | 3009.12 | 4369.00 | 380.00 | 75 |
| 排除 | 14 | +1.197 | 2727.82 | 4311.00 | 246.88 | 69.5 |

两组参考效应各自全为正。原纳入组的等时平均 NW 为 +1.440。固定 180 秒和同组单元，先用典型 C 比较 3 × 60 秒和 4 × 45 秒。两规格均收敛并保留 NW 和 RN 方向，RN 为 -0.927 与 -0.802。通过这场方向门槛后才用较短窗核验排除个体 E。A1821／190910 用 4 × 45 秒得到 NW／RN +0.408／-0.635 bits/s。 固定 180 秒的 45 秒候选窗规则理论上有 21 只达标。较短窗的幅度敏感，保持原 60 秒规格和独立样本口径。

## D · AR 诊断

| 会话／脑区 | AR(5) NW／RN | AR(10) NW／RN | AR(5) 最大绝对残差 ACF | AR(10) 最大绝对残差 ACF | 两阶最小 AR 根 |
|---|---:|---:|---:|---:|---:|
| Mouse12-120809 ADN | +1.358／-1.526 | +0.783／-0.852 | 0.062 | 0.026 | 1.0504 |
| Mouse28-140313 ADN | +0.797／-1.048 | +0.497／-0.523 | 0.069 | 0.039 | 1.0366 |
| Mouse28-140313 PoSub | +0.039／-0.174 | +0.119／-0.058 | 0.058 | 0.042 | 1.0679 |
| 210623 PoSub | +1.041／-0.577 | +0.485／-0.251 | 0.094 | 0.064 | 1.0295 |
| 190910 PoSub | +1.159／-1.018 | +0.585／-0.473 | 0.085 | 0.044 | 1.0370 |
| DL159_7721n1 PAG | +0.337／-0.475 | +0.161／-0.196 | 0.043 | 0.020 | 1.0735 |
| DL108_122320n1 LDT | +0.233／+0.133 | +0.045／+0.108 | 0.056 | 0.035 | 1.0612 |
| e13_1m1_200228 CA1 | +0.619／-0.946 | +0.417／-0.595 | 0.042 | 0.024 | 1.0884 |
| YMV02 V1 | +1.462／— | +0.939／— | 0.071 | 0.046 | 1.0358 |

补充估计的求解正常且 AR 根在单位圆外。残差峰度、残差 ACF、观测和模型 ACF 均保存在 [逐态诊断](../results/validation/audit_comparisons.json)。YMV02 Wake 的 ACF RMSE 从约 0.137 降至 0.047。残差峰度与相关性阈值描述模型偏离，数值收敛和 AR 稳定性是另一组判断。LDT 的 RN 在两个阶数下均为正，单只证据保留为探索性。β=2 固定参考模型流形，AR 系数逐会话和状态估计，生理相变仍需独立实证。

## E · 连续性与长静默

Phase 8 逐会话有完整静默记录的 146 场，以及 Phase 10 的 74 个主要合格记录资产，原群体审计均无大于 10 秒的峰间隙触发。hc-23 的 8/21 场出现额外的支持切段。例如 Rat-26 Session2 从 54 个保留作者 bout 形成 65 个支持段。尖峰时钟经原 Neuralynx .t 核验，但这些切点缺少能够证明采集是否实际中断的独立连续记录证据。恢复切点后 d₂ 会改变多少仍不确定。hc-24 的跑后 CA1 最大连续零计数 **55.28 秒**保留在拟合中。继续用独立时钟界定采集有效范围，并保留原预处理。

## 输入与执行

[抽样脚本](../scripts/atlas/select_validation.py) · [公开源复核](../scripts/atlas/recheck_sources.py) · [固定时长重选窗](../scripts/atlas/build_qc_fixture.py) · [Colab 与本机计数逐数组核验](../results/validation/colab_preparation_receipt.json) · [完整影响表](../results/validation/impact_assessment.json)

- [44775203800](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44775203800) · 2 条输入 · 提交 7cd4091c47a45886b954c28b1c6f87001ffec00c · 输入 SHA-256 9a8db562f74f61caf965bde24f1ee709dbdc31cc290c9507f36e9e425b9881c9。
- [44776271221](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44776271221) · 24 条输入 · 提交 a801459f64c5278686023b3de44fde886d204f48 · 输入 SHA-256 f27a6f735b786b4b7aaab09da33ae6e90870a3d0a82031e0ed661066f24a2cf9。
- [44777991140](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44777991140) · 15 条输入 · 提交 4046763a2c7e7075ee482540faadf45d9dde91aa · 输入 SHA-256 d5749421e0cba399c452c115d0d896bd4d720c1c6bafbd793e4e934ee8c1e53f。
