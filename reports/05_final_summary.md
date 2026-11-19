# 小鼠 PNAS 复现与斑胸草雀 LMAN 重分析

2026 年 11 月 20 日

**13 个小鼠会话均复现了 Wake d₂ 低于 NREM 的方向。斑胸草雀的严格单元配对受阻，现有两只鸟的数据不足以判断 LMAN 的临界性是否随行为状态改变。** 三态结果、重复控制与爆发分析用于诊断。历史两态分析标记为 **SUPERSEDED / NOT VALIDATED**，并在 Git 标签中可复算。[历史标签](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/tree/pre-zebra-reanalysis-20261120)

| 验收 | 状态 | 证据 |
|---|---|---|
| 云端可行性及合成验证 | PASS | Drive 到 GitHub Actions 校验通过，原版 MATLAB AR(3) 示例和 AR(1) 解析值验证通过 |
| 小鼠方向复现 | PASS | 13 个会话的 NREM − Wake 配对均差为 1.601329 bits/s，重现 Fig. 4D 方向 |
| A 状态 | PASS | 四次会话的评分、AP 时钟和探针身份核对一致 |
| B 严格群体 | BLOCKED | 只有 ZF1 的 SWS 达到两个严格单元，无法形成任一严格三态配对 |
| C 历史数值 | PASS | 4 个主要结果与 24 个诊断 d₂ 精确复算 |
| D 爆发检测 | PASS | 密集栅格计算复核并用作者归档时点校验，源码第 97.5 百分位与论文第 99 百分位均保留 |
| E 最终报告 | PASS，诊断范围 | 三态宽口径估计、单元及发放率各 100 次重复、模型和整窗诊断、严格配对限制均与发布结果对照 |

小鼠 V1 的分析使用论文指定的全部 13 个会话、40 ms 分箱和 Yule–Walker 拟合。阶数 5 来自作者小鼠示例。YMV10 的零时长 NREM 区间排除后，保留该记录的敏感性复算差值约为 0.000002 bits/s。[小鼠复现](02_mouse_replication.md) · [Sooter 等人](https://doi.org/10.1073/pnas.2526984123)

斑胸草雀仍显示相反的 SWS − 亮灯活动方向。用同一亮灯录制的非活动区间作参照时，两只鸟的方向也相反，而且各自翻转。单元筛选资格受限、状态选择和伪迹掩码影响配对估计，发放率稀释同时改变方差和时间结构。SWS 60 秒窗的爆发频率与 d₂ 未显示一致关联，发放率和时间匹配后 ZF1 仅四对，ZF2 仅两对。SWS 残差偏离高斯，在 ZF2 尤其明显。亮灯状态来自麦克风，跨录制无法追踪同一单元。[斑胸草雀结果与八项问题](04_zebra_analysis.md)

下一步优先在更多独立的鸟上连续记录同一 LMAN 群体，覆盖多次活动、安静清醒和 SWS 转换。同步视频、麦克风、运动与 LFP，核验单元跨状态稳定性。预注册状态对照、相近录制时刻和发放率的匹配范围、AR 诊断及逐鸟效应，以检验状态差是否可在这些条件下复现。

## 结果与复现

- [云端审计](00_feasibility.md) · [合成验证](01_synthetic.md) · [小鼠复现](02_mouse_replication.md) · [小鼠配对图](../figures/pnas/mouse_wake_nrem.pdf)
- [斑胸草雀数据资格](03_zebra_dataset.md) · [重分析审计](06_zebra_reanalysis_audit.md) · [变更记录](06_zebra_reanalysis_changelog.md) · [三态诊断图](../figures/zebra_finch/lman_state_comparison.pdf)
- [结果 JSON](../results/zebra_finch/lman.json) · [实验索引](../manifests/experiments.json) · [复现步骤](../README.md)
