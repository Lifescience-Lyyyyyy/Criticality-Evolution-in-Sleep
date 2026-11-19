# 睡眠相关神经临界性复现与跨物种探索

2026 年 11 月 20 日

**小鼠的主方向得到复现。斑胸草雀的两个配对录制给出相反方向。** 所有 d₂ 数值均由 GitHub Actions 上的 MATLAB R2024b 和固定版本的 Sooter 工具箱计算。数据准备在 Colab 执行，大体量源文件及运行归档保存在项目 Drive。[工作流](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/workflows/matlab.yml)

| 验收 | 状态 | 已核对的结果 |
| --- | --- | --- |
| 第 0 阶段 云端可行性 | PASS | MATLAB、原始示例、Drive 持久化及跨环境校验通过 |
| 第 1 阶段 合成验证 | PASS | 官方 AR(3) 示例复算一致，SQP 收敛。AR(1) 的估计与解析值符合数值积分精度 |
| 第 2 阶段 小鼠方向复现 | PASS | 全部 13 个会话的 Wake d₂ 低于 NREM。配对均差 1.601329 bits/s，方向与 Fig. 4D 一致 |
| 第 3 阶段 斑胸草雀资格 | PASS | ZF1 和 ZF2 有可信状态评分、配套已排序脉冲及配对 LMAN 探针映射 |
| 第 4 阶段 斑胸草雀计算 | PASS | ZF1 的 SWS d₂ 较高，ZF2 较低。组均差随部分敏感性配置改变符号 |
| 第 5 阶段 跨物种比较 | PASS | 可以比较各自物种内的方向。现有两个鸟样本无法确定物种层面的效应 |

## 解释

小鼠复现使用论文指定的全部 13 个会话、以 40 ms 分箱，并按论文方法对状态片段执行 Yule–Walker 拟合。阶数 5 来自作者鼠数据示例。YMV10 的零时长 NREM 区间被排除，保留该条记录的敏感性复算差值约为 0.000002 bits/s。[小鼠复现](02_mouse_replication.md) · [Sooter 等人](https://doi.org/10.1073/pnas.2526984123)

同一估计器可对审计通过的两只斑胸草雀 LMAN 群体计数得到稳定的数值解。亮灯活动与评分 SWS 在 ZF1 和 ZF2 的 d₂ 排序相反。匹配发放率后，ZF1 的差值接近零，ZF2 的差值进一步为负。鸟的个体差、状态相关发放率、运动或鸣唱、不同录制中的单元取样，以及爆发性计数对 AR 模型的影响，都可能改变配对结果。[斑胸草雀分析](04_zebra_analysis.md) · [Lorenz 等人](https://doi.org/10.1523/JNEUROSCI.1903-24.2025)

跨物种比较只使用各物种内部的状态差方向。小鼠 V1 和斑胸草雀 LMAN 属于功能与解剖不同的回路。小鼠状态为 Wake 与 NREM，鸟的对照为亮灯活动与评分 SWS。即使方向一致，这种设计也只能提供提出演化假设的线索。

一个可检验的假设是，在同一只鸟同一群体的连续记录中，按发放率和单元数匹配后，SWS 与清醒的 d₂ 差仍可出现。验证需要更多独立鸟、同一探针跨多次状态转换并追踪单元、有视频或运动传感器支持的状态评分，以及相近采集条件下可比较的脑区。可将逐鸟状态差和 AR 残差诊断同时设为预先定义的终点。

## 可复现材料

- [云端审计](00_feasibility.md)、[合成验证及首次失败](01_synthetic.md)、[小鼠复现](02_mouse_replication.md)、[斑胸草雀资格](03_zebra_dataset.md) 与 [斑胸草雀结果](04_zebra_analysis.md)
- [合成验证图](../figures/pnas/synthetic_d2.pdf)、[小鼠 Wake–NREM 配对图](../figures/pnas/mouse_wake_nrem.pdf) 与 [斑胸草雀 LMAN 配对图](../figures/zebra_finch/lman_wake_sws.pdf)
- [来源和版本](../manifests/sources.json)、[实验索引](../manifests/experiments.json) 与 [运行步骤](../README.md)
