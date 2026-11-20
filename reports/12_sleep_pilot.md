# Phase 8B · 自发睡眠 ACC 试跑

**S1 PASS。** DANDI `000041@0.261104.0300` 的 BWRat17 ACC 会话 `BWRat17_121912` 通过原始 MATLAB 估计与优化诊断回放。三个参考结果逐值一致，最大绝对差为 0，优化退出为正，AR 根模大于 1。35 个合格分群在三种状态内相同。数据来自无指定任务或刺激的家笼录制。[DANDI 资产](https://api.dandiarchive.org/api/dandisets/000041/versions/0.261104.0300/assets/d89041a6-42b7-4c7f-8da0-d2c65061e04a/) · [MATLAB 运行](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44679472413)

| 作者状态 | 连续片段 | 秒 | d₂（bits/s） | 群体脉冲/s | Fano | 零计数箱 | 残差超额峰度 | ACF RMSE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Awake | 35 | 2573 | 1.946158 | 30.57 | 1.195 | 32.65% | 1.666 | 0.060 |
| Non-REM | 31 | 3665 | 2.030974 | 30.52 | 1.639 | 40.68% | 2.325 | 0.018 |
| REM | 5 | 216 | 1.872957 | 38.45 | 1.164 | 23.43% | 1.419 | 0.053 |

40 ms、AR(5)、β=2、Yule–Walker 的 NREM − Wake 为 +0.084816 bits/s。等时长 1680 秒给出 +0.145743，片段端点各剔除 5 秒时为 +0.130546。AR(10) 的 NREM − Wake 为 +0.279601。参考规格下 REM − NREM 为 −0.158017 bits/s。Wake 与 NREM 的平均计数率几乎一样，但 Fano、静默和残差分布仍有区别。

最大绝对残差 ACF（1–20 滞后）分别为 0.037、0.013、0.040。方差、20 滞后观测与模型 ACF、创新方差、AR 系数与静默段统计逐项归档。参考拟合的残差均重尾，三态 d₂ 仅用于该模型下的比较。此运行显示一个会话的方向与既有小鼠 V1 结果相同，独立动物层面的结论由后续队列给出。

试跑输入为 161,853 字节的分段计数 fixture。Colab、Drive 和 GitHub Actions 的 SHA-256 均为 `c78d79a6434b6ce65535566ac5876ef5bac8e9efc24b5457d9c4917fe903f29a`。发布使用固定 Sooter 提交 `f26c782dca3b1c9ab18cc92c6449591c44e0d507` 与 MATLAB R2024b Update 10。官方 smoke 先于真实试跑通过。

[试跑 d₂ 和残差图](../figures/sleep_dandi/000041/pilot_acc.pdf) · [MATLAB 诊断 JSON](../results/sleep_dandi/000041/pilot.json)
