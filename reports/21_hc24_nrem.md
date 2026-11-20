# hc-24 作者标注 NREM 分析

在第一场同一只大鼠的同一组 NREM 时间箱内，mPFC 的 AR(5) d₂ 比 CA1 高 0.248 bits/s。AR(10) 的差为 0.238，每段两端各裁掉 5 秒后为 0.256。第三场为跑后休息，只有 CA1 达到单元门槛。两个会话来自不同动物，数据只支持描述性估计。

| 会话 | 脑区 | 排序单元 | 参考 NREM 秒数 | AR(5) d₂ | AR(10) d₂ | 裁边 AR(5) d₂ | 群体放电率（Hz） | AR(5) 残差超额峰度 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | CA1 | 5 | 2530.88 | 2.000 | 0.926 | 1.990 | 5.79 | 17.35 |
| 1 | mPFC | 6 | 2530.88 | 2.247 | 1.164 | 2.246 | 6.85 | 6.74 |
| 3 | CA1 | 7 | 4676.80 | 2.742 | 1.492 | 2.738 | 1.32 | 87.47 |

d₂ 单位为 bits/s。用 40 ms 群体计数、Yule–Walker 拟合和 β=2，在固定的 Sooter 工具包和 MATLAB R2024b 中运行。作者标注的 NREM 分段分别传入估计器。Colab CPU 独立重算的 6 组计数均与 MATLAB 输入一致。两个会话的参考拟合都在启用和关闭求解器诊断时逐值复现，绝对差为 0。6 个 AR(5) 和 3 个 AR(10) 拟合均通过求解器、稳定性及时间箱校验。[第一场 Actions](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44712040641) · [第三场 Actions](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44712586860) · [逐拟合诊断](../results/sleep_crcns/hc-24/)

第一场两区使用同一批时间箱，排序单元数和群体放电率仍有差异。d₂ 幅度随 AR 时间跨度改变。脑区差值仅来自一只大鼠，不能外推区域排序。所有丘脑记录的已发布排序单元数低于预设门槛 5。作者评分文件仅提供 NREM，故无法估计 Wake 或 REM 的对比。[数据与时间核查](17_phase9_data_feasibility.md) · [CRCNS 数据集](https://doi.org/10.6080/K0K35RVG) · [原论文](https://doi.org/10.7554/eLife.48881)
