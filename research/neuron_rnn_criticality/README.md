# Neuron-level Recurrent Network Criticality Benchmark

[中文研究报告](REPORT_zh.md) · [模型与参数登记](registry.json) · [条件、门槛与反事实协议](protocol.md) · [机器可读结果](results/)

研究基线仓库为 `5d7d03f1edd4ad0573ef970603103a8274425a1a`。Jercog 原生脉冲源码为 `11603e5a32e71664aee8b0d2945ef7ec40d61b37`。距临界分析复用冻结的 Sooter MATLAB 工具箱 `f26c782dca3b1c9ab18cc92c6449591c44e0d507` 及既有 Gaussian AR 诊断器。Bugnon 和 Krishnan 的状态见报告。旧的 [均场模型基准](../sleep_model_benchmark/) 和 [实证图谱](../criticality_dynamics_atlas/) 可用于独立对照。

## Jercog 复现

以下命令从 WSL 的 Git 仓库根目录运行。Python 分析使用 Python 3.12、NumPy、SciPy、h5py 和 remfile，作图额外使用 Matplotlib。C++ 用 Ubuntu g++ 13.3.0 构建。

```sh
export OPENBLAS_NUM_THREADS=1
A=research/neuron_rnn_criticality
CACHE=/path/to/private-model-cache
UPSTREAM="$CACHE/jercog"
git clone https://github.com/djercog/JercogEtAl2017.git "$UPSTREAM"
git -C "$UPSTREAM" checkout --detach 11603e5a32e71664aee8b0d2945ef7ec40d61b37

python "$A/tests.py"
python "$A/run_jercog.py" --upstream "$UPSTREAM" --out "$CACHE/smoke" --seconds 5 --workers 2 --smoke
python "$A/run_jercog.py" --upstream "$UPSTREAM" --out "$CACHE/baseline" --seconds 100 --workers 6
python "$A/physiology.py" "$CACHE/baseline" "$CACHE/native"
python "$A/analyse_snn.py" analysis "$CACHE/native" "$CACHE/analysis"
python "$A/empirical_sensitivity.py" "$CACHE/empirical_counts.json.gz"
```

[工作流](../../.github/workflows/neuron-rnn-criticality.yml) 支持手动选择批次，也可由含该批次的 `matlab-neuron-rnn-*` 标签触发。实际的作者条件与主要采样对照来自 [运行 45042234880](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/45042234880)。Git 中的轻量输入是逐案例群体计数和源哈希。用各自 ZIP 及匹配的回执归档得到拟合和异常记录。

```sh
python "$A/join_fits.py" "$A/results/matlab_counts.json.gz" \
  "$CACHE/baseline-matlab.zip" "$A/results/baseline_matlab_receipt.json" "$A/results"
python "$A/plot_snn.py" "$A/results" "$A/figures"
```

原始 Spike 文件由模型导出并私有缓存。运行目录的 `simulation_record.npz` 保存排序后的 `neuron_ids`、`neuron_types`、原生 `spike_times_seconds` 和 ragged `spike_times_offsets`。第 i 个单元的尖峰坐标位于 offsets[i] 至 offsets[i+1] 的半开区间。相邻的 `manifest.json` 包含模型、原提交、参数、条件、种子、初始化、实际记录数、仿真时步、原文件和派生文件校验。模型内部状态的附加统计可用 python "$A/analyse_snn.py" epochs "$CACHE/native" "$CACHE/episode_conditioned.json" 重建。另一个缓存 `neuron_dynamics_40ms.npz` 保留全部逐单元计数、逐单元 ACF、发放率、指数 τ 和积分 τ。复现时选用与原运行不同的输出目录。

## 反事实与参数敏感性

完成作者条件的生理验收和 MATLAB 回执之后，使用相同拓扑与噪声种子运行配对脉冲和预先登记的四种参数扰动。

```sh
for mode in pulse parameter; do
  python "$A/intervene.py" "$mode" --upstream "$UPSTREAM" --raw "$CACHE/baseline" \
    --native "$CACHE/native" --matlab-receipt "$A/results/baseline_matlab_receipt.json" \
    --out "$CACHE/interventions-$mode" --workers 5
done
```

参数批次续跑时，按模拟命令、生成器源码和文件哈希验证已完成的原生输出。中断输出移至同批次的 `interrupted/` 留存。

## 交付

| 文件 | 内容 |
|---|---|
| [复现回执](results/physiology.json) · [作者子集校验](results/author_smoke_proof.json) | 三个配对实例每态的生理验收、敏感性、代码与输出摘要 |
| [动力学统计](results/analysis.json) | 按原 ID 对齐的单元统计、群体 ACF、协方差恒等式、对照分解和等时实测 |
| [MATLAB 归档清单](results/baseline_matlab_receipt.json) · [轻量拟合](results/baseline_fits.json) · [完整诊断](results/baseline_full_matlab.json.gz) | 数值通过、拟合残差、配对 d₂ 与校验 |
| [图 1 原生尖峰和 UP/DOWN](figures/01_original_rasters.pdf) · [图 2 时间尺度](figures/02_timescales.pdf) | 两种原作者条件的 E/I 活动和逐单元／群体比较 |
| [图 3 协方差分解](figures/03_covariance_decomposition.pdf) · [图 4 窗内循环移位](figures/04_circular_shift.pdf) | 自项、跨单元项、群体 ACF 与配对 d₂ |
| [片段内动力学](results/episode_conditioned.json) · [图 6 配对脉冲](figures/06_paired_pulse_response.pdf) | 连续状态段的 ACF 与协方差，12 个脉冲与对照的传播、峰值和恢复 |
| [实测采样敏感性](results/empirical_sensitivity.json) · [图 7](figures/07_empirical_sampling.pdf) | 4 只大鼠在固定单元数、全部单元和不同连续时长的同动物对照 |

图中阴影范围在有多个网络时是跨实例最小值至最大值。Figure 7 条件示例 raster 为第一个预先固定的实例和第 12–18 秒。相同实例各条件使用配对的拓扑和噪声设计。多个窗及代理种子嵌套于同一网络，群体图表以网络实例为单位。实测以动物为单位。
