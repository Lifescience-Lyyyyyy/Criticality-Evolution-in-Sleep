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
for mode in pulse parameter parameter-pulse; do
  folder=parameter
  if [ "$mode" = pulse ]; then folder=pulse; fi
  python "$A/intervene.py" "$mode" --upstream "$UPSTREAM" --raw "$CACHE/baseline" \
    --native "$CACHE/native" --matlab-receipt "$A/results/baseline_matlab_receipt.json" \
    --out "$CACHE/interventions-$folder" --workers 5
done
```

参数批次和对应参数脉冲使用同一输出根。汇总的 `parameter_pulses.json` 记录 24 个配对脉冲和未能找到预态窗口的案例。参数批次续跑时，按模拟命令、生成器源码和文件哈希验证已完成的原生输出。中断输出移至同批次的 `interrupted/` 留存。

使用 [参数 MATLAB 批次 45097945176](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/45097945176) 的 ZIP 和回执关联拟合。

```sh
python "$A/join_fits.py" "$A/results/perturb_counts.json.gz" \
  "$CACHE/perturb-matlab.zip" "$A/results/perturb_matlab_receipt.json" "$A/results" --prefix perturb
python "$A/describe_neurons.py" "$A/results/analysis.json" "$A/results/matlab_counts.json.gz" \
  "$A/results/rate_and_events.json"
```

## 发放率敏感性

等率稀释使用原有的 10 单元、60 s 计数，并保持原始时窗。84 条输入在 [MATLAB 批次 45101132727](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/45101132727) 以 AR(5) 和 AR(10) 拟合。群体平均率和三次技术重复的结果见 [等率敏感性摘要](results/rate_sensitivity.json)，求解、残差和完整回执见 [拟合数据](results/rate_fits.json) 和 [回执](results/rate_matlab_receipt.json)。

```sh
python "$A/rate_match.py" prepare "$A/results/matlab_counts.json.gz" \
  "$A/results/baseline_matlab_receipt.json" "$CACHE/rate_counts.json.gz"
python "$A/join_fits.py" "$CACHE/rate_counts.json.gz" \
  "$CACHE/rate-matlab.zip" "$CACHE/rate-receipt.json" "$A/results" --prefix rate
python "$A/rate_match.py" summarize "$A/results/baseline_fits.json" \
  "$A/results/rate_fits.json" "$A/results/rate_sensitivity.json"
```

## Bugnon 原版运行与生理门槛

[模型登记](registry.json)中的复现资料、deNEST 和 NEST 分别独立固定。有效的独立环境为 Python 3.8.20、NumPy 1.23.5、Cython 0.29.37、GSL 2.7、pandas 1.5.3、PyYAML、docopt、tqdm、setuptools 58.0.4 与 libtool。用 NEST 提交的归档树构建，使旧 PyNEST 包名不带 Git 后缀。说明生成器中的 `cgi.escape` 改为 `html.escape`。这项兼容修改只影响文档生成。

```sh
# P 为独立的 ASCII 路径 Python 前缀，NEST_SOURCE 为上述提交的无 Git 归档树。
cmake -S "$NEST_SOURCE" -B "$CACHE/nest-build" \
  -DCMAKE_INSTALL_PREFIX="$P" -DCMAKE_PREFIX_PATH="$P" -DCMAKE_BUILD_TYPE=Release \
  -Dwith-gsl="$P" -Dwith-ltdl=ON -Dwith-openmp=ON -Dwith-mpi=OFF \
  -Dwith-python=ON -DPYTHON_EXECUTABLE="$P/bin/python"
cmake --build "$CACHE/nest-build" -j 4
cmake --install "$CACHE/nest-build"
cmake -S "$BUGNON/extension_module/HTNeuronModule" -B "$CACHE/ht-build" \
  -Dwith-nest="$P/bin/nest-config" -DCMAKE_INSTALL_PREFIX="$P" -DCMAKE_BUILD_TYPE=Release
cmake --build "$CACHE/ht-build" -j 4
cmake --install "$CACHE/ht-build"

source "$P/bin/nest_vars.sh"
export NEST_INSTALL_DIR="$P" NEST_DATA_DIR="$P/share/nest"
export NEST_DOC_DIR="$P/share/doc/nest" NEST_MODULE_PATH="$P/lib/nest"
export LD_LIBRARY_PATH="$P/lib/nest:$P/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$P/lib/python3.8/site-packages:$DENEST"
for state in wake sleep; do
  "$P/bin/python" "$A/run_bugnon.py" "$state" --upstream "$BUGNON" \
    --denest "$DENEST" --nest-source "$NEST_GIT_CHECKOUT" --out "$CACHE/bugnon-$state-native"
done
python "$A/analyse_bugnon.py" --wake "$CACHE/bugnon-wake-native" \
  --sleep "$CACHE/bugnon-sleep-native" --authors "$BUGNON" --out "$CACHE/bugnon-analysis"
```

旧 SLI 启动脚本需要 ASCII 路径。如果 NEST 已安装在含中文的目录，通过 ASCII 软链接加载，并将 `NEST_INSTALL_DIR`、`NEST_DATA_DIR`、`NEST_DOC_DIR`、`NEST_MODULE_PATH`、`LD_LIBRARY_PATH` 与 `PYTHONPATH` 指到别名。每态稳定段的原生尖峰按照原单元 ID、类型和层保存。逐细胞 ACF 和时间尺度写入本地运行目录的 `cortical_neuron_dynamics_40ms.npz`，并在生理回执中记录校验值。Wake 与 Sleep-like 首个稳定 60 s 使用同一组记录单元。 本轮两态分别独立运行，开始正式生理检查以前，两份输出都必须写入完成回执。生理汇总核对执行时的参数树与作者归档，以及所有会话的边界和记录标志。作者在 Sleep-like 的 3–5 s 转换段关闭了记录器，该段单列为缺测。正式的 40 ms 计数从 5 s 开始。主输入使用 23,816 个单元，另固定 500 和 10 单元各三组配对采样。另计算至少 10 个尖峰的单元的 ISI CV，并在作者提供的六个皮层群体内对照中位数。只有三项预登记的生理检查全部通过后，分析才写出该批次的 MATLAB 输入。回执归档同时复核 AR(5) 与 AR(10) 的分段长度、总曝光、均值、方差和前 20 个滞后的观测自相关，并核对拟合设置和下载文件的 SHA-256。

### 论文链接的 mGluR 版本

作者后来在 deNEST 的 `MC_sleep` 分支公开了含 mGluR 电导的网络和独立 NEST 扩展，固定提交为 `d8facf6323aa95151398ee4945623147e33491a8`。该示例的细胞类型电导倍率和 mGluR 输入变化与论文 2.17 节一致。示例内自发 Wake 和 Sleep 的长度分别为 3 s 和 6 s。论文 2.18 节报告的稳定曝光均为 120 s。本地将这套源码和上面的归档配置分别登记。[作者示例](https://github.com/TomBugnon/deNEST/tree/d8facf6323aa95151398ee4945623147e33491a8/example_networks/MC_sleep) · [论文方法](https://onlinelibrary.wiley.com/doi/10.1111/ejn.15836)

复用上面的 ASCII 路径 NEST 环境。令 `PAPER` 指向上述提交的 deNEST checkout，`PAPER_BUILD` 指向用于这个扩展的独立构建目录。

```sh
cmake -S "$PAPER/example_networks/MC_sleep/extension_module/HTNeuronModule" \
  -B "$PAPER_BUILD" -Dwith-nest="$P/bin/nest-config" -DCMAKE_BUILD_TYPE=Release
cmake --build "$PAPER_BUILD" --target htneuronmodule_module --parallel 1
PYTHONPATH="$P/lib/python3.8/site-packages:$PAPER" "$P/bin/python" \
  "$A/bugnon_paper_smoke.py" "$PAPER" "$PAPER_BUILD/htneuronmodule.so" "$CACHE/bugnon_paper_smoke.json"
```

兼容测试按 mGluR 接收端口向细胞发送两个尖峰。电导参数翻倍的单元产生两倍的电导响应，无刺激单元的 mGluR 电导保持零。全网络测试为预热、自发输入、切换和 Sleep 各运行 1 ms，用于核对 5 线程网络的构建和切换后的实际参数。[单细胞回执](results/bugnon_paper_smoke.json) · [全网络回执](results/bugnon_paper_network_check.json)

正式运行保留作者的细胞、连接、外部自发输入和状态切换设置，省略示例中的视觉刺激，并将两态的稳定段分别延长至论文所述的 120 s。Sleep 在 2 s 预热和 3 s 自发输入之后执行 3 s 转换，随后进入稳定段。Wake 和 Sleep 分别从相同的作者种子初始化。

```sh
for state in wake sleep; do
  PYTHONPATH="$P/lib/python3.8/site-packages:$PAPER" "$P/bin/python" \
    "$A/run_bugnon_paper.py" "$state" "$PAPER" "$PAPER_BUILD/htneuronmodule.so" \
    "$CACHE/bugnon-paper-$state-120s" || exit 1
done
```

两态写出完成回执后，可用同一工具验收并准备现有 MATLAB 工作流的输入。

```sh
python "$A/analyse_bugnon.py" --paper --wake "$CACHE/bugnon-paper-wake-120s" \
  --sleep "$CACHE/bugnon-paper-sleep-120s" --out "$CACHE/bugnon-paper-analysis"
```

## 交付

| 文件 | 内容 |
|---|---|
| [复现回执](results/physiology.json) · [作者子集校验](results/author_smoke_proof.json) | 三个配对实例每态的生理验收、敏感性、代码与输出摘要 |
| [动力学统计](results/analysis.json) | 按原 ID 对齐的单元统计、群体 ACF、协方差恒等式、对照分解和等时实测 |
| [MATLAB 归档清单](results/baseline_matlab_receipt.json) · [轻量拟合](results/baseline_fits.json) · [完整诊断](results/baseline_full_matlab.json.gz) | 数值通过、拟合残差、配对 d₂ 与校验 |
| [图 1 原生尖峰和 UP/DOWN](figures/01_original_rasters.pdf) · [图 2 时间尺度](figures/02_timescales.pdf) | 两种原作者条件的 E/I 活动和逐单元／群体比较 |
| [图 3 协方差分解](figures/03_covariance_decomposition.pdf) · [图 4 窗内循环移位](figures/04_circular_shift.pdf) | 自项、跨单元项、群体 ACF 与配对 d₂ |
| [片段内动力学](results/episode_conditioned.json) · [图 6 配对脉冲](figures/06_paired_pulse_response.pdf) | 连续状态段的 ACF 与协方差，12 个脉冲与对照的传播、峰值和恢复 |
| [参数扰动](results/perturbations.json) · [参数拟合](results/perturb_fits.json) · [单元率与事件](results/rate_and_events.json) | 12 个配对参数实例的状态、估计器与事件细节 |
| [参数脉冲](results/parameter_pulses.json) · [图 8](figures/08_parameter_sensitivity.pdf) | 24 个配对脉冲及其恢复、传播和每个网络的参数敏感性 |
| [实测采样敏感性](results/empirical_sensitivity.json) · [图 7](figures/07_empirical_sampling.pdf) | 4 只大鼠在固定单元数、全部单元和不同连续时长的同动物对照 |

图中阴影范围在有多个网络时是跨实例最小值至最大值。Figure 7 条件示例 raster 为第一个预先固定的实例和第 12–18 秒。相同实例各条件使用配对的拓扑和噪声设计。多个窗及代理种子嵌套于同一网络，群体图表以网络实例为单位。实测以动物为单位。
