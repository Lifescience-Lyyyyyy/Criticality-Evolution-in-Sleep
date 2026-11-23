# Criticality–Dynamics Empirical Atlas

[中文报告](REPORT_zh.md) · [分析前锁定的假说](hypotheses.md) · [状态级 CSV](results/atlas.csv)

研究从仓库 `667fbe28fa770d5d3754e9bffa437d579f594c21` 开始。各历史群体计数、归档版本和拟合文件固定在 [inputs.py](inputs.py)。表中每条观测的资产 ID、数据摘要、单元 ID、阶段注释、区间及数据层均可追踪。MATLAB 使用 Sooter 工具箱 `f26c782dca3b1c9ab18cc92c6449591c44e0d507`。

## 验证与重现

在 WSL 仓库根目录执行，Python 3.12。MATLAB 计算在 GitHub Actions 运行，Python 预处理使用 NumPy 2.2.6、SciPy 1.15.3、h5py 3.13.0、remfile 0.1.16。作图额外需要 Matplotlib。

```sh
export OPENBLAS_NUM_THREADS=1
A=research/criticality_dynamics_atlas
CACHE=/path/to/atlas-cache
python "$A/tests.py"
python "$A/analyse.py" population --out "$CACHE/pop"
python "$A/native.py" plan --out "$CACHE/plan-all.json"
```

1. 在 Colab 运行 `colab_extract.py`。将生成的计划命名为 `plan.json`。可以按资产拆分，两批分别为 18 个 Phase 8 资产和 11 个 PoSub 资产。只读取 NWB 排序单元的 ID、尖峰时间和索引。每批输出 `unit-windows.zip`，从私有缓存读取该 ZIP。`native.py` 对 ZIP 逐条校验后提取统计量。同一个命令可重复提供 `--archive`。已有的 Zenodo WT 解包目录传给 `--brainstem`，目录应包含 `WT/<subject>/1k_train.npz`。原始神经信号和逐单元矩阵留在缓存。

   ```sh
   python "$A/native.py" analysis --archive "$CACHE/phase8-unit-windows.zip" \
     --archive "$CACHE/posub-unit-windows.zip" --brainstem "$BRAINSTEM_CACHE" \
     --out "$CACHE/native.json"
   ```

2. 在 [Empirical dynamics atlas MATLAB](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/workflows/dynamics-atlas.yml) 运行工作流。[首轮运行](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/45036974767)使用提交 `d878b72b62ad5e1521d1dc7c70e97f060c3fcc99` 用来验证求解。[配对归档运行](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/45037956569)使用提交 `a03caaff7c1a154e88f6485a2c503bb9cad4744d` 同时保留计算 d₂ 时的 Python 统计。两次运行的 1,833 个主拟合和 606 个替代拟合的 d₂ 完全相同。归档包含拟合、轻量统计和输入 SHA-256。用同一工作流提供的 ZIP 和相应 [归档回执](results/matlab_receipt.json) 汇总。

3. 生成动物级结果和图。

   ```sh
   python "$A/summarize.py" "$CACHE" --artifact "$CACHE/actions-45037956569.zip" \
     --receipt "$A/results/matlab_receipt.json" --native-fits "$A/results/native_d2.json" \
     --out "$A/results"
   python "$A/plot.py"
   python scripts/atlas/validate_atlas_data.py
   ```

## 原生分箱及 PAG 三态复算

使用相同的 NWB／PAG 连续窗口生成 20、40、80、200 ms 的 Yule–Walker 输入。每个窗口只导出低阶自协方差、来源校验值和窗口元数据。MATLAB 沿用固定工具箱的目标函数和优化器。群体计数时序和逐单元矩阵保存在缓存。

```sh
python "$A/native_d2.py" prepare "$CACHE" "$BRAINSTEM_CACHE" "$CACHE/native_d2_input.json"
python "$A/native_d2.py" summarize "$A/native_d2_input.json" "$CACHE/actions-45096650326.zip" \
  "$CACHE/receipt-45096650326.json" "$A/results/native_d2.json"
```

[运行 45096650326](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/45096650326) 由 `matlab-dynamics-atlas-native-*` 标签触发。产物 ID、摘要和 114 条历史拟合的重放误差见结果的 `receipt` 字段。复算通过后，`summarize.py` 将 PAG 同窗 40 ms 的 d₂ 合入主表。NWB 的历史 40 ms 值继续用原来的拟合来源。

## 数据文件

| 文件 | 内容 |
|---|---|
| [native_d2.json](results/native_d2.json) | 原生同窗四箱宽拟合、求解器质控、动物配对区间及运行回执 |
| [atlas.csv](results/atlas.csv) | 数据层、动物、会话、区域、状态、曝光、细胞数和主指标。空值表示缺失 |
| [inventory.json.gz](results/inventory.json.gz) | 每项来源、注释、原计数与可进行的分析 |
| [observations.json.gz](results/observations.json.gz) | 会话统计、单位指标、谱／ACF、AR 模态、采样质控和窗口来源 |
| [paired.json](results/paired.json) · [animal_states.json](results/animal_states.json) | 动物级状态值、配对差、bootstrap 区间和正差数量 |
| [associations.json](results/associations.json) | 状态分层、动物内、配对差、率及曝光调整、p 值与 BH q 值 |
| [acf_profile_associations.json](results/acf_profile_associations.json) | 六个滞后 ACF 轮廓的动物级配对主成分及其与 d₂ 的描述性关系 |
| [simultaneous.json](results/simultaneous.json) · [region_interactions.json](results/region_interactions.json) | 同资产、完全相同时间窗的区域差及区域×状态交互 |
| [continuous_windows.json.gz](results/continuous_windows.json.gz) | 每态 60／120／240 s 连续窗口资格和时间尺度敏感性 |
| [surrogates.json.gz](results/surrogates.json.gz) · [surrogate_effects.json](results/surrogate_effects.json) | 三类代理、分箱、活动条件敏感性与重复种子的差值 |
| [null_models.json](results/null_models.json) | AR(1)、混合模态、已知分支参数、共同输入、门控对照 |
| [MATLAB 回执](results/matlab_receipt.json) · [Phase 8 NWB 回执](results/phase8-unit-windows-receipt.json.gz) · [PoSub NWB 回执](results/posub-unit-windows-receipt.json.gz) | 版本、文件校验及运行信息 |

`neurons` 紧随 `meta.unit_ids` 的次序，每条为 [rate_Hz, Fano, 指数τ_s, 积分τ_s]。`pr_raw_n` 使用全部记录单元，`pr_corr_n` 用非零方差单元。`pr_fixed10` 是 12 次固定 N=10 的原始 PR/N 的中位数。`off_duration` 是完整零计数游程的 p95（秒）。`tau_AR_longest` 是拟合稳定根的最长衰减时间（秒）。`fit` 和 `alternate` 分别保留主规格和替代 AR 规格的数值质控。各层之间存在重复观测，不应直接合并做总体显著性检验。

## 图

1. [动物级总览](figures/01_empirical_overview.pdf)
2. [单元和群体时间尺度](figures/02_timescales.pdf)
3. [d₂ 与时间尺度](figures/03_distance_timescale.pdf)
4. [MR 与已知参数对照](figures/04_branching.pdf)
5. [群体有效维度及抽样](figures/05_dimension.pdf)
6. [自／跨协方差与统计代理](figures/06_collective_covariance.pdf)
7. [波动与零计数游程](figures/07_fluctuations.pdf)
8. [ACF、PSD、knee 和尺度敏感性](figures/08_multiscale.pdf)
9. [脑区×状态](figures/09_region_state.pdf)
10. [理论对照](figures/10_theoretical_nulls.pdf)
11. [关联区间及分母](figures/11_association_matrix.pdf)
12. [七个假说的判断与证据](figures/12_hypothesis_evaluation.pdf)

图 1 的右侧数字以 Wake／NREM／REM 的顺序给出有效动物数，并用状态颜色标出。图 2 上两幅以小提琴展示可估计单元，菱形为动物群体 τ，下面的配对差为动物级估计。图 6 的实线为真实观测，虚线为逐细胞移位。图 8F 使用同一批 PoSub 三态窗口，圆点为动物，菱形为中位数。图 9 左两幅的暖／冷色对应配对正／非正效应比例，数字为独立动物数，第三幅仅用同录同窗。图 11 的灰／青／橙色分别对应 A／B／C 类统计依赖。图 10 全部点来自实际运行的有限样本对照。
