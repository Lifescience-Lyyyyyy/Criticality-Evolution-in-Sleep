# 睡眠中的神经临界性演化

用 Sooter 等人的 MATLAB 工具包复现小鼠 V1 的 Wake–NREM 比较。斑胸草雀 LMAN 重分析按亮灯活动、亮灯非活动和评分 SWS 组织。严格单元配对缺少合格群体，当前的鸟数据结果用于诊断。历史两态分析已由 [v2](reports/04_zebra_analysis.md) 替换，可从 [历史标签](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/tree/pre-zebra-reanalysis-20261120) 复算。

## 执行

- 在 [MATLAB PNAS Replication](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/workflows/matlab.yml) 手动选择 `smoke`、`synthetic`、`mouse` 或 `zebra`。对应的 `matlab-smoke-*`、`matlab-synthetic-*`、`matlab-mouse-*`、`matlab-zebra-*` 标签也能触发。运行完成后下载 `matlab-results`。
- [Drive 初始化](notebooks/00_drive_setup.ipynb) 和 [小鼠数据准备](notebooks/01_mouse_data.ipynb) 在 Colab 执行。原始数据和运行归档位于 `MyDrive/Research/Criticality_Evolution/`。
- 复用 Drive 中的斑胸草雀源文件。若缺少原始文件，在 WSL 运行 `python3 scripts/data/stage_zebra_sources.py --session <已挂载 Drive 的 Colab 会话名>`。缺少四个补充文件时加 `--supplement-only`。在 Colab 运行[数据准备 notebook](notebooks/02_zebra_data.ipynb)，生成 `matlab_inputs/zebra_lman_reanalysis_v2.mat` 和 `manifests/zebra-finch-preparation-v2.json`，以及 `processed/zebra_finch/lman-burst-events-v2.json`。校验 SHA 后更新 Git 计数文件及同名 `.sha256`，再运行 MATLAB `zebra` 阶段。
- 下载新运行产物中的 `zebra_finch_reanalysis.json`，在 Colab 环境从仓库根目录合并并绘图。历史复算清单也可从正式发布 JSON 的 `historical_reproduction` 字段恢复。

```bash
python scripts/analysis/summarize_zebra.py path/to/zebra_finch_reanalysis.json \
  manifests/zebra-finch-preparation.json path/to/zebra-baseline-repeat.json \
  results/zebra_finch/lman.json
python scripts/analysis/plot_zebra.py
```

## 结果

- [中文总报告](reports/05_final_summary.md) · [小鼠 Wake–NREM 图](figures/pnas/mouse_wake_nrem.pdf)
- 斑胸草雀诊断图　[三态 d₂](figures/zebra_finch/lman_state_comparison.pdf) · [发放统计](figures/zebra_finch/lman_firing_statistics.pdf) · [爆发统计](figures/zebra_finch/lman_burst_statistics.pdf) · [SWS 爆发与 d₂](figures/zebra_finch/lman_burst_d2.pdf) · [敏感性](figures/zebra_finch/lman_sensitivity.pdf)
- [来源清单](manifests/sources.json) · [斑胸草雀准备清单](manifests/zebra-finch-preparation.json) · [实验记录和校验值](manifests/experiments.json)

[Sooter 等人的 PNAS 论文](https://doi.org/10.1073/pnas.2526984123) · [原始 MATLAB 工具包](https://github.com/sam-sooter/prox_crit_toolkit) · [Lorenz 等人的斑胸草雀录制](https://doi.org/10.12751/g-node.vsxyi8)
