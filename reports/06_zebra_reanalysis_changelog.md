# 斑胸草雀重分析变更记录

2026 年 11 月 20 日

历史主分支 `ae1ec2cf473cff6efa3a5296c74a766bba0eb37c` 已以 `pre-zebra-reanalysis-20261120` 保留。旧版两态结论标记为 **SUPERSEDED / NOT VALIDATED**。已验证的 MATLAB 分析代码和输入固定在 `514b6e54d12afb3d89be41805dc0857deec186e1`（标签 `matlab-zebra-reanalysis-seeded-20261120`）。发布也保留历史复算与开发验证的工作流记录。鼠、合成、上游工具包及 GitHub Actions 工作流沿用原版本。[历史标签](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/tree/pre-zebra-reanalysis-20261120) · [最终验证](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44582662548)

| 变动 | 路径 |
|---|---|
| 修改 | `README.md`，`config/zebra_finch.json`，`notebooks/02_zebra_data.ipynb` |
| 修改 | `scripts/data/prepare_zebra.py`，`scripts/data/stage_zebra_sources.py`，`scripts/ci/run_zebra_finch.m`，`scripts/analysis/plot_zebra.py` |
| 新增 | `scripts/data/zebra_reanalysis.py`，`scripts/data/zebra_bursts.py`，`scripts/analysis/summarize_zebra.py` |
| 替换 | `data/fixtures/zebra_lman_population.mat`，`data/fixtures/zebra_lman_population.sha256` |
| 修改 | `data/fixtures/README.md`，`manifests/zebra-finch-preparation.json`，`manifests/experiments.json`，`results/zebra_finch/lman.json` |
| 新增 | `manifests/zebra-reanalysis-sources.json` |
| 修改 | `reports/03_zebra_dataset.md`，`reports/04_zebra_analysis.md`，`reports/05_final_summary.md` |
| 新增 | `reports/06_zebra_reanalysis_audit.md`，`reports/06_zebra_reanalysis_changelog.md` |
| 新增 | `figures/zebra_finch/` 下五张图，每张分别为 PDF、PNG、SVG：`lman_state_comparison`、`lman_firing_statistics`、`lman_burst_statistics`、`lman_burst_d2`、`lman_sensitivity` |
| 从主分支移除 | `figures/zebra_finch/lman_wake_sws.pdf`，`figures/zebra_finch/lman_wake_sws.png`，`figures/zebra_finch/lman_wake_sws.svg` |

`lman.json` 使用 `schema_version: 2`。`broad_state_estimates` 是原估计器的全态诊断，`diagnostics` 保存 AR、掩码与严格未配对结果，`resampling` 保存种子和逐次控制，`windows` 保存带录制时间的整窗结果。`state_metrics`、`contrasts` 和 `sws_association` 为相应汇总。`historical_reproduction` 供追溯历史基线。准备清单的 `states` 给出 v2 全态时长及边界；`recording_qc.recording_checks` 还保留原始加载器的录制角色和旧筛选时长字段用于来源核验。[结果 JSON](../results/zebra_finch/lman.json) · [准备清单](../manifests/zebra-finch-preparation.json)
