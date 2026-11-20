# 睡眠中的神经临界性演化

用 Sooter 等人的 MATLAB 工具包复现小鼠 V1 的 Wake–NREM 比较。斑胸草雀 LMAN 重分析按亮灯活动、亮灯非活动和评分 SWS 组织。严格单元配对缺少合格群体，当前的鸟数据结果用于诊断。历史两态分析已由 [v2](reports/04_zebra_analysis.md) 替换，可从 [历史标签](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/tree/pre-zebra-reanalysis-20261120) 复算。

## 执行

- 在 [MATLAB PNAS Replication](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/workflows/matlab.yml) 手动选择 `smoke`、`synthetic`、`mouse`、`zebra` 或 `sleep`。对应的 `matlab-smoke-*`、`matlab-synthetic-*`、`matlab-mouse-*`、`matlab-zebra-*`、`matlab-sleep-*` 标签也能触发。运行完成后下载 `matlab-results`。
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

## Phase 7 数据审计

RA 的 OSF 发布含四只斑胸草雀的双侧排序脉冲。睡眠时段、鸣唱与重放缺少可对齐到神经时钟的标注，RA-0 受阻。三个爬行动物数据源的公开资料也尚不足以开展本项目的群体脉冲 d₂ 分析。

- [RA 可行性](reports/07_ra_feasibility.md) · [Phase 7 科学结论与验证顺序](reports/08_ra_analysis.md) · [RA–LMAN 可比范围](reports/09_ra_lman_comparison.md)
- [爬行动物数据审计](reports/10_reptile_data_audit.md) · [RA 来源清单](manifests/ra-source.json) · [爬行动物来源清单](manifests/reptile-source.json)
- 重跑 OSF 只读清单和文件头检查可用 `python3 scripts/data/inspect_ra.py`。

## Phase 8 自发睡眠

新增大鼠家笼额叶、同录的大鼠 CA1–BLA，以及睡眠与刺激边界明确的小鼠家笼 CA1。各来源以发布版固定，保留作者 Wake 标签。动物是独立比较单位，区域同录对照按 PRE／POST 分开。

- [资格审计](reports/11_sleep_datasets_feasibility.md) · [真实数据试跑](reports/12_sleep_pilot.md) · [区域内结果](reports/13_sleep_within_region.md) · [同期区域对照](reports/14_sleep_cross_region.md) · [NREM 涟漪](reports/15_sleep_events.md) · [科学总结](reports/16_sleep_project_summary.md)
- [逐资产资格 CSV](manifests/sleep_dataset_eligibility.csv) · [来源与版本](manifests/sleep_sources.json) · [运行与 SHA-256](manifests/sleep_experiments.json) · [完整 Drive 归档](https://drive.google.com/drive/folders/1sy_0Tjej6YiUCTJZ-AL9_GUtsnE5fXPL)
- 在 Colab 安装 `remfile==0.1.16`，从仓库根目录运行 `python scripts/data/audit_sleep_dandi.py --out <Drive 项目>/manifests/sleep_dandi`。随后用 `scripts/data/prepare_sleep_dandi.py` 的 `--audit` 和 `--dandiset` 选择合格队列，也可用 `--asset` 限定单会话。新增 `sleep` 工作流阶段读取经过校验的计数与收据。仓库提供 ACC 试跑 fixture，完整批次保存在项目 Drive，对应 Git 标签保留每次运行的准确输入和代码。
- `--paired-regions` 选择同录 CA1–BLA。对应 `--events` 和 `--controls` 分别准备涟漪和重复对照。`--rem-from-fixtures` 从验证过的三套计数 fixture 匹配三态时长。用 `scripts/analysis/summarize_sleep_dandi.py` 汇总 MATLAB JSON，以 `scripts/analysis/plot_sleep_dandi.py` 生成逐动物图。
- 估计使用原始 MATLAB 版本、40 ms、AR(5)、β=2 和 Yule–Walker，0.40 秒跨度及适用的额外控制另存。只读 NWB 排序脉冲及小型注释表，连续区间分别传给估计器。

[Sooter 等人的 PNAS 论文](https://doi.org/10.1073/pnas.2526984123) · [原始 MATLAB 工具包](https://github.com/sam-sooter/prox_crit_toolkit) · [Lorenz 等人的斑胸草雀录制](https://doi.org/10.12751/g-node.vsxyi8)
