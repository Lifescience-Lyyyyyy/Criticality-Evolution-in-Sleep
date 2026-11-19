# 睡眠中的神经临界性演化

用 Sooter 等人的 MATLAB 原始估计器复现小鼠 V1 的清醒与 NREM 比较，再评估同一方法能否用于斑胸草雀。

## 执行

- 在 [MATLAB PNAS Replication](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/workflows/matlab.yml) 手动选择 `smoke`、`synthetic`、`mouse` 或 `zebra`。对应的 `matlab-smoke-*`、`matlab-synthetic-*`、`matlab-mouse-*`、`matlab-zebra-*` Git 标签也可显式触发。结果从运行页的 `matlab-results` 下载。
- 在 Colab 运行 [Drive 初始化 notebook](notebooks/00_drive_setup.ipynb)，将实验数据和归档保存在 `MyDrive/Research/Criticality_Evolution/`。
- 在合成验证通过后运行 [小鼠数据准备 notebook](notebooks/01_mouse_data.ipynb)。Git 中的小鼠输入由公开数据生成，包含 40 ms 计数与分段信息。原始 Figshare 文件归档在 Drive。
- 小鼠复现通过后，在 WSL 用 `python3 scripts/data/stage_zebra_sources.py` 将选定的公开 GIN 文件传到已挂载 Drive 的 Colab 会话。然后用 [斑胸草雀数据准备 notebook](notebooks/02_zebra_data.ipynb) 核验并导出小型 LMAN 计数文件。
- [中文研究报告](reports/05_final_summary.md) 汇总全部阶段、配对图和复现实验入口。数据版本、校验值和运行记录见 [来源清单](manifests/sources.json) 与 [实验索引](manifests/experiments.json)。

主要来源为 [Sooter 等人的 PNAS 论文](https://doi.org/10.1073/pnas.2526984123)、[作者的 MATLAB 工具包](https://github.com/sam-sooter/prox_crit_toolkit) 和 [Lorenz 等人的斑胸草雀录音](https://doi.org/10.12751/g-node.vsxyi8)。
