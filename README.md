# 睡眠中的神经临界性演化

用 Sooter 等人的 MATLAB 原始估计器复现小鼠 V1 的清醒与 NREM 比较，再评估同一方法能否用于斑胸草雀。

## 执行

- 在 [MATLAB PNAS Replication](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/workflows/matlab.yml) 手动选择 `smoke`、`synthetic` 或 `mouse`。对应的 `matlab-smoke-*`、`matlab-synthetic-*`、`matlab-mouse-*` Git 标签也可显式触发。结果从运行页的 `matlab-results` 下载。
- 在 Colab 运行 [Drive 初始化 notebook](notebooks/00_drive_setup.ipynb)，将实验数据和归档保存在 `MyDrive/Research/Criticality_Evolution/`。
- 在合成验证通过后运行 [小鼠数据准备 notebook](notebooks/01_mouse_data.ipynb)。Git 中的小鼠输入由公开数据生成，包含 40 ms 计数与分段信息。原始 Figshare 文件归档在 Drive。
- 验收记录见 [可行性审计](reports/00_feasibility.md)、[合成验证](reports/01_synthetic.md) 和 [小鼠复现](reports/02_mouse_replication.md)。数据版本和校验值见 [来源清单](manifests/sources.json)。

主要来源为 [Sooter 等人的 PNAS 论文](https://doi.org/10.1073/pnas.2526984123)、[作者的 MATLAB 工具包](https://github.com/sam-sooter/prox_crit_toolkit) 和 [Lorenz 等人的斑胸草雀录音](https://doi.org/10.12751/g-node.vsxyi8)。
