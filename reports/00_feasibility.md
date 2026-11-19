# 第 0 阶段可行性审计

审计日期为 2026-11-19，所有 MATLAB 执行均使用 GitHub 托管 runner。

| 检查 | 状态 | 证据 |
| --- | --- | --- |
| MATLAB 执行与 Artifact | PASS | R2024b Update 10 已加载 Optimization Toolbox 和 Parallel Computing Toolbox。作者的示例返回 `exit_status=0` 及 `d₂=0.03465842415367701` bits/s。Artifact 670 字节，已下载并核对 SHA-256。[运行 44552272750](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44552272750) |
| PNAS 官方代码 | PASS | 固定至 `f26c782dca3b1c9ab18cc92c6449591c44e0d507`。从 `demo/demo.mlx` 原样提取首段合成示例与本地函数，执行原始 `src/calc_db.m`。[官方版本](https://github.com/sam-sooter/prox_crit_toolkit/tree/f26c782dca3b1c9ab18cc92c6449591c44e0d507) |
| 小鼠数据访问 | PASS（清单及格式说明） | Figshare v1 包含论文使用的全部 13 个会话，总计 930,220,584 字节。说明文件给出以秒为单位的 `st` 与 `WAKE`、`NREM`、`REM` 区间，以及单元标签 `clu`。实际会话质控安排在合成验收后。[Figshare v1](https://doi.org/10.6084/m9.figshare.32012181.v1) |
| 斑胸草雀数据访问 | PASS（元数据） | GIN 仓库包含 BIDS 风格的会话清单、派生分群文件及行为状态评分。实际配对资格要验证脑区映射、时钟与可用状态时长。[GIN](https://gin.g-node.org/E4-Leblois/npx_sleep_2024) |
| Google Drive 持久化 | PASS | Colab CPU 启动后经用户授权挂载。项目创建了计划目录，原子写入 153 字节并读回。后续独立下载得到相同 SHA-256。挂载报告可用 230,258,921,472 字节。 |
| 跨环境传输 | PARTIAL | Drive → Colab → 工作区已核验同一份 153 字节输入。runner 上的校验步骤已准备。 |

GitHub 仓库为公开仓库。MathWorks 的自动 CI 许可已支持上述实际运行。原始任务超时为 120 分钟，本次 runner 在安装前报告根卷可用 91 GiB。工作流可手动触发，也可通过阶段标签启动。[MathWorks 许可说明](https://github.com/matlab-actions/setup-matlab#licensing)

官方 `demo_crit_db.m` 设置 `b=6`，同时示例要求 Yule–Walker 拟合计算误差条。工具包不支持这个误差条组合。实际冒烟测试使用官方 Live Script 的首段，按 `b=2`、AR(3)、10,000 个采样点、步长 1 秒执行，随机种子固定为 2026。

原论文 Fig. 4D 对同一状态的连续区间计算长度加权自相关，并用 Yule–Walker 拟合。区间之间不拼接。Fig. 4H 的滑窗用 MLE，复现 Fig. 4D 时以论文方法为准。[PNAS Methods §5.3](https://pmc.ncbi.nlm.nih.gov/articles/PMC13732201/)

DANDI 000166 的版本 `0.250624.0434` 有 19 个 NWB 文件，总计约 787 GB。小鼠复现优先使用作者发布的 Figshare 文件。斑胸草雀归档约 1.2 TiB，将按会话挑选排序结果、脑区映射和状态文件。实际输入需保留来源和校验值，超出 Git 体量限制的持久文件保存在 Drive。[DANDI](https://dandiarchive.org/dandiset/000166) · [斑胸草雀数据 DOI](https://doi.org/10.12751/g-node.vsxyi8)

官方冒烟测试的验收状态为 **PASS**。仅使用合成输入的估计器验证可以继续。真实数据计算需先通过各自的输入与传输验收。
