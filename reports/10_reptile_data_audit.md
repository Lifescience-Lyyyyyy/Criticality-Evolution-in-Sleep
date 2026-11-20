# Phase 7 · 爬行动物睡眠数据审计

**REP-0 完成。三个候选均未达到 READY。** 审计查到一条可公开访问的整夜宽带记录、一套需向作者请求的多通道数据，以及一套含 LFP 与作图衍生量的公开数据。原 Sooter 估计器在本项目中的输入是连续群体脉冲计数。三个公开入口均缺少满足这一规格的已核验区域群体脉冲与对齐睡眠时段。[来源与文件清单](../manifests/reptile-source.json) · [Sooter 等人](https://doi.org/10.1073/pnas.2526984123)

| 研究 | 分类 | 实际公开模态 | 动物与时长 | 状态及单元可用性 |
|---|---|---|---|---|
| Shein-Idelson 等，2016 | **PARTIAL** | 一只 Pogona 的单通道 DVR Neuralynx 宽带记录、日志、记录表与状态分析代码 | 公开示例 Lizard14，单次约 17.285 h。论文 Fig. 4B 同类状态曲线来自 4 只动物。 | 可依据作者 δ/β 程序识别 SWS／REM 型周期。分类单元所需多通道及连续脉冲时间未发布。 |
| Fenk 等，2023 | **RESTRICTED** | 公开代码含 NeuroNexus／Neuropixels 读取、睡眠特征和 spike-train 工具。神经数据需合理请求。 | 论文 Extended Data Fig. 3 的双侧屏状核分析为 4 只 Pogona，各 9 h；其他实验有各自样本。 | 原研究有 SWS／REMₚ 及排序结果；公开代码库未提供对齐的连续个体脉冲和状态表。 |
| Yamaguchi、Hatori 等，2024 | **LFP_ONLY** | DVR LFP 作图切片、δ/β、状态周期汇总，以及行为和体色表；MAT/XLSX。 | 论文总计 10 只 Pogona，实验跨 6–8 天；公开示例波形约 15 或 20 分钟。 | LFP 支持作者的周期分期方法。项目清单与数组审计均未见脉冲或分群文件。 |

## Shein-Idelson 等 · 单通道宽带保留

原研究网页已迁移至 Edmond DOI `10.17617/3.RNVKOY`。Dataverse API v1.0 清单含四个公开文件，许可 CC0 1.0。核验了 395 字节 README 和 218,086 字节分析代码压缩包的 MD5。通过 HTTP Range 读取 3,132,168,090 字节示例 ZIP 的中央目录，再抽查其中的 `runMe.m`、`sleepData.xlsx`、`data/CSC29.ncs` 头。示例只含一条连续采样通道，32 kHz，硬件滤波 0.1–9,000 Hz。这是宽带数据，示例文件内没有波形分群或已排序的脉冲时间表。按采样记录估算 17.285 h，与作者表一致。论文中 Fig. 4 的单元图依赖单独的 `spikeSorting.mat`，分析代码引用了该文件，公开 ZIP 的成员清单中没有它。要运行等同的群体计数分析，应先取得多单元排序导出、脑区映射和同一时钟的分期及缺口信息。用一条原始通道重新排序最多支持另行评估的单动物探索。[研究论文](https://pubmed.ncbi.nlm.nih.gov/27126045/) · [迁移后数据](https://doi.org/10.17617/3.RNVKOY) · [公开 API 清单](https://edmond.mpg.de/api/datasets/:persistentId/?persistentId=doi:10.17617/3.RNVKOY)

## Fenk 等 · 有代码，数据需请求

论文提供双侧屏状核以及相关 Imc、其他脑区和操作的结果。方法包括 IronClust 人工后处理，Neuropixels 支路使用 Kilosort2 和 Phy；这些说明了研究中具备单元排序。Data availability 明确写明数据按合理请求提供。查阅了公开 GitLab 代码树和 README，固定于 `7cb3a6ffa0cd7fc0359b74d2c9201114d83902dc`。代码包括原始格式读取、睡眠特征与脉冲工具，实际读取的数据不在仓库内。需要先获得同一动物的连续排序时间戳、单位隔离、脑区、两侧对时和 SWS／REMₚ 边界。若作者可提供排序与时钟资料，初始分析可从导出的脉冲表开始。[Nature 论文与数据声明](https://www.nature.com/articles/s41586-023-05827-w#data-availability) · [作者代码](https://gitlab.mpcdf.mpg.de/mpibr/laur/inter-hemispheric-rem)

## Yamaguchi、Hatori 等 · LFP 与出版衍生量

Mendeley v1 的公开文件 API 返回一个 16,580,730 字节的 `dataForFigures.zip`。先用 Range 检查中央目录，随后在审计中对该小型作图归档核验 SHA-256、27 个成员和 26 个 MAT 的顶层变量。Fig. 3、4 的 `repTrace` 含命名为 `raw_day*` 的代表性 LFP 和 δ/β 等值，波形为 1,200,001 或 900,001 个采样点。作者采集通带为 0.1–100 Hz，采样率 1 kHz。Fig. 2 的 `repDeltaBeta` 是派生特征，另有 REMS／SWS 周期与累计比例的汇总。这些资料可支持原论文的 LFP 节律描述。公开归档没有神经宽带阵列、分群 ID、校验过的脉冲时间或每只动物全程的同钟分期文件。资料许可为 CC BY 4.0。[Mendeley 数据](https://doi.org/10.17632/v7kkwgnyn5.1) · [文件 API](https://data.mendeley.com/public-api/datasets/v7kkwgnyn5/files?folder_id=root&version=1) · [论文方法](https://doi.org/10.1093/pnasnexus/pgad481)

以上分类适用于本次检查到的公开文件。爬行动物的排序、分箱与 d₂ 数值阶段停在资格审计。
