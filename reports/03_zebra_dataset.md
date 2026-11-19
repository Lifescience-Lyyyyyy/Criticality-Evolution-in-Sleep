# 第 3 阶段斑胸草雀数据资格 · PASS

ZF1 和 ZF2 有可配对的 LMAN 排序脉冲、作者状态评分，以及匹配的 Neuropixels 探针通道表。两只鸟适合进行探索性比较。资格要求为每种状态至少 120 秒完整评分时长、至少两个可用 LMAN 单元、评分时长与 AP 时钟一致，并可核对配对录制的探针映射。[Lorenz 等人 Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12634191/) · [公开数据](https://doi.org/10.12751/g-node.vsxyi8)

| 鸟 | SWS 录制 | SWS 时长（秒） | SWS 单元 | 亮灯活动录制 | 活动时长（秒） | 活动单元 | 资格 |
| --- | --- | ---: | ---: | --- | ---: | ---: | --- |
| ZF1 | ses-23072019 | 3918 | 30 | ses-23072019b | 1368 | 35 | 通过 |
| ZF2 | ses-11112020 | 2997 | 32 | ses-11112020b | 1106 | 31 | 通过 |
| ZF3 | ses-28022020 | — | — | ses-28022020b | — | — | 公开睡眠派生目录有分群元数据，缺配套时刻及分群数组 |
| ZF4 | ses-14032022 | — | — | — | — | — | 已核查的目录缺成对排序数组 |
| ZF5 | ses-23102020 | — | — | ses-21102020 | — | — | 单独亮灯录制仅发布 LFP。睡眠录制的附表有亮灯标记，但未确认可作为作者标注清醒的配对样本 |
| ZF6–7 | 原作者记录 | — | — | — | — | — | 公开参与者清单标为 LFP only |
| ZF8–10 | — | — | — | — | — | — | 作者因可能残留麻醉效应排除 |

SWS 使用作者最终评分 `final_states_new=3` 且无麦克风标记的秒。觉醒侧只使用作者指定清醒录制中 `soundDetected=1` 的亮灯活动秒。睡眠侧的深夜声音、活动及未分类秒未当作清醒。两个状态来自不同录制，无法确认逐单元追踪，也无法排除录制时间、群体取样及活动本身的影响。[作者的状态分析方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC12634191/) · [评分与处理代码](https://gin.g-node.org/E4-Leblois/npx_sleep_2024/src/master/code)

ZF1 的 BIDS 会话说明与 AP 创建时刻及配套评分文件有对调迹象。用作者的 `initPaths.m`、清醒评分文件名以及 AP 时长交叉确定角色。光照或日期单独不足以定义状态。正式使用的 LMAN 通道范围来自作者的 `birdParams.mat`，ZF1 为 240–280，ZF2 为 213–260。逐录制按实际 AP 采样率将排序索引换算为秒。每对 AP 探针通道表一致。[会话说明](https://gin.g-node.org/E4-Leblois/npx_sleep_2024/src/master/dataset/sub-ZF1) · [作者的加载函数](https://gin.g-node.org/E4-Leblois/npx_sleep_2024/src/master/code/UtilsMatlab/inputOutput)

GIN 服务器对 Colab 请求返回 HTTP 403。使用 WSL 的公开下载端点取得 21 个文件，共 224,335,888 字节，并在前后核对 GIN 提交。文件传入用户的 Drive，按 WSL 生成的 SHA-256 逐个读回校验，随后在 Colab 进行数据质控和格式转换。完整 GIN 归档约 1.2 TiB。分析使用公开的已排序脉冲。

四个录制的排序索引按时间排列，时刻与分群数组行数一致。组内同一采样点同一单元的重复只在 ZF1 SWS 发现 3 次，计数准备时各保留一次。区分同一采样点不同单元的脉冲。四个录制内全单元相邻脉冲的最大间隔为 0.011–0.170 秒。ZF1 睡眠排序文件末次脉冲约比 AP 时长早 10 秒，完整录制的个体单元稳定性仍缺乏独立证据。

[四个候选录制的 QC 和文件清单](../manifests/zebra-finch-audit.json) · [Colab 数据准备](../notebooks/02_zebra_data.ipynb)
