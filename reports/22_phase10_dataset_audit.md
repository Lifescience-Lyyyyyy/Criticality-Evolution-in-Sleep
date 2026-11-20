# Phase 10 数据审计

基线提交为 9e0ad8e。核验数据版本、来源指纹、时钟、脑区、干预和有效暴露后再选试点。门槛为每脑区至少 5 个排序单元、每个比较状态至少 120 秒 40 ms 完整时间箱，连续片段至少 1 秒。Wake 指可核验休息背景中的作者 Wake。开放场、训练和光刺激段排除在自然睡眠主比较之外。对未经验证的清醒类别统一标记为未分出静息 Wake。

| 数据来源 | 版本 | 发布动物 | 核验会话 | 会话资格 | 可用于状态比较的动物 |
|---|---|---:|---:|---|---:|
| 000056 ADN／PoSub | 0.261106.1949 | 7 | 42 | BLOCKED_METADATA 1／BLOCKED_QUALITY 2／READY_3STATE 39 | 7 |
| 000939 PoSub | 0.261016.2310 | 31 | 31 | BLOCKED_QUALITY 2／READY_2STATE 4／READY_3STATE 25 | 29 |
| 脑干 Neuropixels | Ver2 | 24 | 28 | PARTIAL 22／READY_3STATE 6 | 6 |
| 000978 CA1／PFC | 0.240511.0307 | 8 | 9 | BLOCKED_METADATA 9 | 0 |
| 000987 RSC | 0.261003.2342 | 3 | 14 | PARTIAL 14 | 0 |

[逐会话资格及来源指纹](../manifests/phase10_eligibility.csv) · [源版本和许可证](../manifests/phase10_sources.json)

## 000056 小鼠 ADN 和 PoSub

新版 42 场来自 7 只小鼠。逐簇比对单元脑区与电极位置，检查 20 kHz 原始簇的离散时钟、重复或逆序脉冲、记录支持和分期。论文报告三只动物同步录制 PoSub，Mouse20 的来源标注超出论文所述的植入范围，其名义 PoSub 簇未使用。未明确区分的脑区亦未作为 ADN 或 PoSub。Mouse12-120808 的簇 0／1 分别对应同一 shank 的独立作者簇 2／3，均严格递增并位于独立 LFP 记录内。该场重复标记的两类开放环境均在休息笼之外，实测重叠不涉及保留的休息片段。[逐簇脑区核验](../results/phase10/audit/000056_unit_region_validation.csv) · [新版转换](https://github.com/catalystneuro/buzsaki-lab-to-nwb/pull/101) · [版本数据](https://doi.org/10.48324/dandi.000056/0.261106.1949) · [原论文](https://doi.org/10.1038/nn.3968)

Mouse28-140314 和 Mouse32-140822 用 DatSegments.txt 的 1250 Hz 子记录边界作采集时间基准。两路 LED 坐标共享相机时钟，其发布记录跨度在两帧内与独立采集边界一致。脉冲在 20 kHz 时间格上，评分位于同一边界内。Mouse32-140820 缺可验证分期。Mouse17-130125 的合格 ADN 单元少于 5，Mouse17-130130 的脉冲延伸到已验证的记录边界之外。主分析使用 2019 年新版标记的 SleepScoreMaster 三态，限于 sleep_pre 和 sleep_post。论文版 Wake 为开放场觅食或迷宫，另用论文版 REM/SWS 做睡眠内敏感性核验，保留各源子记录边界。[完整逐会话审计](../results/phase10/audit/000056.json)

## 000939 小鼠 PoSub

逐一读取 31 只动物的 LFP 时钟、20 kHz 排序脉冲、home_cage 和作者阶段区间。SleepScoreMaster 与 TheStateEditor 使用 LFP 和伪 EMG。显式探索和光刺激段从主记录排除，刺激段两侧设 60 秒防护。超出独立记录边界的整条评分排除。与明确探索或干预重叠超过 1 秒的 REM／NREM 评分条目排除，并逐条列时间和标签。Wake 只使用家笼内已评分片段。分期的时钟、行为域矛盾和不足暴露按会话记录。整场不可用于比较的记录为 A3705-200306、A3711-200810b。A3705 的原始电生理和 LFP 均在 11168.136 秒结束，最后脉冲为 16198.560 秒。A3711 对应的边界为 9499.292 秒和 12176.531 秒。两场均无独立来源支持时轴修复。[完整时间异常表](../results/phase10/audit/000939_temporal_anomalies.csv) · [原始 MATLAB 配套数据](https://doi.org/10.6084/m9.figshare.24921252) · [版本数据](https://doi.org/10.48324/dandi.000939/0.261016.2310) · [原论文](https://doi.org/10.1038/s41593-024-01588-5)

## 脑干 Neuropixels

通过范围请求读取 22,619,797,342 字节归档的中央目录，确认 28 个 1k_train.npz 及配套排序表、脑区映射、评分、EEG／EMG。只为未干预 WT 子队列读取高时间分辨率数组。每个神经元的 1 ms 二进制计数与作者 2.5 秒汇总率逐箱比对，由 EEG／EMG 样本数和评分文件确认时间原点和末尾范围。40 ms 输入来自 1 ms 静态数组，保留睡眠片段与 KCUT 排除范围。版本和 ZIP 成员 CRC、提取文件 SHA-256 均列入审计。DL158、DL159 与 DL172 的 1 ms 阵列和 EEG 仅在末尾分别相差 50、108、30 个采样点。DL194 的 2.5 秒表多一条末尾不完整的记录。使用完整共同箱复核每个单元，全部保留的共同箱均与作者汇总相等，记录原始长度和舍弃的末尾。其他队列按光刺激协议标记，表中区分作者初始单元计数和经 1 ms 核查的 WT 单元。带不同颜色或重复后缀的同一动物只计一次。脑区由簇峰值通道与作者解剖图谱配对。[每只动物／核团／干预表](../results/phase10/audit/brainstem_nuclei_by_session.csv) · [作者加载代码](https://github.com/tortugar/Npx) · [版本数据](https://doi.org/10.5281/zenodo.19462601) · [原论文](https://doi.org/10.1038/s41593-026-02314-z)

## 000978 大鼠 CA1 和 PFC

9 个 NWB 档案，受试标识去重为 8 只。归档有排序脉冲、LFP 和视频速度。区段表缺三态评分，电极登记的逐单元映射为 unknown。作者示例读取本地 remeps、theta/delta 和 tetinfo 等 MAT 变量。公开脚本提供算法和参数，仓库文件清单没有匹配的逐动物派生数据。要验证重建，仍需可与 NWB 同钟对照的作者标注和 CA1／PFC 四极管映射。[作者代码](https://github.com/JadhavLab/REMHFOs) · [逐档案记录](../results/phase10/audit/000978.json) · [版本数据](https://doi.org/10.48324/dandi.000978/0.240511.0307) · [原论文](https://doi.org/10.7554/eLife.110795.2)

## 000987 小鼠 RSC

逐一检查 14 场、3 只小鼠。单元位置为 RSC。发布页提到睡眠分期，实见的睡眠阶段表只有 SWS，另列 UP/DOWN 和实验 sleep/wake 片段。Wake 指跟踪到的场地探索，家笼内 Wake 和 REM 的第二个可比状态缺失。按 PARTIAL 登记。此版本支持继续探索 SWS 单状态，以公开元数据目前无法做要求的自发状态配对。[实见分期表](../results/phase10/audit/000987.json) · [版本数据](https://doi.org/10.48324/dandi.000987/0.261003.2342) · [原论文](https://doi.org/10.1016/j.cub.2024.05.048)
