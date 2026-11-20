# Phase 9 数据资格审计

## 采集和分析门槛

先读取论文、官方数据说明和文件目录，再做逐会话只读检查。状态对比要求可核对的排序脉冲、状态与神经信号共用时钟、明确脑区以及同一会话内同一组单元。预设每区至少 5 个可用单元，每个待比较状态至少 120 秒有效 40 ms 箱，单连续片段至少 1 秒。保留作者状态名称和边界。跨区只用共同标注时间。[逐会话资格表](../manifests/phase9_dataset_eligibility.csv) · [来源及文件指纹](../manifests/phase9_sources.json)

## 原始 PNAS 大鼠 V1

状态为 `ACCESS_REQUEST_REQUIRED`。论文公开了 Methods、补充结果、估计代码和小鼠数据。原始大鼠的可重算输入需要联系作者。申请范围为现成排序脉冲或分段 40 ms 群体计数及同期阶段评分。[核查与申请条件](18_pnas_rat_v1_access.md)

## CRCNS hc-24 · 丘脑、mPFC、CA1

状态为 `ACCESS_REQUIRED`。官方说明列三场、三只大鼠，时长分别约 52、69 和 93 分钟，整个 MATLAB 档案约 101 MB。登录页要求 CRCNS 账号。按用户当前安排，本轮完成文档审计后暂缓档案内核验。官方说明区分 `cl-` 排序单元和来自一根四极管的多单元活动，排序文件为五列，第五列为秒。LFP 文档采样率为 600 Hz。论文用 mPFC δ 功率和 CA1 涟漪频带功率识别部分 NREM。三态标签、每区单元数、同步范围和可用时长须在实际档案中核验。丘脑有单元记录的公开说明无法证明 CA1 或 mPFC 具有足量单元。[CRCNS 文档](https://crcns.org/files/data/hc-24/crcns_hc-24_data_description.pdf) · [Varela 与 Wilson 的论文](https://doi.org/10.7554/eLife.48881)

在浏览器打开官方下载页，拥有账号时直接登录，首次使用者按 CRCNS 页面办理注册。下载档案后将文件交给当前项目的本地输入目录。需要核验的材料包括 `Primary Data`、`Derived Data` 和 `Sample Scripts`，尤其是时间轴、单元区域目录、人工或算法评分定义、连续记录边界。账号信息保留在浏览器中。论文使用需引用论文和 DOI，作者欢迎事先沟通。[官方登录与下载](https://download.crcns.org/hc-24) · [CRCNS 使用说明](https://crcns.org/data-sets/hc/hc-24/about-hc-24)

## CRCNS pfc-8 · mPFC 和纹状体

16 个已发布会话目录均可匿名访问。可见动物 ID 有 R205、R209、R268、R287、R294 共 5 个，官方说明则描述采集了 6 只。保留实见 ID，缺少的第 6 只无法对应到发布文件。文件目录总计 33,116,450,597 字节，其中所有 LFP 共 17,585,969,664 字节。审计读取 XML、事件文本和 CLU 整数标签，并对每会话的一份排序 RES 的开头和末尾取样。XML 的 19,531.25 Hz 脉冲计数时钟与 1,250 Hz LFP 的时长及毫秒事件边界相容。每会话均有 POST，R268-20140930 无标为 PRE 的休息片段。所有 `MIX`、`ShiftTasks` 等任务段单独标识并排除。[数据门户](https://crcns.org/data-sets/pfc/pfc-8/about-pfc-8) · [正式论文](https://doi.org/10.1016/j.cub.2021.10.007)

**16 个会话的状态对比均为 `BLOCKED`。** `.cat.evt` 只提供实验段边界，`.all.evt` 和 `.cue.evt` 是数值和任务提示事件。发布目录中没有三态评分边界。XML 的电极分组支持逐会话统计候选排序簇，但通道备注使用固定的 mPFC／腹侧纹状体模板，未给出可逐四极管复核的 dmSTR 身份。可用的 LFP 和位置提供重新评分的输入。需要实验者评分，或有独立可验证的分期标准及验证样本，才能用来比较 Wake、NREM、REM。还需索取已核验的四极管到 vSTR／dmSTR 对照表。使用这些数据筹备发表时，应按 CRCNS 要求联系 Sidney Wiener 或 Michaël Zugaro 讨论作者安排。当前工作只进行审计。[数据和使用条件](https://crcns.org/files/data/pfc-8/crcns_pfc-8_data_description.pdf)

## CRCNS hc-23 · MEC 浅层 II／III

公开数据提供 6 只大鼠的 29 场 Day 1 Overnight。采用 Figshare 处理版 v1，216,394,648 字节，已核对发布 MD5 和复算 SHA-256。用 HTTP 字节范围读取配套原始发布包内的 10,978 份仅含排序脉冲的 Overnight \`.t\`，逐文件校验 ZIP CRC。原发布包为 493,578,348 字节。该核验读取 24,063,943 字节，包括目录和所需范围。[处理版](https://doi.org/10.6084/m9.figshare.7973339.v1) · [排序脉冲档案](https://doi.org/10.6084/m9.figshare.7966679.v1) · [数据说明](https://doi.org/10.6084/m9.figshare.7966664.v1)

作者先用两位评分者的视频一致意见确定睡眠，再用 LFP 的 θ／δ 功率划分 REM 与 NREM。用于分析的是 \`Day1/Overnight/bout/remTimes\` 和 \`nRemTimes\` 的秒级区间。按四极管、簇号和网格类型核对处理版脉冲与原始 \`.t\` 的 100 μs 整数时钟。记录保留各 \`Sleep\` 片段，在片段内合并同态重叠区间，剔除异态共同区间及无已核验脉冲支持的部分。80 个片段的异态共同区间共 243.5 秒。筛除无效或无法明确对应的单元身份，并使用同一会话跨所有保留片段都有证据的单元集合。[作者分期与单位稳定性方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC7412059/)

21 场有足量稳定单元与 NREM／REM，覆盖全部 6 只动物，标为 \`PARTIAL\`，可估计 REM − NREM。8 场的可用单元不足 5 个，标为 \`BLOCKED\`。公开 Overnight 档案由多个不连续的睡眠片段组成，Wake 及整夜原始 LFP 未发布。虽然 README 泛称附 \`.nvt\`，Day 1 Overnight 的归档索引中未见该格式。可核验的是作者标签与公开排序脉冲的时间一致性，无法独立重做 LFP 分期。Open Field 和第二天记录未纳入自发睡眠对比。许可为 CC BY 4.0。[逐会话审计](../manifests/phase9_dataset_eligibility.csv) · [Figshare 公开集合](https://doi.org/10.6084/m9.figshare.c.4496375)
