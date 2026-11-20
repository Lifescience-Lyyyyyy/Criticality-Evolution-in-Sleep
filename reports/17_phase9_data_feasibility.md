# Phase 9 数据资格审计

## 采集和分析门槛

先读取论文、官方数据说明和文件目录，再做逐会话只读检查。状态对比要求可核对的排序脉冲、状态与神经信号共用时钟、明确脑区以及同一会话内同一组单元。预设每区至少 5 个可用单元，每个待比较状态至少 120 秒有效 40 ms 箱，单连续片段至少 1 秒。保留作者状态名称和边界。跨区只用共同标注时间。[逐会话资格表](../manifests/phase9_dataset_eligibility.csv) · [来源及文件指纹](../manifests/phase9_sources.json)

## CRCNS hc-24 · 丘脑、mPFC、CA1

已核验用户提供的官方归档及数据说明。两个文件匹配发布的 MD5，归档的 99 个条目通过 ZIP CRC 检查。三场各来自一只大鼠，本报告以 Rat-1 至 Rat-3 对应会话编号。所有会话的 LFP 时间轴连续，按 600 Hz 采样。每个会话的作者评分文件只含 NREM 区间，Wake 和 REM 无可用标签。[数据说明](https://crcns.org/files/data/hc-24/crcns_hc-24_data_description.pdf) · [原论文](https://doi.org/10.7554/eLife.48881)

按来源文件名映射脑区和单元，只使用作者排序的 cl- 簇。脉冲取第五列的秒级时间，与 LFP 共用时轴。落在 LFP 支持之外的脉冲排除，含非递增时间戳的簇按既定时间 QC 排除。第一场通过筛选的 CA1 5 个和 mPFC 6 个单元可进行同一 NREM 时段的描述性比较。第二场只有 2 个排序丘脑单元，低于 5 个单元门槛。第三场为跑后休息，CA1 发布的 8 个簇有 7 个通过检验，满足单区 NREM 分析。公开材料中的多单元信号单独记录。[逐会话资格表](../manifests/phase9_dataset_eligibility.csv)

## CRCNS hc-23 · MEC 浅层 II／III

公开数据提供 6 只大鼠的 29 场 Day 1 Overnight。采用 Figshare 处理版 v1，216,394,648 字节，已核对发布 MD5 和复算 SHA-256。用 HTTP 字节范围读取配套原始发布包内的 10,978 份仅含排序脉冲的 Overnight \`.t\`，逐文件校验 ZIP CRC。原发布包为 493,578,348 字节。该核验读取 24,063,943 字节，包括目录和所需范围。[处理版](https://doi.org/10.6084/m9.figshare.7973339.v1) · [排序脉冲档案](https://doi.org/10.6084/m9.figshare.7966679.v1) · [数据说明](https://doi.org/10.6084/m9.figshare.7966664.v1)

作者先用两位评分者的视频一致意见确定睡眠，再用 LFP 的 θ／δ 功率划分 REM 与 NREM。用于分析的是 \`Day1/Overnight/bout/remTimes\` 和 \`nRemTimes\` 的秒级区间。按四极管、簇号和网格类型核对处理版脉冲与原始 \`.t\` 的 100 μs 整数时钟。记录保留各 \`Sleep\` 片段，在片段内合并同态重叠区间，剔除异态共同区间及无已核验脉冲支持的部分。80 个片段的异态共同区间共 243.5 秒。筛除无效或无法明确对应的单元身份，并使用同一会话跨所有保留片段都有证据的单元集合。[作者分期与单位稳定性方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC7412059/)

21 场有足量稳定单元与 NREM／REM，覆盖全部 6 只动物，标为 \`PARTIAL\`，可估计 REM − NREM。8 场的可用单元不足 5 个，标为 \`BLOCKED\`。公开 Overnight 档案由多个不连续的睡眠片段组成，Wake 及整夜原始 LFP 未发布。虽然 README 泛称附 \`.nvt\`，Day 1 Overnight 的归档索引中未见该格式。可核验的是作者标签与公开排序脉冲的时间一致性，无法独立重做 LFP 分期。Open Field 和第二天记录未纳入自发睡眠对比。许可为 CC BY 4.0。[逐会话审计](../manifests/phase9_dataset_eligibility.csv) · [Figshare 公开集合](https://doi.org/10.6084/m9.figshare.c.4496375)
