# Sooter 等人的原始大鼠 V1 数据获取

**状态为 `ACCESS_REQUEST_REQUIRED`。** PNAS 的数据声明把公开小鼠数据和原始大鼠记录分开。大鼠源文件为超过 10 TB 的原始宽带电信号，需向 Keith B. Hengen 申请。现有公开材料未提供能与 Fig. 4D 大鼠点逐一对应的排序脉冲或分段群体计数、动物与会话清单及同步睡眠评分。本项目的大鼠 V1 复现暂停在原数据获取门槛。[论文与数据声明](https://pmc.ncbi.nlm.nih.gov/articles/PMC13732201/)

## 已核对的公开材料

- 正式论文 §5.1 描述家笼自由行为的 10–14 天大鼠 V1 记录、25 kHz 电信号和同步视频，并引 Xu 等人的采集和排序方法。§5.3 和 Fig. 4 的说明给出 40 ms 群体脉冲箱。Fig. 4D 按会话和状态分别计算各连续段 ACF，加权合并后用 Yule–Walker 拟合。Fig. 4H 的 60 秒滑窗采用 MLE。两者各有独立的复现规格。[论文 Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC13732201/)
- 通过 Europe PMC 取得 PNAS 补充 PDF。§3 和图 S6–S7 添加了大鼠雪崩分布及模型比较，图 S8 给出时间分箱和 AR 阶数的敏感性。补充包中的图像和 PDF 缺少可重算的数据序列与 Fig. 4D 资产映射。[期刊补充附件](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13732201/supplementaryFiles)
- 官方工具包固定在 `f26c782dca3b1c9ab18cc92c6449591c44e0d507`。固定提交与当前上游树一致，含原估计器、模拟示例、小鼠示例 `demo/ephys_data/YMV15.mat`。论文指向的 Figshare 项目 `32012181.v1` 也明确为小鼠数据。[官方工具包](https://github.com/sam-sooter/prox_crit_toolkit/tree/f26c782dca3b1c9ab18cc92c6449591c44e0d507) · [小鼠 Figshare](https://doi.org/10.6084/m9.figshare.32012181.v1)
- Xu 等人的相关 V1 论文同样注明宽带需申请。其公开代码及附带图源表是独立论文的分析材料，目前没有能将其中记录与 Sooter Fig. 4D 的动物、会话和状态片段逐一对应的清单。[Xu 等人数据可用性](https://www.nature.com/articles/s41593-023-01536-9#data-availability) · [Xu 等人代码及图源表](https://github.com/hengenlab/Sleep_restores_criticality/tree/0c6509f1969cbb6e9d07d7c0b2ff0bbc0328f347)

## 申请与恢复条件

优先申请 Fig. 4D 用的排序脉冲及单元 ID，或按动物、会话、连续片段提供的 40 ms 非重叠群体计数。两种格式均需同期 Wake／NREM／REM 标注、时间原点、片段边界、单位筛选和 Fig. 4D 对应关系。若公开补充仍不足以确认大鼠的 AR 阶数、加权和排除规则，随数据附短方法说明或原分析脚本。若只有逐段计数，可复现状态级估计，单元数和单元质量对照需要逐单元数据。

作者授权后先核对许可、文件哈希、上下游时钟、标签和片段，按 Fig. 4D 规格恢复复现。可独立做项目预设的 AR(5)／AR(10) 比较。原始大鼠 V1 记录和其他来源的 V1 或额叶结果保持各自来源标识。[待发送的英文申请草稿](pnas_rat_v1_data_request_email.md)
