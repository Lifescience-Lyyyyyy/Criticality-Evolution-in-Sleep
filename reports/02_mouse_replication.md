# 第 2 阶段小鼠复现 · PASS

在论文指定的全部 13 个小鼠 V1 会话中，Wake 的 d₂ 均低于 NREM。Wake 和 NREM 的会话均值分别为 0.799720 和 2.401049 bits/s，NREM 减 Wake 的配对均差为 1.601329 bits/s，中位配对差为 1.529422 bits/s。[小鼠 MATLAB 运行](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44557104540)

## 数据与处理

本次分析采用作者发布的 Figshare v1 和论文列出的会话。13 个原始文件共 930,220,584 字节，在 Colab 下载并核对 MD5，再保存到 Drive。记录长于 30 分钟的 Wake 和 NREM 均通过检查。每个状态的群体脉冲按 40 ms 非重叠窗口计数，在原始 20 kHz 采样网格上处理。每个评分区间作为独立段输入，末端按作者示例的 `histc` 规则计数。传给 MATLAB 的压缩计数文件为 3,061,865 字节，写入 Drive 后独立下载，GitHub runner 校验同一 SHA-256。[Figshare v1](https://doi.org/10.6084/m9.figshare.32012181.v1) · [准备清单](../manifests/mouse-preparation.json)

按论文 §5.3，Fig. 4D 的状态汇总使用分段 Yule–Walker。这里的阶数 5 取自作者的鼠数据示例，对应 0.20 秒历史。我们复现的是主方向，图中数值对应本项目的计数实现和参数。以会话为配对单位。[PNAS 论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13732201/) · [原始工具箱](https://github.com/sam-sooter/prox_crit_toolkit/tree/f26c782dca3b1c9ab18cc92c6449591c44e0d507)

正式出版版的补充 PDF 未能通过公开端点读取。主分析的具体设置根据正式 Methods 与原始鼠数据示例确定，未主张图 4D 的逐点数值复现。

## 质量控制

未发现正时长状态区间重叠。YMV10 有一条 NREM 记录的起止均为 14,164 秒，主分析排除了这个零时长区间。保留它时，YMV10 NREM d₂ 从 2.5822250036 变为 2.5822270447 bits/s，配对方向一致。YMV02 的零时长区间属于本次比较未使用的 REM。

26 次主分析均输出有限、非负的 d₂，拟合 AR 过程均通过稳定性和残差标准差检查。只读取求解器诊断的复算与原始估计一致，所有 SQP 退出标志为 1。MATLAB 计算用时 25.50 秒。验收通过，准许开始斑胸草雀的配对资格核查。

[配对图 PDF](../figures/pnas/mouse_wake_nrem.pdf) · [逐会话结果与求解器诊断](../results/pnas/mouse.json)
