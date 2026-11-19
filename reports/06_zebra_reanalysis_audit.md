# 斑胸草雀重分析前审计

## 历史基线

重分析前的主分支为 `ae1ec2cf473cff6efa3a5296c74a766bba0eb37c`，已保存为 `pre-zebra-reanalysis-20261120`。旧版的斑胸草雀准备清单、结果、三份正文报告和 `lman_wake_sws` 图可通过该标签取回。旧版原始计数文件的 SHA-256 为 `9671a3a9fdf52d46635a7015552d7e8c340f432a32c455d8b30af66e7d0ad110`。旧版 MATLAB 使用全部非噪声分群，主体为 MUA，包括两类状态与各一次随机匹配。原配置的云端复算已触发，闸门 C 以结果文件逐项对照为准。[历史标签](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/tree/pre-zebra-reanalysis-20261120) · [历史复算](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44578774922)

**闸门 C PASS。** 从历史标签启动的 GitHub Actions 复算文件已独立下载并校验。按鸟、状态和变体对齐后，4 个历史主估计和 24 个优化诊断估计全部逐数相等，最大绝对差为 0。[复算记录](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/runs/44578774922)

## 状态与时钟

| 录制 | AP 元数据起始时刻 | 采样率（Hz） | AP 时长（秒） | 作者代码指向的角色 |
| --- | --- | ---: | ---: | --- |
| ZF1 23072019 | 2019-07-23 15:58:32 | 30,000 | 5373.330 | `initPaths` 载入的睡眠评分 |
| ZF1 23072019b | 2019-07-23 09:32:32 | 30,000 | 1903.792 | 声音检测程序的 `ZF1_awake` |
| ZF2 11112020 | 2020-11-11 20:42:53 | 30,000.331353 | 3855.720 | `initPaths` 载入的睡眠评分 |
| ZF2 11112020b | 2020-11-11 18:20:27 | 30,000.331353 | 7594.662 | 声音检测程序的 `ZF2_awake` |

ZF1 的会话 TSV 所示条件与程序路径、AP 时刻和评分长度交叉核验后可解决。亮灯数据的 `soundDetected` 标记手工设定阈值以上的麦克风 RMS，每秒有超过阈值即记 1。亮灯活动使用 1，亮灯非活动使用 0。亮灯非活动包含安静、休息或可能小睡的区间。声音来源可包括鸣唱和活动声。作者 Fig. 2F 源码将 `final_states_new=3` 映射为 SWS，4 为中间态，5 为 REM，99 表示未分类。SWS 同时要求 `new_mic=0`。不扩展 REM 或中间态。[Lorenz 等人](https://pmc.ncbi.nlm.nih.gov/articles/PMC12634191/) · [状态映射](https://gin.g-node.org/E4-Leblois/npx_sleep_2024/src/master/code/fig2/Figure2F_PieChart.m) · [声音检测源码](https://gin.g-node.org/E4-Leblois/npx_sleep_2024/src/master/code/fig3/SoundAndArtefactDetection/main_soundDetectionScript.m)

作者发布的 `ZF1_adResults.csv`、`ZF2_adResults.csv` 可定位亮灯录制的运动伪迹。写出代码使用零基 `sec` 列和 MATLAB 一基索引标记 `floor(artifactTimes)`，存在偏移一秒的风险。重分析使用发布秒和后移一秒的并集进行排除，独立保存两个时间解释。跨状态间断保持为独立序列。缺少 AP 原始运动伪迹掩码，作者的电压伪迹表来自 LFP。[运动伪迹源码](https://gin.g-node.org/E4-Leblois/npx_sleep_2024/src/master/code/fig3/SoundAndArtefactDetection/movArtDet.m)

## 单元审计与闸门 B

| 录制 | LMAN `good` | `mua` | `noise` | 其他 |
| --- | ---: | ---: | ---: | ---: |
| ZF1 睡眠 | 2 | 28 | 4 | 0 |
| ZF1 亮灯 | 1 | 34 | 0 | 0 |
| ZF2 睡眠 | 0 | 32 | 6 | 0 |
| ZF2 亮灯 | 0 | 31 | 1 | 0 |

各表还提供发放率、幅度与 Kilosort 污染值，完整记录级生理隔离和跨录制单元稳定性证据缺失。严格群体要求人工标记 `good`、正发放率、估计污染不超过 10%，并核查 1.5 ms 内的不应期违例比例。睡眠和亮灯每种状态至少需要两个有活动的合格单元。仅凭原始标签已可确定严格群体的三种配对均 **BLOCKED**。宽口径敏感性支路明确复用作者公开代码的 `group != noise` 且发放率大于零，允许 MUA。共享探针表支持脑区映射核验，无法据此跟踪单元身份。[作者加载源码](https://gin.g-node.org/E4-Leblois/npx_sleep_2024/src/master/code/UtilsMatlab/inputOutput/loadSpikes.m)

## 分析和验证要求

- 云端历史复算通过前不启动新 MATLAB 估计。优先复用 Drive 内已校验的 21 个文件，追加小型公开掩码与派生事件资料。
- 原始估计器仍以 40 ms、AR(5)、β=2、Yule–Walker 为参照，0.40 s AR(10) 作为单独敏感性。严格筛选和作者宽口径均记录；后者不代替受阻的严格主分析。以鸟为独立单位。
- 作者的群体爆发值使用连续两脉冲的瞬时发放率、同步窗内的爆发单元比例、Gaussian 平滑和逐录制第 99 百分位峰高。需核对公开代码中的 2,500 Hz 栅格、50 ms 窗、10 ms 步长及窗内超过一个 IFR 栅格点的判据。爆发检测在正式估计前用独立计算和原作者派生资料交叉核验。[爆发检测源码](https://gin.g-node.org/E4-Leblois/npx_sleep_2024/src/master/code/fig6/burstAnalysis.m) · [同步窗源码](https://gin.g-node.org/E4-Leblois/npx_sleep_2024/src/master/code/UtilsMatlab/computation/get_nActiveUnitsInSlidingWindow.m)
- 审计发现旧脚本在探针映射断言失败分支对错误文字误用了单目加号。旧结果的配对探针表一致，因此未进入该分支。修订入口需保留清晰的原始诊断。
