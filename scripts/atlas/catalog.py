"""Anatomical crosswalk and dataset context for the Phase 11 read-only export."""
BASE = "https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep"
SOURCE_COMMIT = "a834a85ce765dfb2fb77f267ffd3410c64a1f339"
PIN = "f26c782dca3b1c9ab18cc92c6449591c44e0d507"
STATES = ("Wake", "NREM", "REM")
PAIRS = {"NW": ("Wake", "NREM"), "RN": ("NREM", "REM"), "RW": ("Wake", "REM")}
FIT_FILES = (
    "results/sleep_dandi/000552/mouse-ca1.json",
    "results/sleep_dandi/000041/cortex.json",
    "results/sleep_dandi/000061/cross.json",
    "results/sleep_crcns/hc-23/cohort.json",
    "results/sleep_crcns/hc-24/pilot.json",
    "results/sleep_crcns/hc-24/post-run.json",
    *("results/phase10/fits/" + name + ".json" for name in (
        "a_cohort", "posub-adn-supplement", "c_cohort", "c-rescue_cohort", "quiet-wake"
    )),
)
V1_FILE = "results/pnas/mouse.json"
# display abbreviation, Chinese name, English name, system, Allen CCF acronym, description
MOUSE = (
    ("V1", "初级视觉皮层", "Primary visual cortex", "感觉皮层", "VISp", "处理由外侧膝状体传来的视觉输入。"),
    ("CA1", "海马 CA1 区", "Hippocampal field CA1", "海马系统", "CA1", "海马的主要输出区之一，参与空间表征和记忆。"),
    ("ADN", "丘脑前背核", "Anterodorsal thalamic nucleus", "丘脑", "AD", "头方向网络的丘脑节点，与海马旁区域相连。"),
    ("PoSub", "后下托", "Postsubiculum", "海马旁系统", "POST", "头方向网络的皮层节点，关联空间定向与环境线索。"),
    ("PAG", "导水管周围灰质", "Periaqueductal gray", "中脑", "PAG", "环绕中脑导水管，与防御、自主调节和睡眠状态相关。"),
    ("P", "脑桥（作者 P 标签）", "Pons", "脑干", "P", "脑干的脑桥分区。P 是源数据的区域标签，覆盖范围较广。"),
    ("DR", "中缝背核", "Dorsal raphe nucleus", "中脑", "DR", "脑干中缝核团之一，包含广泛投射的血清素能神经元。"),
    ("LDT", "背外侧被盖核", "Laterodorsal tegmental nucleus", "脑桥被盖", "LDT", "含多类神经元，参与觉醒与睡眠调节。"),
    ("MRN", "中脑网状核", "Midbrain reticular nucleus", "中脑", "MRN", "中脑网状结构的一部分。"),
    ("PRNc", "尾侧脑桥网状核", "Caudal pontine reticular nucleus", "脑桥", "PRNc", "脑桥网状结构的尾侧分区。"),
    ("PRNr", "脑桥网状核（作者 PRNr 标签）", "Pontine reticular nucleus", "脑桥", "PRNr", "保留作者 PRNr 标签，CCF 使用 PRNr 对应分区。"),
    ("RPO", "脑桥中缝核", "Nucleus raphe pontis", "脑桥", "RPO", "位于脑桥正中附近的中缝核团。"),
    ("CS", "中央上中缝核", "Superior central nucleus raphe", "脑桥", "CS", "中缝核团之一。"),
    ("SCdw", "上丘运动深白层", "Superior colliculus, motor related, deep white layer", "中脑", "SCdw", "上丘运动相关区的深白层。"),
    ("SCdg", "上丘运动深灰层", "Superior colliculus, motor related, deep gray layer", "中脑", "SCdg", "上丘运动相关区的深灰层。"),
    ("MB", "中脑（作者 MB 标签）", "Midbrain", "中脑", "MB", "保留作者大区标签，无法细分为具体核团。"),
    ("RSC", "压后皮层", "Retrosplenial cortex", "联合皮层", "RSP", "参与空间定向与记忆。现有审计只确认慢波睡眠资料。"),
)
RAT = (
    ("Frontal", "额叶皮层队列", "Frontal cortex", "皮层"),
    ("mPFC", "内侧前额叶皮层", "Medial prefrontal cortex", "皮层"),
    ("ACC", "前扣带皮层", "Anterior cingulate cortex", "皮层"),
    ("OFC", "眶额皮层", "Orbitofrontal cortex", "皮层"),
    ("MotorCtx", "运动皮层", "Motor cortex (source label)", "皮层"),
    ("CA1", "海马 CA1 区", "Hippocampal field CA1", "海马系统"),
    ("BLA", "基底外侧杏仁核", "Basolateral amygdala", "杏仁核"),
    ("MEC II/III", "内侧内嗅皮层 II／III 层", "Medial entorhinal cortex layers II/III", "海马旁系统"),
)
# Descriptions summarize the named source report rather than reinterpret annotations.
DATASETS = {
    "PNAS-V1": ("mouse", "Sooter et al. · Mouse V1", "02_mouse_replication", "已排序皮层尖峰", "作者发布 WAKE 和 NREM 区间", "沿用作者 WAKE 标注。", "动物与会话的对应关系未核实。仅按会话描述，REM 无已复现结果。"),
    "000552": ("mouse", "Preconfigured dynamics in the hippocampus are guided by embryonic birthdate and rate of neurogenesis", "13_sleep_within_region", "海马胞外尖峰与 LFP", "NWB processing/behavior/SleepStates", "作者家笼 Awake。排除任务及光刺激前后 60 秒。", "36 个合格会话中有 33 个提供可用 REM。"),
    "000056": ("mouse", "Peyrache et al. · Anterodorsal thalamus and postsubiculum", "24_phase10_multiregion_results", "多单元胞外记录、LFP、行为视频", "新版 NWB Buzcode 分期。论文分期另列。", "主分析采用休息场景的作者 Wake。排除开放场和干预。视频低移动 Wake 另列。", "部分 ADN 和 PoSub 来自同一动物或同时记录。两种分期使用同一批记录。"),
    "000939": ("mouse", "Duszkiewicz et al. · Postsubiculum", "24_phase10_multiregion_results", "PoSub 胞外尖峰和 LFP", "NWB 作者睡眠状态表", "确认家笼时间并避开干预及其前后 60 秒。", "31 个会话中 2 个存在时钟冲突，已排除。4 个合格会话缺可用 REM。"),
    "zenodo:19462601": ("mouse", "High-density electrophysiological recordings in midbrain and pons during sleep", "24_phase10_multiregion_results", "Neuropixels、EEG／EMG", "remidx 作者评分，1 ms 尖峰与 EEG／EMG 时钟独立核对", "WT 无刺激对照的作者 Wake。", "不同核团共享部分动物。部分核团仅一只。P、MB 是来源大区标签。"),
    "000041": ("rat", "Network Homeostasis and State Dynamics of Neocortical Sleep", "13_sleep_within_region", "额叶胞外尖峰与 LFP", "NWB 作者 SleepStates", "作者家笼 Wake。", "汇总包含多个额叶子区。子区与队列共享动物。"),
    "000061": ("rat", "Reactivations of emotional memory in the hippocampus–amygdala system during sleep", "14_sleep_cross_region", "CA1／BLA 同步胞外记录", "NWB 作者 SleepStates，按 epoch 区分训练前后", "保留 PRE／POST 自发时段。排除跑道、训练和刺激。", "CA1 与 BLA、PRE 与 POST 共享动物。"),
    "hc-23": ("rat", "Entorhinal cortex superficial-layer recordings", "19_phase9_sleep_results", "MEC 浅层胞外尖峰与 LFP", "作者发布的 Overnight 状态标注", "来源缺 Wake。", "仅支持 REM−NREM。"),
    "hc-24": ("rat", "Simultaneous recordings in rat hippocampus, prefrontal cortex, and thalamus", "21_hc24_nrem", "胞外尖峰与 600 Hz LFP", "DerivedData/SleepTimes.mat:nremSleep", "来源未发布可用 Wake／REM。", "仅展示核验过的 NREM 估计。丘脑单元数量未通过门槛。"),
}
DANDI_PAPERS = {
    "000041": "https://doi.org/10.1016/j.neuron.2016.03.036",
    "000061": "https://doi.org/10.1038/nn.4637",
    "000552": "https://doi.org/10.1038/s41593-022-01138-x",
}
CONTEXTS = {
    "HOME": "家笼", "REST": "休息 · 参考分期", "REST_PAPER": "休息 · 原论文评分",
    "REST_QUIET_EXPLORATORY": "视频低移动 Wake · 探索", "PRE": "训练前",
    "POST": "训练后", "POST_RUN": "跑后睡眠", "OVERNIGHT": "过夜状态", "PUBLISHED": "作者公开会话",
}
