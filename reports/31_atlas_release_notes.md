# 睡眠临界性脑图谱 · Phase 11.1 发布说明

图谱以小鼠 Allen CCFv3 为主视图，大鼠单独呈现。三态模式也显示脑图颜色所代表的差值。结果卡片标出证据类别、独立动物数和当前对比的样本数。V1 的计量单位仍为会话，缺少状态时显示缺失。

000939 PoSub 将原等时、新成对等时与固定三态对照并列。样本数随 NW、RN、RW 切换，参考数值继续单独显示。排序敏感性面板可选 Mouse28-140313 PoSub、DL108 的 LDT 以及本轮合格补充案例，直接比较各自的全单元参考值与所选子集。确认的数据排除和模型拟合告警分别标识。数据版本、发布 DOI、分析提交和 MATLAB 运行保留可追溯链接。

## 构建与本机打开

在 WSL 的现有仓库根目录内执行。用 Node.js 22 以及已有 Python 环境。

    python3 scripts/atlas/build_atlas_data.py
    python3 scripts/atlas/validate_atlas_data.py
    cd website
    npm ci
    npx playwright install chromium
    npm run build
    npm test
    npm run preview -- --host 127.0.0.1

按终端给出的本机地址访问。将 dist 目录放到静态 HTTP 服务上也可使用。构建使用相对资源路径，允许从 GitHub 项目子路径载入。所有绘图数据均随静态包存放，点击链接查看原始档案才与其网站通信。

## GitHub Pages

原运行 44779662303 的静态构建通过，配置 Pages 的步骤返回 404。因此要确认公网是否已开放，先核查仓库 [Settings → Pages](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/settings/pages)。若尚未配置，请仓库管理员在 Build and deployment 的 Source 选择 **GitHub Actions**。然后打开 [Build and deploy Sleep Criticality Atlas](https://github.com/Lifescience-Lyyyyyy/Criticality-Evolution-in-Sleep/actions/workflows/atlas-pages.yml)，从 main 运行现有工作流。发布时使用 Deploy 日志给出的地址及部署提交核验版本。[GitHub 官方配置步骤](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site)

[本轮科学审计](29_preliminary_validation.md) · [初步科学汇总](30_preliminary_research_summary.md) · [解剖来源及 Phase 11 口径](27_sleep_criticality_atlas_data.md)
