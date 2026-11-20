# 睡眠临界性脑图谱 · 数据与复现

图谱默认显示小鼠，另一入口仅查询大鼠。源结果基线为 a834a85ce765dfb2fb77f267ffd3410c64a1f339。Phase 11 补充核验使用独立目录及拟合收据。历史数据、拟合与报告可以继续按各自提交、工具箱版本和 GitHub Actions 运行复现。

## 导出内容

脚本读取仓库的原模型拟合 JSON、动物与会话对比 CSV、独立重抽样对照汇总及数据源清单。生成文件由各物种的脑区索引、每脑区的会话数据、来源与运行索引、解剖切片组成。每项有效拟合都能追溯至物种、脑区、数据集、记录场景、动物或会话、公开资产 ID、状态和规格。数据表保留公开资产路径和校验值、原拟合文件 SHA-256、Actions 运行、分析提交、MATLAB 版本和估计器提交。时钟区间与确切单元 ID 通过原拟合文件追溯，网页只加载交互需要的字段。

默认数据来源在主要记录场景中优先选独立动物最多的来源，再比较会话数和稳定标识。训练前优先显示，训练后可另选。PoSub 首先显示 000939，000056 可手动切换。两种分期和低移动 Wake 在 000056 的记录场景中分别保留。同一动物的多个会话先在会话内配对，再按动物等权平均。V1 的 13 场以会话为单位，独立动物数未知。000061 两个脑区及训练前后共享动物。000056、脑干分区以及大鼠额叶总组与各子区也共享动物。

各状态采用各自通过门槛的有效会话。配对图只连同会话都具有所选状态的值。三态图有同会话三态时使用完整三态会话，否则仅展示缺失状态与可用状态的描述分布。地区颜色使用参考 AR(5) 的平均有符号效应，范围为 ±1.5 bits/s，越界颜色饱和。脑区索引以少于五只已核实动物或未知动物数描述有限证据。这是浏览提示，原拟合和动物均保留。RSC 保留解剖入口，现有合格结果缺状态对比。

## 解剖来源

使用 Allen Institute 成年小鼠 CCFv3 October 2017 的 25 μm 标注体和解剖层级，按 AP、SI、LR 阅读体素。界面以 50 μm 采样生成 SVG 轮廓，按选中区域的最大截面积自动选择矢状或冠状切片。矢状面在左半球及中线中选取，冠状面左至右显示，两者均为背侧向上。界面 AP、LR 数值以该 CCF 体素原点计，比例尺来自 50 μm 采样。冷暖填色为 d₂ 差值，实验电极尚未映射到这些体素。图形用作成年小鼠解剖参考。V1 映射 VISp，ADN 映射 AD，PoSub 映射 POST，RSC 映射 RSP，脑干细核团保留源数据缩写。详细 CCF 结构 ID 与两个源文件的 SHA-256 见 anatomy.json。图谱依 Allen Institute 的研究或非商业用途条款分发并在图上署名。[Allen 技术文档](https://brain-map.org/support/documentation/api-allen-brain-connectivity-atlas) · [Wang et al., 2020](https://doi.org/10.1016/j.cell.2020.04.007) · [使用条款](https://alleninstitute.org/legal/terms-of-use)

## 本地运行

在仓库根目录执行。需要 Python 3.11+ 和 Node.js 22.12+。

    python3 scripts/atlas/build_atlas_data.py
    python3 scripts/atlas/validate_atlas_data.py
    cd website
    npm ci
    npm run dev

Vite 会打印开发地址。验证生产产物可运行 `npm run build` 和 `npm run preview`。文件生成在 `website/dist`，可由任意静态 HTTP 服务器提供。浏览过程中无需访问 DANDI、Figshare、Zenodo、CRCNS 或 Drive。点击原始来源链接才会访问对应站点。

如需重新生成已提交的解剖轮廓，从 `scripts/atlas/build_anatomy.py` 顶部的固定来源下载两个文件，核验脚本中的校验值。用具备 NumPy、Matplotlib 和 ContourPy 的环境运行。

    python scripts/atlas/build_anatomy.py annotation_25.nrrd structure_graph_1.json

## 验证与部署

`validate_atlas_data.py` 对照已发布文件的哈希，重新生成全部 JSON 并逐字节比对，逐项匹配 Phase 10 的会话和动物效应，对照 Phase 8/9 已存动物结果与 AR(10) 差值，并核对 V1 与 hc-24 的特殊口径。浏览器测试阻断站点以外的请求，覆盖解剖区点击及键盘选择、物种和数据源切换、状态差值、REM 缺失、动物与会话层级、原文链接、AR 控制和移动布局。可在 `website/` 中安装 Chromium 后执行。

    npx playwright install chromium
    npm run build
    npm test

`atlas-pages.yml` 在 main 上校验数据、构建、测试、上传并发布静态产物。Vite 的相对资源路径兼容仓库子路径。在 GitHub 的 Settings → Pages 将 Source 设为 GitHub Actions，再从 Actions 运行 Build and deploy Sleep Criticality Atlas。各 Phase 1–10 的 MATLAB 入口继续使用原工作流。Phase 11 补充回放有独立工作流。[GitHub Pages 自定义工作流](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages)
