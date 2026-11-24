import { EVIDENCE_LABELS } from "./components/PreliminaryDetails";
import { StrictMode, useEffect, useState } from "react";
import { createRoot } from "react-dom/client";
import OverviewPlot from "./components/OverviewPlot";
import Brain, { Chip } from "./components/Brain";
import { EffectPlot, StatePlot, Unavailable } from "./components/Charts";
import { Provenance, ReplicateDetail, Robustness } from "./components/Details";
import {
  type Anatomy, type Atlas, type Cohort, type Contrast, type Region, type Replicate, type Species, type View,
  DELTA, NO_PAIR, PAIRS, SPECIES, STATES, effectColor, json, paired, signed,
} from "./data";
import "./style.css";

const REPO = "https://github.com/Lifescience-Lyyyyyy/Criticality-in-Sleep";
const VIEWS: Record<View, string> = { NW: "Wake / NREM", RN: "NREM / REM", RW: "Wake / REM", all: "三态" };

function RegionIndex({ regions, region, contrast, select }: {
  regions: Region[]; region: Region; contrast: Contrast; select: (id: string) => void;
}) {
  const [query, setQuery] = useState("");
  const list = regions.filter(r => [r.id, r.name, r.english, r.system, r.allen].some(x => x?.toLowerCase().includes(query.trim().toLowerCase())));
  return <section className="card region-index" aria-labelledby="region-list-heading">
    <div className="index-heading">
      <h2 id="region-list-heading">脑区索引</h2>
      <label className="search">
        <svg viewBox="0 0 24 24" width="16" height="16" aria-hidden="true"><circle cx="10" cy="10" r="6.5" /><path d="m15 15 5 5" /></svg>
        <input aria-label="搜索脑区" type="search" placeholder="搜索脑区、缩写或系统" value={query} onChange={e => setQuery(e.target.value)} />
      </label>
    </div>
    <div className="region-grid">
      {list.map(r => {
        const c = r.cohorts[0];
        const status = c == null ? "无估计" : c.animals == null || c.animals < 5 ? "有限证据" : "已分析";
        const statusKey = c == null ? "empty" : c.animals == null || c.animals < 5 ? "limited" : "ready";
        return <button type="button" key={r.id} className="region-item" aria-label={r.name + " " + r.id + "，" + status}
          aria-pressed={r.id === region.id} onClick={() => select(r.id)}>
          <span className="index-dot" style={{background: effectColor(c?.effects.reference[contrast].mean)}} />
          <span className="region-label"><b>{r.id}</b><span>{r.name}</span></span>
          <small data-status={statusKey}>{status}</small>
        </button>;
      })}
      {list.length === 0 && <p className="muted">没有匹配的脑区。</p>}
    </div>
  </section>;
}

function RegionResult({ region, atlas, cohort, isLoading, failure, view, setView, selectDataset, selectContext }: {
  region: Region; atlas: Atlas; cohort?: Cohort; isLoading: boolean; failure: boolean; view: View;
  setView: (v: View) => void; selectDataset: (id: string) => void; selectContext: (id: string) => void;
}) {
  const [selected, choose] = useState<string | null>(null);
  const evidence = cohort ? atlas.preliminary?.evidence[cohort.id] : undefined;
  const contrast = view === "all" ? "NW" : view;
  const states = view === "all" ? STATES : PAIRS[view];
  const rows = cohort ? paired(cohort.records, "reference", states) : [];
  const summary = cohort?.effects.reference[contrast];
  const groups = [...new Set(region.cohorts.map(c => c.dataset))];
  const currentDataset = cohort?.dataset ?? groups[0];
  const choices = region.cohorts.filter(c => c.dataset === currentDataset);
  const groupLabel = cohort?.animals == null ? "会话" : "动物";
  const connect = rows.length > 0;
  const present = cohort?.states ?? [];
  let distribution = rows;
  if (!connect && cohort && view === "all") {
    const byId = new Map<string, Replicate>();
    for (const state of present) for (const x of paired(cohort.records, "reference", [state])) {
      const base = byId.get(x.id);
      byId.set(x.id, base ? {
        ...base, values: {...base.values, ...x.values},
        sessions: [...new Map([...base.sessions, ...x.sessions].map(s => [s.id, s])).values()],
        count: new Set([...base.sessions, ...x.sessions].map(s => s.id)).size,
      } : x);
    }
    distribution = [...byId.values()];
  }
  const focused = distribution.find(r => r.id === selected);
  return <section className="card result-card" id="region-result" aria-labelledby="selected-heading" aria-busy={isLoading}>
    <header>
      <div className="result-context"><span className="abbreviation">{region.id}</span><span>{SPECIES[region.species]} · {region.system}</span></div>
      <h2 id="selected-heading">{region.name}</h2>
      <p className="region-en">{region.english}</p>
      {region.about && <p className="region-about">{region.about}</p>}
    </header>
    {groups.length > 0 && <div className="selectors">
      <label><span>数据来源</span><select aria-label="数据来源" value={currentDataset} onChange={e => selectDataset(e.target.value)}>
        {groups.map(id => <option key={id} value={id}>{id === "PNAS-V1" ? "PNAS · Mouse V1" : id === "zenodo:19462601" ? "Zenodo 19462601 · 脑干 WT" : id + " · " + atlas.datasets[id].title.split(" · ")[0]}</option>)}
      </select></label>
      {choices.length > 1 && <label><span>记录场景</span><select aria-label="记录场景" value={cohort?.id ?? choices[0].id} onChange={e => selectContext(e.target.value)}>
        {choices.map(c => <option key={c.id} value={c.id}>{c.contextLabel}</option>)}
      </select></label>}
    </div>}
    {evidence && <div className="evidence-box">
      <p className="evidence-tag" data-category={evidence.category}>{EVIDENCE_LABELS[evidence.category]}</p>
      <p className="muted evidence-note">{evidence.reason}</p>
    </div>}
    {cohort && <dl className="region-metrics">
      <div><dt>数据来源</dt><dd>{groups.length}</dd></div>
      <div><dt>独立动物</dt><dd>{cohort.animals ?? <small>未核实</small>}</dd></div>
      <div><dt>记录会话</dt><dd>{cohort.sessions}</dd></div>
      <div><dt>可用状态</dt><dd className="state-chips">{present.map(s => <Chip key={s}>{s}</Chip>)}</dd></div>
    </dl>}
    <div className="comparison-tabs segmented" aria-label="状态对比">
      {(Object.keys(VIEWS) as View[]).map(v => <button type="button" aria-pressed={v === view} key={v}
        onClick={() => {choose(null); setView(v);}}>{VIEWS[v]}</button>)}
    </div>
    {failure ? <Unavailable detail="结果文件加载失败。请刷新页面重试。" />
      : isLoading ? <p className="blank-chart">正在载入会话数据…</p>
      : !cohort ? <Unavailable detail="现有记录未获得可展示的状态估计。" />
      : view !== "all" && (!summary?.n || !rows.length) ? <Unavailable detail={"有效状态：" + present.join("、") + "。"} />
      : <>
        <div className="chart-intro">
          <div><h3>{view === "all" ? "三态 d₂" : DELTA[contrast]}</h3>
            <p>{connect ? "同场记录配对" : "各可用状态的描述分布"} · {distribution.length} {groupLabel}{connect ? " · " + rows.reduce((total, row) => total + row.count, 0) + " 场" : ""}</p>
          </div>
          {view !== "all" && <div className="big-effect" data-sign={summary?.mean == null ? "none" : summary.mean > 0 ? "pos" : summary.mean < 0 ? "neg" : "zero"}><strong>{signed(summary?.mean)}</strong><span>bits/s</span></div>}
        </div>
        {view === "all" && !connect && <p className="missing-states">{present.length < 3 ? "缺少可用状态：" + STATES.filter(s => !present.includes(s)).join("、") : "缺少同会话三态支持。"} {NO_PAIR}。</p>}
        <StatePlot rows={distribution} states={states} selected={selected} choose={choose} connect={connect} />
        {view !== "all" && summary && <div className="effects-band">
          <EffectPlot rows={rows} selected={selected} choose={choose} summary={summary} />
          <p>{summary.positive} 正向 · {summary.negative} 负向
            {summary.n - summary.positive - summary.negative > 0 && " · " + (summary.n - summary.positive - summary.negative) + " 零差"}
          </p>
        </div>}
        <label className="replicate-selector"><span>查看{groupLabel}与各次记录</span>
          <select value={selected ?? ""} onChange={e => choose(e.target.value || null)} aria-label={"选择" + groupLabel}>
            <option value="">选择{groupLabel}</option>
            {distribution.map(r => <option key={r.id} value={r.id}>{r.id} · {r.count} 场{r.effect == null ? "" : " · Δ " + signed(r.effect)}</option>)}
          </select>
        </label>
        {focused && <ReplicateDetail key={focused.id} chosen={focused} cohort={cohort} atlas={atlas} />}
      </>}
  </section>;
}

function Overview({ regions, region, contrast, select }: {
  regions: Region[]; region: Region; contrast: Contrast; select: (id: string) => void;
}) {
  return <section className="card lower-card" id="coverage">
    <header><h2>跨脑区覆盖</h2><p className="muted">每行显示默认数据源。共享动物涉及多个脑区时分别呈现。</p></header>
    <OverviewPlot regions={regions} active={region.id} select={id => { select(id); document.getElementById("atlas")?.scrollIntoView(); }} />
    <div className="table-scroll"><table className="data-table">
      <thead><tr><th>脑区</th><th>数据源</th><th className="num">动物</th><th className="num">会话</th><th>状态</th><th className="num">{DELTA[contrast]}</th></tr></thead>
      <tbody>{regions.map(r => {
        const d = r.cohorts[0], result = d?.effects.reference[contrast];
        const sign = result?.mean == null ? "none" : result.mean > 0 ? "pos" : result.mean < 0 ? "neg" : "zero";
        return <tr key={r.id} aria-selected={r.id === region.id}><th scope="row"><a href="#atlas" onClick={() => select(r.id)}>{r.id} · {r.name}</a></th>
          <td>{d?.dataset ?? "—"}</td><td className="num">{d?.animals ?? (d ? "未知" : "—")}</td><td className="num">{d?.sessions ?? "—"}</td>
          <td>{d?.states.join(" / ") ?? "—"}</td><td className="num effect-cell" data-sign={sign}>{result?.mean != null ? signed(result.mean) + " bits/s" : "—"}</td>
        </tr>;
      })}</tbody>
    </table></div>
  </section>;
}

function Method() {
  return <section className="method card lower-card" id="method">
    <h2>如何读这张图谱</h2>
    <div className="method-grid">
      <div><h3>模型距离</h3><p>d₂ 是估计的动力学模型到指定临界模型流形的信息论距离。值越低，越接近固定 β=2 的参考流形。AR 系数按会话与状态分别估计。β=2 是参考模型约束，生理相变需要独立证据。</p></div>
      <div><h3>先配对，再按动物汇总</h3><p>同会话按状态配对，然后等权汇总同一动物的会话。V1 只有可核实的会话身份，图谱保留其会话口径。三态连线只使用同场具备全部三态的记录。</p></div>
      <div><h3>默认数据源</h3><p>先选主要记录场景，再按独立动物数、会话数和来源标识排序。训练前后分别保留。多来源只在同一个来源内求均值。地图着色为参考 AR(5) 的有符号平均差。</p></div>
    </div>
    <div className="links">
      <a href={REPO + "/blob/main/reports/27_sleep_criticality_atlas_data.md"}>图谱数据说明</a>
      <a href={REPO + "/blob/main/reports/30_preliminary_research_summary.md"}>初步研究汇总</a>
      <a href={REPO + "/blob/main/reports/29_preliminary_validation.md"}>本轮验证</a>
      <a href={REPO}>项目仓库</a>
    </div>
  </section>;
}

function App() {
  const [atlas, setAtlas] = useState<Atlas | null>(null);
  const [anatomy, setAnatomy] = useState<Anatomy | null>(null);
  const [anatomyError, setAnatomyError] = useState(false);
  const [regionLists, setLists] = useState<Record<Species, Region[]> | null>(null);
  const [error, setError] = useState(false);
  const [species, setSpecies] = useState<Species>("mouse");
  const [regionId, setRegion] = useState("ADN");
  const [choice, setChoice] = useState<string | null>(null);
  const [view, setView] = useState<View>("NW");
  const [geoContrast, setGeoContrast] = useState<Contrast>("NW");
  const [loaded, setLoaded] = useState<{path: string; cohorts: Cohort[]} | null>(null);
  const [failed, setFailed] = useState<string | null>(null);
  useEffect(() => {
    const controller = new AbortController(), signal = controller.signal;
    Promise.all([json<Atlas>("atlas.json", signal), json<Region[]>("mouse.json", signal), json<Region[]>("rat.json", signal)])
      .then(([index, mouse, rat]) => { if (!signal.aborted) {setAtlas(index); setLists({mouse, rat});} })
      .catch(() => { if (!signal.aborted) setError(true); });
    json<Anatomy>("anatomy.json", signal).then(setAnatomy).catch(() => { if (!signal.aborted) setAnatomyError(true); });
    return () => controller.abort();
  }, []);
  const regions = regionLists?.[species] ?? [];
  const region = regions.find(r => r.id === regionId) ?? regions[0];
  const current = region?.cohorts.find(c => c.id === choice) ?? region?.cohorts[0];
  const path = region?.path;
  useEffect(() => {
    if (!path) return;
    const controller = new AbortController(), signal = controller.signal;
    json<{cohorts: Cohort[]}>(path, signal).then(r => {
      if (!signal.aborted) {setFailed(null); setLoaded({path, cohorts: r.cohorts});}
    }).catch(() => { if (!signal.aborted) setFailed(path); });
    return () => controller.abort();
  }, [path]);
  const cohort = path && loaded?.path === path ? loaded.cohorts.find(c => c.id === current?.id) : undefined;
  const unavailable = Boolean(path && failed === path);
  const isLoading = Boolean(path && !unavailable && loaded?.path !== path);
  const select = (id: string) => { setRegion(id); setChoice(null); };
  const switchSpecies = (s: Species) => {
    if (s === species) return;
    setSpecies(s);
    select(s === "mouse" ? "ADN" : "Frontal");
  };
  const updateView = (v: View) => { setView(v); if (v !== "all") setGeoContrast(v); };
  return <>
    <header className="topbar">
      <a href="#atlas" className="brand">
        <svg viewBox="0 0 36 36" width="34" height="34" aria-hidden="true">
          <rect width="36" height="36" rx="9" /><path d="M6 19h6l3-8 4 15 3-8h8" />
        </svg>
        <span>Sleep Criticality Atlas</span>
      </a>
      <nav aria-label="页面导航"><a href="#atlas">脑区图谱</a><a href="#coverage">覆盖概览</a><a href="#method">方法与来源</a></nav>
    </header>
    <main>
      <div className="intro" id="atlas">
        <div className="intro-copy">
          <span className="eyebrow">Allen CCFv3 · AR(5) · β = 2</span>
          <h1>睡眠临界性脑图谱</h1>
          <p>Exploring Neural Proximity to Criticality Across Sleep States</p>
        </div>
        <div className="segmented species-switch" aria-label="物种">
          {(["mouse", "rat"] as const).map(s => <button type="button" key={s} onClick={() => switchSpecies(s)}
            aria-pressed={species === s}>{SPECIES[s]}</button>)}
        </div>
      </div>
      {error ? <Unavailable detail="图谱数据加载失败。请刷新页面重试。" /> : !atlas || !region
        ? <p className="blank-chart">正在载入图谱…</p>
        : <>
          <div className="hero-grid">
            <div className="hero-left">
              {species === "mouse"
                ? <Brain data={anatomy} loadingError={anatomyError} regions={regions} region={region} current={current}
                  contrast={geoContrast} select={select} />
                : <section className="card rat-heading"><h2>大鼠数据探索</h2><p>按脑区查看大鼠数据源。鼠种、动物与记录身份在每个队列中保留。</p></section>}
              <RegionIndex key={species} regions={regions} region={region} contrast={geoContrast} select={id => {
                select(id);
                if (window.matchMedia("(max-width: 790px)").matches) document.getElementById("region-result")?.scrollIntoView();
              }} />
            </div>
            <RegionResult key={current?.id ?? region.id} region={region} atlas={atlas} cohort={cohort}
              isLoading={isLoading} failure={unavailable} view={view} setView={updateView}
              selectDataset={id => setChoice(region.cohorts.find(c => c.dataset === id)?.id ?? null)}
              selectContext={setChoice} />
          </div>
          {cohort && <div className="lower-grid">
            <Robustness key={cohort.id} cohort={cohort} atlas={atlas} region={region} contrast={geoContrast} />
            <Provenance dataset={atlas.datasets[cohort.dataset]} cohort={cohort} atlas={atlas}
              states={view === "all" ? STATES : PAIRS[view]} />
          </div>}
          <Overview regions={regions} region={region} contrast={geoContrast} select={select} />
          <Method />
        </>}
    </main>
    <footer><span>Sleep Criticality Atlas</span><a href={REPO}>GitHub · 复现代码与数据说明</a></footer>
  </>;
}

createRoot(document.getElementById("root")!).render(<StrictMode><App /></StrictMode>);
