import { PairedDuration, SortingCompare } from "./PreliminaryDetails";
import { useState } from "react";
import {
  type Atlas, type Cohort, type Contrast, type Dataset, type Replicate, type Region, type State,
  amount, average, DELTA, PAIRS, signed, SPEC, SPECIES, STATES, WARNING,
} from "../data";

const SOURCE = "https://github.com/Lifescience-Lyyyyyy/Criticality-in-Sleep";
const AUDIT = SOURCE + "/blob/main/reports/26_targeted_validation.md";
const STATUS: Record<string, string> = {
  CONFIRMED_IMPACT: "已确认有影响", NO_OBSERVED_IMPACT: "抽样未观察到影响",
  INCONCLUSIVE: "证据不足", NOT_APPLICABLE: "不适用",
};

export function SourceLinks({ dataset }: { dataset: Dataset }) {
  return <div className="links">
    <a href={dataset.url}>原始数据</a>
    <a href={dataset.paper}>原始论文</a>
    <a href={dataset.report}>项目报告</a>
    <a href={dataset.manifest}>来源清单</a>
  </div>;
}

export function Provenance({ dataset, cohort, atlas, states }: {
  dataset: Dataset; cohort: Cohort; atlas: Atlas; states: State[];
}) {
  const matching = cohort.records.filter(s => states.every(state => s.fits.reference[state]));
  const scope = matching.length ? states : cohort.states;
  const reference = (matching.length ? matching : cohort.records).flatMap(s => Object.entries(s.fits.reference)
    .filter(([state]) => scope.includes(state as State)).map(([, f]) => f!));
  const runs = [...new Set(reference.map(f => f.provenance))];
  const scores = [...new Set(reference.flatMap(f => f.score ?? []))];
  const exposure = reference.map(f => f.seconds);
  const units = reference.map(f => f.units);
  return <section className="card lower-card" id="provenance">
    <header><h2>数据与来源</h2><p className="muted">{dataset.title}</p></header>
    <dl className="metadata-grid">
      <div><dt>数据集</dt><dd>{dataset.id} · {dataset.version}</dd></div>
      <div><dt>物种与队列</dt><dd>{SPECIES[dataset.species]} · {cohort.contextLabel}</dd></div>
      <div><dt>该数据源覆盖</dt><dd>{dataset.animals ?? "未知"} 只 · {dataset.sessions} 场</dd></div>
      <div><dt>记录</dt><dd>{dataset.recording}</dd></div>
      <div><dt>分期</dt><dd>{scores.join("；") || dataset.scoring}</dd></div>
      <div><dt>Wake 定义</dt><dd>{cohort.context === "REST_QUIET_EXPLORATORY"
        ? "双 LED 均可见且速度 ≤2 cm/s，限定在作者家笼 Wake。" : dataset.wake}</dd></div>
      <div><dt>每次拟合的有效时长</dt><dd>{exposure.length ? amount(Math.min(...exposure), 1) + "–" + amount(Math.max(...exposure), 1) + " 秒" : "—"}</dd></div>
      <div><dt>群体单元数</dt><dd>{units.length ? Math.min(...units) + "–" + Math.max(...units) : "—"}</dd></div>
      <div><dt>参考规格</dt><dd>40 ms · AR(5) · β=2 · Yule–Walker</dd></div>
      <div><dt>访问</dt><dd>{dataset.access}</dd></div>
    </dl>
    <p className="source-note">{dataset.limitations}</p>
    <SourceLinks dataset={dataset} />
    <div className="lineage">
      {runs.map(id => {
        const p = atlas.provenance[id];
        return <div key={id}><a href={p.run}>MATLAB 运行 {id}</a><a href={p.file}>逐拟合原始输出</a>
          <a href={SOURCE + "/commit/" + p.commit}>分析提交 <code>{p.commit.slice(0, 8)}</code></a></div>;
      })}
    </div>
  </section>;
}

export function ReplicateDetail({ chosen, cohort, atlas }: { chosen: Replicate; cohort: Cohort; atlas: Atlas }) {
  return <div className="session-drill">
    <h3>{cohort.animals == null ? "会话" : "动物"} {chosen.id} <span className="muted">· {chosen.count} 场符合当前状态条件</span></h3>
    {chosen.sessions.map(s => <details key={s.id} name="animal-session">
      <summary>{s.session}</summary>
      <p className="asset-caption">{s.asset}</p>
      <div className="fit-list">
        {STATES.flatMap(state => {
          const f = s.fits.reference[state];
          if (!f) return [];
          const provenance = atlas.provenance[f.provenance];
          return [<div key={state}>
            <b>{state} · {amount(f.value, 8)} bits/s</b>
            <span>{f.units} 单元 · {amount(f.seconds, 2)} s · AR 根 {amount(f.root, 4)}</span>
            {f.warnings.length > 0 && <span>{f.warnings.map(code => WARNING[code]).join("，")}</span>}
            <a href={provenance.file}>原拟合</a>
          </div>];
        })}
      </div>
    </details>)}
  </div>;
}

export function Robustness({ cohort, atlas, region, contrast }: { cohort: Cohort; atlas: Atlas; region: Region; contrast: Contrast }) {
  const [requested, setRequested] = useState("");
  const [a, b] = PAIRS[contrast];
  const descriptive = cohort.states.length === 1 ? cohort.states[0] : null;
  const refStates = descriptive ? [descriptive] : [a, b];
  const refs = cohort.records.filter(s => refStates.every(state => s.fits.reference[state]));
  const controls = cohort.controls.filter(c => c.contrast === contrast);
  const choices = [...new Set([
    ...["AR10", "equal_duration", "edge_trim", "strict_units"].filter(spec => refs.some(s => refStates.every(state => s.fits[spec]?.[state]))),
    ...controls.map(x => x.model),
  ])];
  const spec = choices.includes(requested) ? requested : choices[0];
  const grouped = new Map<string, {references: number[]; values: number[]}>();
  const feed = (key: string, reference: number, value: number) => {
    if (!grouped.has(key)) grouped.set(key, {references: [], values: []});
    const group = grouped.get(key)!; group.references.push(reference); group.values.push(value);
  };
  const repeat = controls.filter(c => c.model === spec);
  for (const r of repeat) feed(r.subject, r.reference, r.value);
  if (!repeat.length) for (const s of refs) {
    const ref = s.fits.reference, alternative = s.fits[spec];
    if (descriptive && alternative?.[descriptive]) feed(s.subject ?? s.session,
      ref[descriptive]!.value, alternative[descriptive]!.value);
    else if (!descriptive && alternative?.[a] && alternative[b]) feed(s.subject ?? s.session,
      ref[b]!.value - ref[a]!.value, alternative[b]!.value - alternative[a]!.value);
  }
  const results = [...grouped].sort(([x], [y]) => x.localeCompare(y)).map(([id, x]) =>
    ({ id, reference: average(x.references), value: average(x.values) }));
  const fits = refs.flatMap(s => refStates.map(state => s.fits.reference[state]!));
  const audit = atlas.audit.filter(r => r.dataset === cohort.dataset && (!r.region || r.region === region.id));
  const modelsWarn = fits.filter(f => f.warnings.length);
  const diagnosticsMissing = fits.filter(f => f.kurtosis == null || f.residual == null);
  const preliminary = atlas.preliminary;
  const sortCases = preliminary?.sorting.filter(c => c.cohort === cohort.id) ?? [];
  const sourceErrors = audit.filter(r => r.category === "exclusion");
  const warnings = [...new Set(modelsWarn.flatMap(f => f.warnings.map(w => WARNING[w])))];
  return <section className="card lower-card" id="quality">
    <header><h2>质量与稳健性</h2><p className="muted">{descriptive ? descriptive + " d₂ · 同记录对照" : DELTA[contrast] + " · 同记录对照"}</p></header>
    {choices.length > 0 ? <>
      <label className="select-caption" htmlFor="control-select">分析规格</label>
      <select id="control-select" value={spec} onChange={e => setRequested(e.target.value)}>
        {choices.map(c => <option key={c} value={c}>{SPEC[c]}</option>)}
      </select>
      <div className="method-comparison">
        <div><span>{repeat.length ? "配平时间窗的参考" : "同会话参考 AR(5)"}</span><strong>{signed(average(results.map(x => x.reference)))}</strong></div>
        <div><span>{SPEC[spec]}</span><strong>{signed(average(results.map(x => x.value)))}</strong></div>
      </div>
      <p className="muted">{descriptive
        ? results.length + " 个同记录的描述估计 · bits/s"
        : results.filter(x => Math.sign(x.reference) === Math.sign(x.value)).length + "/" + results.length + " 个配对单位保留方向 · bits/s"}
        {repeat.length > 0 && " · 对照采用各动物的重复中位差"}</p>
      <details className="qc-detail"><summary>逐个体参考与对照</summary>
        <table className="data-table"><thead><tr><th>个体</th><th>参考</th><th>{SPEC[spec]}</th></tr></thead>
          <tbody>{results.map(x => <tr key={x.id}><th scope="row">{x.id}</th><td>{signed(x.reference, 6)}</td><td>{signed(x.value, 6)}</td></tr>)}</tbody>
        </table>
        {repeat.length > 0 && <div className="links"><a href={repeat[0].run}>对照 Actions</a><a href={repeat[0].source}>对照输出</a></div>}
      </details>
    </> : <p className="muted">这个状态对比缺少可并排比较的对照。</p>}
    {preliminary?.pairwise.cohort === cohort.id && <PairedDuration data={preliminary.pairwise} contrast={contrast} />}
    {sortCases.length > 0 && <SortingCompare cases={sortCases} contrast={contrast} />}
    <div className="qc-notes">
      {sourceErrors.length > 0 && <div className="source-error">
        <strong>已核实的数据排除</strong>
        {sourceErrors.map(r => <p key={r.session}>{r.session} · {r.inference}</p>)}
      </div>}
      {fits.length > 0 && <div><b>模型诊断</b><p>{modelsWarn.length}/{fits.length} 项参考拟合有模型充分性提示。
        {warnings.join("，")}{warnings.length > 0 ? "。这些提示未用于排除。" : ""}
        {diagnosticsMissing.length > 0 && " 部分先前输出未保存完整残差诊断。"}</p></div>}
      <details className="qc-detail"><summary>针对性审计 {audit.length > 0 && "(" + audit.length + ")"}</summary>
        {audit.length > 0 ? <ul className="audit-list">{audit.map((r, i) => <li key={i}>
          <b>{STATUS[r.status]} · {r.concern}</b>
          <span>{r.session} · {r.inference}</span>
          {r.dataset === "PNAS-V1" && r.category === "model" && r.original?.NW != null && r.audited?.NW != null &&
            <span>同一会话的 NREM−Wake · AR(5) {signed(r.original.NW)} · AR(10) {signed(r.audited.NW)} bits/s</span>}
          {r.evidence && <a href={r.evidence}>核验依据</a>}
        </li>)}</ul> : <p>该脑区未列入针对性抽样。</p>}
        <a href={AUDIT}>Phase 11 完整验证报告</a>
      </details>
    </div>
  </section>;
}
