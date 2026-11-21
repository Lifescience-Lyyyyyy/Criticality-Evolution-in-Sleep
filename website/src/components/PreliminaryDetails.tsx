import { useState } from "react";
import { DELTA, signed, type Contrast, type PreliminaryAtlas, type SortingExample, type Stats } from "../data";

export const EVIDENCE_LABELS: Record<string, string> = {
  ROBUST_PRELIMINARY: "较稳健的初步证据",
  SENSITIVE: "对分析规格敏感",
  EXPLORATORY: "探索性",
  UNRESOLVED: "尚待核实",
  INVALIDATED: "当前比较无效",
};

function SampleLine({ title, stats }: {title: string; stats: Stats}) {
  return <div className="prelim-stat">
    <span>{title}</span>
    <strong>{signed(stats.mean)} <small>bits/s</small></strong>
    <span>{stats.n} 只 · {stats.sessions} 场 · {stats.positive} 正 / {stats.negative} 负</span>
  </div>;
}

export function PairedDuration({ data, contrast }: { data: PreliminaryAtlas["pairwise"]; contrast: Contrast }) {
  const row = data.contrasts[contrast];
  return <section className="prelim-panel" aria-label="成对等时对照">
    <h3>成对等时对照</h3>
    <p className="muted">{DELTA[contrast]} · 每态 3 × 60 秒 · AR(5)</p>
    <div className="prelim-grid">
      <SampleLine title="原等时规则" stats={row.legacy} />
      <SampleLine title="新成对规则" stats={row.corrected} />
      <SampleLine title="三态固定等时" stats={row.three_state_fixed} />
    </div>
    <p className="muted">共同 {row.matched_animals} 只的均差变化 {signed(row.corrected_minus_legacy_same_animals)} bits/s。
      {row.corrected.bootstrap95 && (" 新成对规则的动物 bootstrap 95% 区间 [" + row.corrected.bootstrap95.map(x => signed(x)).join(", ") + "]。")}
    </p>
    <div className="links"><a href={data.provenance.action}>MATLAB 运行</a><a href={data.source}>逐动物结果和排除依据</a></div>
  </section>;
}

const SORT_LABELS: Record<string,string> = {
  lratio_q75: "低 L-ratio q75", lratio_q50: "低 L-ratio q50",
  contamination_q75: "低污染率 q75", contamination_q50: "低污染率 q50",
  isi_q50: "低短 ISI 比例 q50",
};

export function SortingCompare({ cases, contrast }: {cases: SortingExample[]; contrast: Contrast}) {
  const [session, setSession] = useState(cases[0].session);
  const [variant, setVariant] = useState("");
  const item = cases.find(c => c.session === session) ?? cases[0];
  const chosen = item.variants.find(v => v.id === variant) ?? item.variants[item.variants.length-1];
  const a = item.reference.effects[contrast];
  const b = chosen.effects[contrast];
  return <section className="prelim-panel" aria-label="排序敏感性对比">
    <h3>排序敏感性对比</h3>
    <p className="muted">同一会话 · AR(5) · 原分期与时间箱</p>
    <div className="prelim-grid">
      <label>会话<select aria-label="会话" value={item.session} onChange={e => {setSession(e.target.value); setVariant("");}}>
        {cases.map(c => <option value={c.session} key={c.session}>{c.session}</option>)}
      </select></label>
      <label>排序子集<select aria-label="排序子集" value={chosen.id} onChange={e => setVariant(e.target.value)}>
        {item.variants.map(v => <option value={v.id} key={v.id}>{SORT_LABELS[v.id] ?? v.id}</option>)}
      </select></label>
    </div>
    <dl className="prelim-sort-values" aria-live="polite">
      <div><dt>参考全组 · {item.reference.units} 单元</dt><dd>{signed(a, 5)} bits/s</dd></div>
      <div><dt>所选子集 · {chosen.units} 单元</dt><dd>{signed(b, 5)} bits/s</dd></div>
    </dl>
    <p className="muted">{DELTA[contrast]} · {a == null || b == null ? "缺少可配对状态" :
      "差值变化 " + signed(b-a,5) + " bits/s · " + (Math.sign(a) === Math.sign(b) ? "方向保持" : "符号改变")}</p>
    <div className="links"><a href={item.run}>MATLAB 运行</a><a href={item.source}>逐状态与单元来源</a></div>
  </section>;
}
