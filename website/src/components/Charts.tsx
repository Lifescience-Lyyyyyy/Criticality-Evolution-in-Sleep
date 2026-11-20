import { amount, signed, STATE_COLOR, type Replicate, type State, type Stats, average, NO_PAIR } from "../data";

type Props = { rows: Replicate[]; states: State[]; selected: string | null; choose: (id: string) => void };
const tick = (v: number, max: number) => 198 - (v / max) * 155;

export function StatePlot({ rows, states, selected, choose, connect }: Props & { connect: boolean }) {
  const points = rows.flatMap(r => states.flatMap(s => r.values[s] ?? []));
  const max = Math.max(.5, Math.ceil(Math.max(0, ...points) * 2) / 2);
  const x = (s: State) => 104 + states.indexOf(s) * (states.length === 3 ? 126 : 212);
  const ticks = Array.from({length: 5}, (_, i) => max / 4 * i);
  return <svg role="group" aria-label={connect ? "同会话配对状态 d₂" : "按状态展示 d₂"} className="measure-chart" viewBox="0 0 460 246">
    <text x="35" y="22" className="plot-axis-title">d₂ (bits/s)</text>
    {ticks.map(v => <g key={v}>
      <path d={"M 47 " + tick(v, max) + " H 414"} className="plot-grid" />
      <text x="39" y={tick(v, max) + 4} textAnchor="end" className="plot-axis">{amount(v, 2).replace(/0$/, "")}</text>
    </g>)}
    <path d="M47 43 V198 H414" className="axis-line" />
    {states.map(s => <text key={s} x={x(s)} y="225" textAnchor="middle" className="plot-state">{s}</text>)}
    {rows.map((r, n) => {
      const offset = Math.sin(n * 2.24) * (connect ? 3 : 21);
      const pairs = states.filter(s => r.values[s] != null).map(s => ({s, x: x(s) + offset, y: tick(r.values[s]!, max)}));
      const chosen = r.id === selected;
      return <g key={r.id} aria-label={r.id + " " + pairs.map(p => p.s + " " + String(r.values[p.s])).join(", ")}
        data-selected={chosen} className="chart-replicate" tabIndex={0} role="button"
        onClick={() => choose(r.id)}
        onKeyDown={e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); choose(r.id); } }}>
        <title>{r.id} · {pairs.map(p => p.s + " " + String(r.values[p.s])).join(" / ")} · {r.count} 场</title>
        {connect && pairs.length > 1 && <polyline points={pairs.map(p => p.x + "," + p.y).join(" ")} />}
        {pairs.map(p => <circle key={p.s} cx={p.x} cy={p.y} r={chosen ? 5 : 3.4} fill={STATE_COLOR[p.s]} />)}
      </g>;
    })}
    {states.map(s => {
      const values = rows.flatMap(r => r.values[s] ?? []);
      if (!values.length) return null;
      return <path key={s} d={"M " + (x(s) - 11) + " " + tick(average(values), max) + " h22"} className="mean-dash" />;
    })}
  </svg>;
}

export function EffectPlot({ rows, selected, choose, summary }: Omit<Props, "states"> & { summary: Stats }) {
  if (summary.mean == null) return null;
  const limit = Math.max(.25, Math.ceil(Math.max(...rows.map(r => Math.abs(r.effect!))) * 4) / 4);
  const x = (v: number) => 45 + ((v + limit) / (2 * limit)) * 370;
  return <svg viewBox="0 0 460 79" className="effect-chart" aria-label="逐个体配对差值（bits/s）">
    <text x="7" y="15" className="plot-axis-title">Δd₂</text>
    <line x1={x(-limit)} x2={x(limit)} y1="34" y2="34" className="axis-line" />
    <line x1={x(0)} x2={x(0)} y1="14" y2="51" className="zero-line" />
    {[-limit, 0, limit].map(v => <text key={v} x={x(v)} y="66" textAnchor="middle" className="plot-axis">{signed(v, 2)}</text>)}
    {rows.map((r, n) => <circle key={r.id} cx={x(r.effect!)} cy={34 + Math.sin(n * 2.24) * 10}
      r={selected === r.id ? 5 : 3.4} className={r.effect! < 0 ? "negative-point" : "positive-point"}
      stroke={selected === r.id ? "#153b47" : "white"} strokeWidth={selected === r.id ? 1.8 : .6}
      role="button" tabIndex={0} aria-label={r.id + " Δd₂ " + signed(r.effect, 7)}
      onClick={() => choose(r.id)} onKeyDown={e => { if (e.key === "Enter") choose(r.id); }}>
      <title>{r.id} · Δd₂ {r.effect} bits/s</title>
    </circle>)}
    <path d={"M " + x(summary.mean) + " 18 V 48"} className="mean-dash" pointerEvents="none" />
  </svg>;
}

export function Unavailable({ detail }: { detail?: string }) {
  return <div className="blank-chart"><strong>{NO_PAIR}</strong>{detail && <p>{detail}</p>}</div>;
}
