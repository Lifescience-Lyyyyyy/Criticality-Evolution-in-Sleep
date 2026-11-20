import { signed, type Region } from "../data";

export default function OverviewPlot({ regions, active, select }: {
  regions: Region[]; active: string; select: (id: string) => void;
}) {
  const points = regions.flatMap(region => {
    const cohort = region.cohorts[0];
    const nw = cohort?.effects.reference.NW, rn = cohort?.effects.reference.RN;
    return nw?.mean == null || rn?.mean == null ? [] : [{region, nw, rn, x: nw.mean, y: rn.mean}];
  });
  const lim = Math.max(1.5, Math.ceil(Math.max(...points.flatMap(p => [Math.abs(p.x), Math.abs(p.y)])) * 2) / 2);
  const x = (v: number) => 74 + ((v + lim) / (2 * lim)) * 512;
  const y = (v: number) => 28 + ((lim - v) / (2 * lim)) * 168;
  return <figure className="coverage-plot">
    <svg className="overview-chart" viewBox="0 0 615 246" aria-label="各区域 NREM 减 Wake 与 REM 减 NREM">
      <text x="75" y="16" className="plot-axis-title">REM − NREM (bits/s)</text>
      {[-lim, -lim/2, 0, lim/2, lim].map(v => <g key={v}>
        <line x1="74" x2="586" y1={y(v)} y2={y(v)} className={v === 0 ? "zero-line" : "plot-grid"} />
        <line y1="28" y2="196" x1={x(v)} x2={x(v)} className={v === 0 ? "zero-line" : "plot-grid"} />
        <text x="63" y={y(v) + 4} className="plot-axis" textAnchor="end">{signed(v, 2)}</text>
        <text x={x(v)} y="213" className="plot-axis" textAnchor="middle">{signed(v, 2)}</text>
      </g>)}
      <text x="586" y="237" className="plot-axis-title" textAnchor="end">NREM − Wake (bits/s)</text>
      {points.map(p => <g key={p.region.id} role="button" tabIndex={0} className="overview-point"
        aria-label={p.region.name + " " + p.region.id + " " + signed(p.x) + " " + signed(p.y)}
        aria-pressed={active === p.region.id} onClick={() => select(p.region.id)} onKeyDown={e => {
          if (e.key === "Enter" || e.key === " ") { e.preventDefault(); select(p.region.id); }
        }}>
        <title>{p.region.name} · NW {p.x}（n={p.nw.n}）· RN {p.y}（n={p.rn.n}）</title>
        <circle cx={x(p.x)} cy={y(p.y)} r={active === p.region.id ? 6.25 : 4.5} />
      </g>)}
    </svg>
    <figcaption>每点为一个脑区的默认来源。横纵坐标的配对样本可能不同。选择点可切换脑区。</figcaption>
  </figure>;
}
