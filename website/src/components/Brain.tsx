import { useState } from "react";
import { type Anatomy, type Plane, type Region, type CohortSummary, type Contrast, effectColor, signed } from "../data";

export function Chip({children}: {children: React.ReactNode}) {
  return <span className="chip">{children}</span>;
}

export default function Brain({ data, loadingError, regions, region, current, contrast, select }: {
  data: Anatomy | null; loadingError: boolean; regions: Region[]; region: Region;
  current?: CohortSummary; contrast: Contrast; select: (id: string) => void;
}) {
  const [plane, setPlane] = useState<Plane>("sagittal");
  const [hover, setHover] = useState<string | null>(null);
  const view = data?.views[plane];
  const position = view?.positions[region.id];
  const slice = position == null ? undefined : view?.slices[String(position)];
  const regionToShow = regions.find(r => r.id === hover) ?? region;
  const pointer = slice?.labels[regionToShow.id] ?? slice?.labels[region.id];
  const value = (r: Region) =>
    (r.id === region.id ? current : r.cohorts[0])?.effects.reference[contrast].mean;
  const orientation = plane === "sagittal" ? ["前", "后"] : ["左", "右"];
  const layers = [...regions].sort((a, b) => (slice?.labels[b.id]?.area ?? 0) - (slice?.labels[a.id]?.area ?? 0));
  return <section className="card brain-card" aria-labelledby="brain-heading">
    <header className="card-header">
      <div><h2 id="brain-heading">小鼠脑区</h2><div className="muted">Allen CCFv3 · 解剖参考切片</div></div>
      <div className="segmented" aria-label="解剖视图">{(["sagittal", "coronal"] as const).map(p =>
        <button type="button" aria-pressed={plane === p} key={p}
          onClick={() => { setPlane(p); setHover(null); }}>{p === "sagittal" ? "矢状面" : "冠状面"}</button>)}</div>
    </header>
    <div className="brain-frame">
      <div className="anatomy-position">
        <span>{plane === "sagittal" ? "SAGITTAL" : "CORONAL"}</span>
        {position !== undefined && data ? <span>{plane === "sagittal" ? "LR" : "AP"} {(position * data.voxelMm).toFixed(2)} mm</span> : null}
      </div>
      {slice ? <svg className="brain-image" viewBox={[-12, -9, slice.w + 24, slice.h + 18].join(" ")} aria-label="交互小鼠脑切片">
        <path d={slice.outline} fill="#e6ebed" stroke="#8397a4" strokeWidth=".7" fillRule="evenodd" />
        <g fillRule="evenodd" stroke="#a8b6bd" strokeWidth=".32" aria-hidden="true">
          {slice.broad.map(({ d, fill }) => <path key={d} d={d} fill={fill} />)}
        </g>
        <g fillRule="evenodd">
          {layers.map(r => {
            const shape = slice.labels[r.id];
            if (!shape) return null;
            const selected = r.id === region.id, focused = hover === r.id;
            return <path key={r.id} d={shape.d} fill={effectColor(value(r))} className="region-shape"
              data-center={shape.center.join(",")} data-selected={selected} data-hover={focused} aria-label={r.name + " " + r.id}
              role="button" aria-pressed={selected} tabIndex={0} strokeWidth={selected || focused ? 1 : .42}
              onPointerEnter={() => setHover(r.id)} onPointerLeave={() => setHover(null)}
              onFocus={() => setHover(r.id)} onBlur={() => setHover(null)}
              onClick={() => {setHover(null); select(r.id);}}
              onKeyDown={e => {if (e.key === "Enter" || e.key === " ") {e.preventDefault(); setHover(null); select(r.id);}}}>
              <title>{r.name} · {r.english} ({r.id})</title>
            </path>;
          })}
        </g>
        {pointer ? <g className="region-callout" pointerEvents="none" aria-hidden="true">
          <circle cx={pointer.center[0]} cy={pointer.center[1]} r="3.2" />
          <circle cx={pointer.center[0]} cy={pointer.center[1]} r="1.5" />
        </g> : null}
        <g className="scale-bar" aria-hidden="true">
          <path d={"M " + (slice.w - 24) + " " + (slice.h - 2) + " h20"} />
          <text x={slice.w - 14} y={slice.h + 5} textAnchor="middle">1 mm</text>
        </g>
      </svg> : <p className="blank-chart">{loadingError ? "解剖切片暂不可用。可通过脑区索引查看结果。" : "正在载入解剖切片…"}</p>}
      <div className="anatomy-bottom">
        <span>{orientation[0]} <span className="axis-mark">—</span> {orientation[1]}</span>
        <span className="map-focus"><b>{regionToShow.id}</b> {regionToShow.name}{value(regionToShow) != null ? " · " + signed(value(regionToShow)) : ""}</span>
        <span>背侧 ↑</span>
      </div>
    </div>
    <div className="map-legend">
      <div><span className="scale-swatches" aria-hidden="true">{[-1.5, -.75, 0, .75, 1.5].map(x => <i key={x} style={{background: effectColor(x)}} />)}</span>
        <span>≤−1.5　0　≥+1.5 bits/s</span></div>
      <span><i className="missing-key"></i> 无该对比</span>
    </div>
    {data && <small className="atlas-credit">参考解剖切片，电极未配准。
      {" "}<a href={data.sourceInfo}>Allen Institute · CCFv3 (2017)</a>
      {" · "}<a href={data.paper}>Wang et al.</a>{" · "}<a href={data.license}>使用条款</a>
    </small>}
  </section>;
}
