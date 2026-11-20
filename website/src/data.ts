export type Species = "mouse" | "rat";
export type State = "Wake" | "NREM" | "REM";
export type Contrast = "NW" | "RN" | "RW";
export type View = Contrast | "all";
export const STATES: State[] = ["Wake", "NREM", "REM"];
export const PAIRS: Record<Contrast, [State, State]> = {
  NW: ["Wake", "NREM"], RN: ["NREM", "REM"], RW: ["Wake", "REM"],
};
export const DELTA: Record<Contrast, string> = {
  NW: "NREM − Wake", RN: "REM − NREM", RW: "REM − Wake",
};
export const WARNING: Record<string, string> = {
  kurtosis: "残差超额峰度超过描述阈值", residualACF: "残差存在短滞后自相关", modelACF: "观测与拟合 ACF 偏差较大",
};
export const SPEC: Record<string, string> = {
  reference: "参考 AR(5)", AR10: "AR(10)", equal_duration: "等时长",
  edge_trim: "状态边界裁剪", strict_units: "严格单元", rate_thinning: "配平平均放电率",
  unit_count: "配平单元数",
};
export const STATE_COLOR: Record<State, string> = { Wake: "#61768c", NREM: "#af6b28", REM: "#567c88" };
export const SPECIES: Record<Species, string> = { mouse: "小鼠", rat: "大鼠" };
export const NO_PAIR = "无已核验的状态对比";
export type Stats = { n: number; sessions: number; mean: number | null; positive: number; negative: number };
export type CohortSummary = {
  id: string; dataset: string; context: string; contextLabel: string;
  animals: number | null; sessions: number; states: State[]; specs: string[];
  effects: Record<string, Record<Contrast, Stats>>; unitRange: number[];
};
export type Fit = {
  value: number; order: number; units: number; seconds: number; root: number; exit: number; solver: number;
  kurtosis: number | null; residual: number | null; rmse: number | null; silence: number | null;
  rate: number | null; warnings: string[]; score?: string; provenance: string;
};
export type Session = {
  id: string; subject: string | null; session: string; asset: string; digest: string | Record<string, string>;
  fits: Record<string, Partial<Record<State, Fit>>>;
};
export type Control = {
  model: string; contrast: Contrast; subject: string; asset: string; reference: number;
  value: number; repeats: number; range: number[]; run: string; source: string;
};
export type Cohort = CohortSummary & { records: Session[]; controls: Control[] };
export type Region = {
  id: string; species: Species; name: string; english: string; system: string;
  allen: string | null; about: string | null; cohorts: CohortSummary[]; path: string | null;
};
export type Dataset = {
  id: string; species: Species; title: string; version: string; url: string; paper: string;
  report: string; manifest: string; recording: string; scoring: string; wake: string;
  limitations: string; access: string; animals: number | null; sessions: number;
};
export type Provenance = {
  run: string; file: string; path: string; sha256: string; toolkit: string; commit: string; matlab: string;
  binSeconds: number; beta: number; method: string; units: string;
};
export type Audit = {
  dataset: string; region?: string; session: string; concern: string; status: string;
  inference: string; evidence?: string; category?: string;
  original?: Partial<Record<Contrast, number>>; audited?: Partial<Record<Contrast, number>>;
};
export type Atlas = {
  schema: number; baseline: string; datasets: Record<string, Dataset>;
  provenance: Record<string, Provenance>; audit: Audit[];
};
export type Plane = "sagittal" | "coronal";
export type Slice = {
  w: number; h: number; outline: string; broad: { d: string; fill: string }[];
  labels: Record<string, { d: string; center: [number, number]; area: number }>;
};
export type Anatomy = {
  title: string; source: string; sourceInfo: string; ontology: string; paper: string; license: string;
  voxelMm: number; views: Record<Plane, { positions: Record<string, number>; slices: Record<string, Slice> }>;
};
export type Replicate = {
  id: string; values: Partial<Record<State, number>>;
  effect: number | null; count: number; sessions: Session[];
};

export async function json<T>(path: string, signal?: AbortSignal): Promise<T> {
  const res = await fetch(import.meta.env.BASE_URL + "data/" + path, { signal });
  if (!res.ok) throw new Error(path + " · " + res.status);
  return res.json() as Promise<T>;
}

export function average(numbers: number[]) {
  return numbers.reduce((sum, n) => sum + n, 0) / numbers.length;
}

export function paired(records: Session[], variant: string, states: State[]): Replicate[] {
  const groups = new Map<string, Session[]>();
  for (const r of records) {
    const id = r.subject ?? r.session;
    const group = groups.get(id);
    if (group) group.push(r);
    else groups.set(id, [r]);
  }
  return Array.from(groups).sort(([a], [b]) => a.localeCompare(b)).flatMap(([id, group]) => {
    const matches = group.filter(r => states.every(s => r.fits[variant]?.[s] !== undefined));
    if (!matches.length) return [];
    const values: Partial<Record<State, number>> = {};
    states.forEach(s => { values[s] = average(matches.map(r => r.fits[variant][s]!.value)); });
    const effect = states.length === 2 ? values[states[1]]! - values[states[0]]! : null;
    return [{ id, values, effect, count: matches.length, sessions: matches }];
  });
}

export function effectColor(value: number | null | undefined) {
  if (value == null) return "#d9e1e6";
  const strength = Math.min(Math.abs(value) / 1.5, 1);
  const end = value < 0 ? [66, 112, 152] : [186, 102, 39];
  const base = [243, 237, 222];
  return "rgb(" + base.map((b, i) => Math.round(b + (end[i] - b) * strength)).join(",") + ")";
}
export const amount = (x: number | null | undefined, digits = 3) => x == null ? "—" : x.toFixed(digits);
export const signed = (x: number | null | undefined, digits = 3) =>
  x == null ? "—" : (x > 0 ? "+" : x < 0 ? "−" : "") + Math.abs(x).toFixed(digits);
