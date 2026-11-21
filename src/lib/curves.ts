// Traces the critical curves and caustics of the lens in lens.ts. This runs
// when the site is built, on a much finer grid than a browser could afford at
// page load, and returns SVG path data.
import { HALF, deflect } from './lens';

const N = 2200; // grid nodes per side
const STEP = (2 * HALF) / (N - 1);

const MIN_LOOP = 0.007; // smallest loop kept, in image widths

type Point = [number, number];

/** Ramer–Douglas–Peucker simplification of an open polyline. */
function simplify(points: Point[], tolerance: number): Point[] {
  if (points.length < 3) return points;
  const last = points[points.length - 1];
  if (points[0][0] === last[0] && points[0][1] === last[1]) {
    // A closed loop has no chord to measure from, so split it at the point
    // farthest from its start and simplify the two halves.
    let far = 1;
    let best = -1;
    for (let i = 1; i < points.length - 1; i++) {
      const d = Math.hypot(points[i][0] - points[0][0], points[i][1] - points[0][1]);
      if (d > best) {
        best = d;
        far = i;
      }
    }
    return [...simplify(points.slice(0, far + 1), tolerance).slice(0, -1), ...simplify(points.slice(far), tolerance)];
  }
  const keep = new Uint8Array(points.length);
  keep[0] = keep[points.length - 1] = 1;
  const stack: [number, number][] = [[0, points.length - 1]];
  while (stack.length) {
    const [a, b] = stack.pop()!;
    const [ax, ay] = points[a];
    const [bx, by] = points[b];
    const dx = bx - ax;
    const dy = by - ay;
    const len = Math.hypot(dx, dy) || 1e-12;
    let worst = 0;
    let at = -1;
    for (let i = a + 1; i < b; i++) {
      const d = Math.abs(dy * (points[i][0] - ax) - dx * (points[i][1] - ay)) / len;
      if (d > worst) {
        worst = d;
        at = i;
      }
    }
    if (worst > tolerance && at > 0) {
      keep[at] = 1;
      stack.push([a, at], [at, b]);
    }
  }
  return points.filter((_, i) => keep[i]);
}

function toPath(lines: Point[][], tolerance: number): string {
  return lines
    .map((line) => {
      const pts = simplify(line, tolerance);
      return 'M' + pts.map(([x, y]) => `${x.toFixed(4)} ${y.toFixed(4)}`).join('L');
    })
    .join('');
}

function trace() {
  // Source-plane position of every grid node.
  const bx = new Float64Array(N * N);
  const by = new Float64Array(N * N);
  const a = [0, 0];
  for (let j = 0; j < N; j++) {
    const y = -HALF + j * STEP;
    for (let i = 0; i < N; i++) {
      const x = -HALF + i * STEP;
      deflect(x, y, a);
      bx[j * N + i] = x - a[0];
      by[j * N + i] = y - a[1];
    }
  }
  // Determinant of the lensing Jacobian; it changes sign on a critical curve.
  const det = new Float64Array(N * N);
  for (let j = 0; j < N; j++) {
    const jm = Math.max(j - 1, 0);
    const jp = Math.min(j + 1, N - 1);
    for (let i = 0; i < N; i++) {
      const im = Math.max(i - 1, 0);
      const ip = Math.min(i + 1, N - 1);
      const ex = 1 / ((ip - im) * STEP);
      const ey = 1 / ((jp - jm) * STEP);
      const axx = (bx[j * N + ip] - bx[j * N + im]) * ex;
      const ayx = (by[j * N + ip] - by[j * N + im]) * ex;
      const axy = (bx[jp * N + i] - bx[jm * N + i]) * ey;
      const ayy = (by[jp * N + i] - by[jm * N + i]) * ey;
      det[j * N + i] = axx * ayy - axy * ayx;
    }
  }

  // Marching squares. Each crossing lies on a grid edge, which gives it an id
  // that the two cells sharing the edge agree on, so segments can be chained.
  const pos = new Map<number, Point>();
  const links = new Map<number, number[]>();
  const crossing = (i0: number, j0: number, i1: number, j1: number): number => {
    const p0 = j0 * N + i0;
    const p1 = j1 * N + i1;
    const d0 = det[p0];
    const d1 = det[p1];
    if (d0 < 0 === d1 < 0) return -1;
    const id = 2 * p0 + (j1 === j0 ? 0 : 1); // edge leaving node p0 rightwards or downwards
    if (!pos.has(id)) {
      const t = d0 / (d0 - d1);
      pos.set(id, [-HALF + (i0 + t * (i1 - i0)) * STEP, -HALF + (j0 + t * (j1 - j0)) * STEP]);
    }
    return id;
  };
  const link = (u: number, v: number) => {
    (links.get(u) ?? links.set(u, []).get(u)!).push(v);
    (links.get(v) ?? links.set(v, []).get(v)!).push(u);
  };
  for (let j = 0; j < N - 1; j++) {
    for (let i = 0; i < N - 1; i++) {
      const ids = [
        crossing(i, j, i + 1, j),
        crossing(i + 1, j, i + 1, j + 1),
        crossing(i, j + 1, i + 1, j + 1),
        crossing(i, j, i, j + 1),
      ].filter((id) => id >= 0);
      if (ids.length === 2) link(ids[0], ids[1]);
      else if (ids.length === 4) {
        link(ids[0], ids[1]);
        link(ids[2], ids[3]);
      }
    }
  }

  // Chain the segments into polylines, open ends first and then closed loops.
  const used = new Set<number>();
  const chains: number[][] = [];
  const walk = (start: number) => {
    const chain = [start];
    used.add(start);
    let cur = start;
    for (;;) {
      const next = (links.get(cur) ?? []).find((n) => !used.has(n));
      if (next === undefined) break;
      used.add(next);
      chain.push(next);
      cur = next;
    }
    // Close the loop when the walk ended next to where it began.
    if (chain.length > 2 && (links.get(cur) ?? []).includes(start)) chain.push(start);
    return chain;
  };
  for (const [id, n] of links) if (n.length === 1 && !used.has(id)) chains.push(walk(id));
  for (const id of links.keys()) if (!used.has(id)) chains.push(walk(id));

  const critical: Point[][] = [];
  const caustic: Point[][] = [];
  for (const chain of chains) {
    const line = chain.map((id) => pos.get(id)!);
    // Skip loops too small to read on screen (the cores of member galaxies).
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const [x, y] of line) {
      if (x < x0) x0 = x;
      if (x > x1) x1 = x;
      if (y < y0) y0 = y;
      if (y > y1) y1 = y;
    }
    if (Math.hypot(x1 - x0, y1 - y0) < MIN_LOOP) continue;
    critical.push(line);
    // A caustic is the critical curve seen through the lens equation.
    caustic.push(
      line.map(([x, y]) => {
        deflect(x, y, a);
        return [x - a[0], y - a[1]] as Point;
      }),
    );
  }
  return {
    critical: toPath(critical, STEP * 0.6),
    caustic: toPath(caustic, STEP * 0.4),
    loops: critical.length,
  };
}

let cached: ReturnType<typeof trace> | undefined;
export function lensCurves() {
  return (cached ??= trace());
}
