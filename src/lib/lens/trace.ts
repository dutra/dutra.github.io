// Critical curves and caustics from the tabulated deflection field, quick
// enough to redo whenever the model or the source redshift changes.
//
// The determinant of the lensing Jacobian changes sign on a critical curve, so
// marching squares over its sign gives the curve as short segments. The same
// segments carried through the lens equation are the caustics.
import type { LensField } from './field';

export interface Curves {
  /** Segments as x0, y0, x1, y1 in width units: lens plane and source plane. */
  critical: Float32Array;
  caustic: Float32Array;
}

/** Curves smaller than this many grid cells across are dropped: they are the unresolved cores of member galaxies. */
const MIN_CELLS = 7;

let det = new Float32Array(0);
let label = new Int32Array(0);

export function traceCurves(field: LensField, eff: number): Curves {
  const { gw, gh, h, x0, y0, a } = field;
  const n = gw * gh;
  if (det.length !== n) {
    det = new Float32Array(n);
    label = new Int32Array(n);
  }

  const e = eff / (2 * h);
  for (let j = 0; j < gh; j++) {
    const up = j > 0 ? -2 * gw : 0;
    const down = j < gh - 1 ? 2 * gw : 0;
    const ky = up && down ? e : 2 * e;
    for (let i = 0; i < gw; i++) {
      const p = 2 * (j * gw + i);
      const left = i > 0 ? -2 : 0;
      const right = i < gw - 1 ? 2 : 0;
      const kx = left && right ? e : 2 * e;
      const axx = (a[p + right] - a[p + left]) * kx;
      const ayx = (a[p + right + 1] - a[p + left + 1]) * kx;
      const axy = (a[p + down] - a[p + up]) * ky;
      const ayy = (a[p + down + 1] - a[p + up + 1]) * ky;
      det[j * gw + i] = (1 - axx) * (1 - ayy) - axy * ayx;
    }
  }

  // Cells the curve passes through, grouped into connected pieces so the tiny ones can be dropped.
  const crossed = (c: number) => {
    const s = det[c] < 0;
    return s !== det[c + 1] < 0 || s !== det[c + gw] < 0 || s !== det[c + gw + 1] < 0;
  };
  label.fill(0);
  const keep: number[] = [];
  const stack: number[] = [];
  let next = 1;
  for (let j = 0; j < gh - 1; j++) {
    for (let i = 0; i < gw - 1; i++) {
      const c = j * gw + i;
      if (label[c] || !crossed(c)) continue;
      const id = next++;
      const first = keep.length;
      let i0 = i, i1 = i, j0 = j, j1 = j;
      label[c] = id;
      stack.push(c);
      while (stack.length) {
        const cur = stack.pop()!;
        keep.push(cur);
        const ci = cur % gw;
        const cj = (cur / gw) | 0;
        if (ci < i0) i0 = ci;
        if (ci > i1) i1 = ci;
        if (cj < j0) j0 = cj;
        if (cj > j1) j1 = cj;
        for (let dj = -1; dj <= 1; dj++) {
          for (let di = -1; di <= 1; di++) {
            const ni = ci + di;
            const nj = cj + dj;
            if (ni < 0 || nj < 0 || ni >= gw - 1 || nj >= gh - 1) continue;
            const nc = nj * gw + ni;
            if (!label[nc] && crossed(nc)) {
              label[nc] = id;
              stack.push(nc);
            }
          }
        }
      }
      if (Math.hypot(i1 - i0, j1 - j0) < MIN_CELLS) keep.length = first;
    }
  }

  const critical = new Float32Array(keep.length * 8);
  const caustic = new Float32Array(keep.length * 8);
  let len = 0;
  // The four edges of a cell: node offsets (in nodes) of each end.
  const EI = [0, 1, 1, 0, 0, 1, 0, 0];
  const EJ = [0, 0, 0, 1, 1, 1, 0, 1];
  const px = [0, 0, 0, 0];
  const py = [0, 0, 0, 0];
  const qx = [0, 0, 0, 0];
  const qy = [0, 0, 0, 0];
  for (const c of keep) {
    const ci = c % gw;
    const cj = (c / gw) | 0;
    let found = 0;
    for (let k = 0; k < 4; k++) {
      const ia = ci + EI[2 * k];
      const ja = cj + EJ[2 * k];
      const ib = ci + EI[2 * k + 1];
      const jb = cj + EJ[2 * k + 1];
      const na = ja * gw + ia;
      const nb = jb * gw + ib;
      const da = det[na];
      const db = det[nb];
      if (da < 0 === db < 0) continue;
      const t = da / (da - db);
      const x = x0 + (ia + t * (ib - ia)) * h;
      const y = y0 + (ja + t * (jb - ja)) * h;
      px[found] = x;
      py[found] = y;
      // The deflection varies linearly along the edge, like the interpolation used everywhere else.
      qx[found] = x - eff * (a[2 * na] + t * (a[2 * nb] - a[2 * na]));
      qy[found] = y - eff * (a[2 * na + 1] + t * (a[2 * nb + 1] - a[2 * na + 1]));
      found++;
    }
    for (let k = 0; k + 1 < found; k += 2) {
      critical[len] = px[k];
      critical[len + 1] = py[k];
      critical[len + 2] = px[k + 1];
      critical[len + 3] = py[k + 1];
      caustic[len] = qx[k];
      caustic[len + 1] = qy[k];
      caustic[len + 2] = qx[k + 1];
      caustic[len + 3] = qy[k + 1];
      len += 4;
    }
  }
  return { critical: critical.subarray(0, len), caustic: caustic.subarray(0, len) };
}
