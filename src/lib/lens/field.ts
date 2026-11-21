// The deflection field of the model, tabulated on a grid covering the whole
// photograph. Everything interactive reads from this table: the WebGL
// renderer samples it per pixel, and the curve tracer and readouts use it on
// the CPU, so the physics lives in one place (model.ts).
//
// The member galaxies need the full grid but never move, so their field is
// built once. The cluster halos are smooth, so they are tabulated on a grid
// four times coarser and interpolated, which makes moving a slider cheap.
import { BCG, IMAGE_ASPECT, makeHaloDeflector, makeMemberDeflector, type Model } from './model';

const COARSE = 4;

export class LensField {
  /** Grid nodes across and down. */
  readonly gw: number;
  readonly gh: number;
  /** Lens-plane position of node (0, 0) and the node spacing, in width units. */
  readonly x0 = -BCG.x;
  readonly y0 = -BCG.y * IMAGE_ASPECT;
  readonly h: number;
  /** Deflection for D_ls/D_s = 1, interleaved (ax, ay) per node. */
  readonly a: Float32Array;

  private members: Float32Array; // field of the members for memberSigmaRef
  private memberSigmaRef = 0;
  private memberShape = '';
  private cw: number;
  private ch: number;
  private coarse: Float32Array;

  constructor(gw: number) {
    // Both counts are 1 more than a multiple of COARSE so the coarse grid lands on fine nodes.
    this.gw = Math.ceil((gw - 1) / COARSE) * COARSE + 1;
    this.h = 1 / (this.gw - 1);
    this.gh = Math.ceil(IMAGE_ASPECT / this.h / COARSE) * COARSE + 1;
    this.a = new Float32Array(2 * this.gw * this.gh);
    this.members = new Float32Array(2 * this.gw * this.gh);
    this.cw = (this.gw - 1) / COARSE + 1;
    this.ch = (this.gh - 1) / COARSE + 1;
    this.coarse = new Float32Array(2 * this.cw * this.ch);
  }

  /** Tabulates the member galaxies. Slow (tens of millions of evaluations), so it yields between slices. */
  async buildMembers(model: Model, onProgress?: (fraction: number) => void): Promise<void> {
    const deflect = makeMemberDeflector(model);
    const out = [0, 0];
    const { gw, gh, h, x0, y0, members } = this;
    let last = performance.now();
    for (let j = 0; j < gh; j++) {
      const y = y0 + j * h;
      let p = 2 * j * gw;
      for (let i = 0; i < gw; i++) {
        deflect(x0 + i * h, y, out);
        members[p++] = out[0];
        members[p++] = out[1];
      }
      if (performance.now() - last > 24) {
        onProgress?.(j / gh);
        await new Promise((r) => setTimeout(r));
        last = performance.now();
      }
    }
    this.memberSigmaRef = model.memberSigma;
    this.memberShape = `${model.memberCut}/${model.memberCore}`;
  }

  /** Recomputes the halos and combines them with the members. Fast enough to call while a slider moves. */
  update(model: Model): void {
    if (this.memberShape !== `${model.memberCut}/${model.memberCore}`) {
      throw new Error('LensField: call buildMembers() first, and again if the member radii change');
    }
    const halo = makeHaloDeflector(model);
    const out = [0, 0];
    const { gw, gh, h, x0, y0, a, members, cw, ch, coarse } = this;
    const step = h * COARSE;
    for (let j = 0, p = 0; j < ch; j++) {
      for (let i = 0; i < cw; i++) {
        halo(x0 + i * step, y0 + j * step, out);
        coarse[p++] = out[0];
        coarse[p++] = out[1];
      }
    }
    // The members' field scales with the square of their velocity dispersion.
    const k = (model.memberSigma / this.memberSigmaRef) ** 2;
    const inv = 1 / COARSE;
    for (let j = 0; j < gh; j++) {
      const cj = Math.min((j / COARSE) | 0, ch - 2);
      const fy = (j - cj * COARSE) * inv;
      const r0 = 2 * cj * cw;
      const r1 = r0 + 2 * cw;
      let p = 2 * j * gw;
      for (let i = 0; i < gw; i++) {
        const ci = Math.min((i / COARSE) | 0, cw - 2);
        const fx = (i - ci * COARSE) * inv;
        const c00 = r0 + 2 * ci;
        const c01 = r1 + 2 * ci;
        const w00 = (1 - fx) * (1 - fy);
        const w10 = fx * (1 - fy);
        const w01 = (1 - fx) * fy;
        const w11 = fx * fy;
        a[p] = k * members[p] + w00 * coarse[c00] + w10 * coarse[c00 + 2] + w01 * coarse[c01] + w11 * coarse[c01 + 2];
        a[p + 1] =
          k * members[p + 1] + w00 * coarse[c00 + 1] + w10 * coarse[c00 + 3] + w01 * coarse[c01 + 1] + w11 * coarse[c01 + 3];
        p += 2;
      }
    }
  }

  /** Deflection at any lens-plane point for D_ls/D_s = 1, by bilinear interpolation. */
  sample(x: number, y: number, out: number[]): void {
    const { gw, gh, h, a } = this;
    const gx = Math.min(Math.max((x - this.x0) / h, 0), gw - 1.001);
    const gy = Math.min(Math.max((y - this.y0) / h, 0), gh - 1.001);
    const i = gx | 0;
    const j = gy | 0;
    const fx = gx - i;
    const fy = gy - j;
    const p = 2 * (j * gw + i);
    const q = p + 2 * gw;
    out[0] = (1 - fy) * ((1 - fx) * a[p] + fx * a[p + 2]) + fy * ((1 - fx) * a[q] + fx * a[q + 2]);
    out[1] = (1 - fy) * ((1 - fx) * a[p + 1] + fx * a[p + 3]) + fy * ((1 - fx) * a[q + 1] + fx * a[q + 3]);
  }

  /**
   * Determinant of the lensing Jacobian at the node nearest (x, y) for a
   * source with lensing efficiency `eff`. The magnification there is 1/det;
   * its sign is the parity of the image.
   */
  detAt(x: number, y: number, eff: number): number {
    const { gw, gh, h, a } = this;
    const i = Math.min(Math.max(Math.round((x - this.x0) / h), 1), gw - 2);
    const j = Math.min(Math.max(Math.round((y - this.y0) / h), 1), gh - 2);
    const p = 2 * (j * gw + i);
    const e = eff / (2 * h);
    const axx = (a[p + 2] - a[p - 2]) * e;
    const ayx = (a[p + 3] - a[p - 1]) * e;
    const axy = (a[p + 2 * gw] - a[p - 2 * gw]) * e;
    const ayy = (a[p + 2 * gw + 1] - a[p - 2 * gw + 1]) * e;
    return (1 - axx) * (1 - ayy) - axy * ayx;
  }

  /**
   * Effective Einstein radius about the BCG, in width units: the radius at
   * which the mean convergence inside it is 1. By the divergence theorem that
   * is where the azimuthal average of the radial deflection equals the radius.
   */
  einsteinRadius(eff: number): number {
    const out = [0, 0];
    const angles = 72;
    let prev = 0;
    let prevR = 0;
    const max = Math.min(-this.x0, 1 + this.x0, -this.y0, IMAGE_ASPECT + this.y0);
    for (let r = 0.004; r < max; r += 0.004) {
      let sum = 0;
      for (let k = 0; k < angles; k++) {
        const c = Math.cos((2 * Math.PI * k) / angles);
        const s = Math.sin((2 * Math.PI * k) / angles);
        this.sample(r * c, r * s, out);
        sum += out[0] * c + out[1] * s;
      }
      const excess = (eff * sum) / angles - r; // positive inside the Einstein radius
      if (excess < 0) return prevR ? prevR + ((r - prevR) * prev) / (prev - excess) : 0;
      prev = excess;
      prevR = r;
    }
    return max;
  }

  /**
   * Total magnification of a round source of radius `size` at (sx, sy): the
   * area of the lens plane that maps into it, over its own area.
   */
  totalMagnification(sx: number, sy: number, size: number, eff: number): number {
    const { gw, gh, h, a } = this;
    const r2 = size * size;
    let count = 0;
    for (let j = 0, p = 0; j < gh; j++) {
      const y = this.y0 + j * h - sy;
      for (let i = 0; i < gw; i++, p += 2) {
        const dy = y - eff * a[p + 1];
        if (dy > size || dy < -size) continue;
        const dx = this.x0 + i * h - sx - eff * a[p];
        if (dx * dx + dy * dy < r2) count++;
      }
    }
    return (count * h * h) / (Math.PI * r2);
  }
}
