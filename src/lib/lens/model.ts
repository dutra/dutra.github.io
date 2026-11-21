// The cluster mass model: one or two cluster-scale halos plus the member
// galaxies, as in parametric cluster lens models.
//
// Parameters are physical (velocity dispersions in km/s, radii in arcsec).
// Positions and deflections are worked out in units of the photograph's width,
// measured from the brightest cluster galaxy, x to the right and y down, which
// is what the renderers want. All deflections here are for D_ls/D_s = 1; scale
// them by the lensing efficiency of a source to get what that source sees.
//
// This is an illustration tied to the photograph, not a fit to the cluster.
import membersData from '../../data/lens-members.json';
import { sisEinsteinRadius } from './cosmology';

/** Redshift of Abell S1063. */
export const Z_LENS = 0.348;
/** Width of the photograph on the sky: ESA gives 2.31′ across 6920 px, and the site's crop keeps 6300 of them. */
export const ARCSEC_PER_WIDTH = (2.31 * 60 * 6300) / 6920;
/** Position of the brightest cluster galaxy, as a fraction of the photograph. */
export const BCG = { x: 0.5749, y: 0.5018 };
/** Photograph height divided by its width. */
export const IMAGE_ASPECT = 3360 / 3200;

export interface Halo {
  /** Velocity dispersion, km/s. */
  sigma: number;
  /** Core radius, arcsec. */
  core: number;
  /** Axis ratio of the mass, short over long (below 1). */
  q: number;
  /** Direction of the long axis, degrees from the x axis (y points down). */
  angle: number;
  /** Centre, arcsec from the BCG. */
  x: number;
  y: number;
}

export interface Model {
  halo: Halo;
  /** A second cluster-scale clump, as in a merging cluster. */
  halo2: Halo & { on: boolean };
  /** Velocity dispersion of a member as bright as the reference galaxy, km/s. */
  memberSigma: number;
  /** Truncation and core radius of that member, arcsec. */
  memberCut: number;
  memberCore: number;
}

export const DEFAULT_MODEL: Model = {
  halo: { sigma: 1160, core: 6.3, q: 0.7, angle: -30, x: 0, y: 0 },
  halo2: { on: false, sigma: 600, core: 5, q: 0.8, angle: -30, x: 30, y: -24 },
  memberSigma: 180,
  memberCut: 3.8,
  memberCore: 0.1,
};

export const cloneModel = (m: Model): Model => ({ ...m, halo: { ...m.halo }, halo2: { ...m.halo2 } });

type Vec = Float64Array | Float32Array | number[];

const widths = (arcsec: number) => arcsec / ARCSEC_PER_WIDTH;

/** A halo with its constants worked out, in width units. */
interface PreparedHalo {
  x: number;
  y: number;
  s: number;
  q: number;
  cos: number;
  sin: number;
  ecc: number;
  norm: number;
}

function prepareHalo(h: Halo): PreparedHalo {
  const q = Math.min(Math.max(h.q, 0.2), 0.995);
  const ecc = Math.sqrt(1 - q * q);
  const phi = (h.angle * Math.PI) / 180;
  // Kormann et al. (1994) normalisation: the mass inside an isodensity ellipse
  // matches a singular isothermal sphere of the same velocity dispersion.
  const b = widths(sisEinsteinRadius(h.sigma)) / Math.sqrt(q);
  return {
    x: widths(h.x),
    y: widths(h.y),
    s: widths(h.core),
    q,
    cos: Math.cos(phi),
    sin: Math.sin(phi),
    ecc,
    norm: (b * q) / ecc,
  };
}

/** Adds the deflection of a cored isothermal ellipsoid (Kormann et al. 1994; Keeton & Kochanek 1998) to `out`. */
function addHalo(h: PreparedHalo, x: number, y: number, out: Vec): void {
  const dx = x - h.x;
  const dy = y - h.y;
  const xr = h.cos * dx + h.sin * dy;
  const yr = -h.sin * dx + h.cos * dy;
  const psi = Math.sqrt(h.q * h.q * (h.s * h.s + xr * xr) + yr * yr);
  const axr = h.norm * Math.atan((h.ecc * xr) / (psi + h.s));
  const ayr = h.norm * Math.atanh((h.ecc * yr) / (psi + h.q * h.q * h.s));
  out[0] += h.cos * axr - h.sin * ayr;
  out[1] += h.sin * axr + h.cos * ayr;
}

/** Member galaxies: position in width units and brightness relative to the reference galaxy. */
export const MEMBERS = membersData.members.map((m) => ({
  x: m.x - BCG.x,
  y: (m.y - BCG.y) * IMAGE_ASPECT,
  l: m.l,
}));

/**
 * Deflection of the cluster-scale halos alone at (x, y), written into `out`.
 * Smooth everywhere, so it can be tabulated on a coarse grid.
 */
export function makeHaloDeflector(model: Model): (x: number, y: number, out: Vec) => void {
  const halos = [prepareHalo(model.halo)];
  if (model.halo2.on) halos.push(prepareHalo(model.halo2));
  return (x, y, out) => {
    out[0] = 0;
    out[1] = 0;
    for (const h of halos) addHalo(h, x, y, out);
  };
}

/**
 * Deflection of all member galaxies at (x, y), written into `out`.
 *
 * Each is a truncated isothermal sphere (dual pseudo-isothermal, Elíasdóttir et
 * al. 2007). Velocity dispersion scales as L^(1/4), so the lensing strength
 * and the radii scale as L^(1/2). The result is proportional to memberSigma².
 */
export function makeMemberDeflector(model: Model): (x: number, y: number, out: Vec) => void {
  const b0 = widths(sisEinsteinRadius(model.memberSigma));
  const cut0 = widths(model.memberCut);
  const core0 = widths(model.memberCore);
  const n = MEMBERS.length;
  const mx = new Float64Array(n);
  const my = new Float64Array(n);
  const amp = new Float64Array(n);
  const cut = new Float64Array(n);
  const core = new Float64Array(n);
  MEMBERS.forEach((m, i) => {
    const k = Math.sqrt(m.l);
    mx[i] = m.x;
    my[i] = m.y;
    cut[i] = cut0 * k;
    core[i] = core0 * k;
    amp[i] = (b0 * k * cut[i]) / (cut[i] - core[i]);
  });
  return (x, y, out) => {
    let ax = 0;
    let ay = 0;
    for (let i = 0; i < n; i++) {
      const dx = x - mx[i];
      const dy = y - my[i];
      const r2 = dx * dx + dy * dy;
      if (r2 < 1e-14) continue;
      const c = core[i];
      const t = cut[i];
      const f = (amp[i] * (Math.sqrt(c * c + r2) - c - Math.sqrt(t * t + r2) + t)) / r2;
      ax += f * dx;
      ay += f * dy;
    }
    out[0] = ax;
    out[1] = ay;
  };
}

/** Full deflection of the model at (x, y) for D_ls/D_s = 1. */
export function makeDeflector(model: Model): (x: number, y: number, out: Vec) => void {
  const halo = makeHaloDeflector(model);
  const members = makeMemberDeflector(model);
  const tmp = [0, 0];
  return (x, y, out) => {
    members(x, y, tmp);
    halo(x, y, out);
    out[0] += tmp[0];
    out[1] += tmp[1];
  };
}
