// A self-interacting dark matter (SIDM) halo over time, and its strength as a
// gravitational lens. Pure functions, no DOM.
//
// Model: the parametric SIDM halo of Yang, Nadler, Yu & Zhong (2024, JCAP;
// arXiv:2305.16176). An NFW halo is followed through gravothermal evolution as
// a function of τ = t / t_c, where t_c is the core-collapse time:
//
//   Eq. (4)  ρ(r) = ρ_s / { [(r^β + r_c^β)^(1/β) / r_s] (1 + r / r_s)² },  β = 4
//
//   Eq. (5)  t_c = (150 / C) · 1 / [(σ/m) ρ_s,0 r_s,0] · 1 / √(4π G ρ_s,0),  C = 0.75
//            (Balberg, Shapiro & Inagaki 2002; Pollack, Spergel & Steinhardt 2015;
//             Essig, McDermott, Yu & Zhong 2019)
//
//   Eq. (6)  ρ_s / ρ_s,0 = 2.033 + 0.7381 τ + 7.264 τ⁵ − 12.73 τ⁷ + 9.915 τ⁹
//                          + (1 − 2.033) ln(τ + 0.001) / ln 0.001
//            r_s / r_s,0 = 0.7178 − 0.1026 τ + 0.2474 τ² − 0.4079 τ³
//                          + (1 − 0.7178) ln(τ + 0.001) / ln 0.001
//            r_c / r_s,0 = 2.555 √τ − 3.632 τ + 2.131 τ² − 1.415 τ³ + 0.4683 τ⁴
//
// The fits come from an isolated N-body SIDM halo followed up to its collapse,
// so they are compared with the simulation only over 0 ≤ τ ≤ 1. They reduce to the
// initial NFW halo at τ = 0 (the log terms are exactly 1 − 2.033 and 1 − 0.7178
// there). Like the paper, we stop at τ = TAU_MAX = 1 to avoid extrapolation (the
// polynomials make the central density jump by about ×2.35 between τ = 1 and 1.1): past t_c the
// halo is held at its τ = 1 profile and counts as collapsed.
//
// Units throughout: M☉, kpc, km/s, Gyr; σ/m in cm²/g.

import { C_KMS, H0, OMEGA_M, comoving, efficiency } from './lens/cosmology';

/** Newton's constant in kpc (km/s)² / M☉. */
export const G = 4.30091e-6;
/** One kpc / (km/s) in Gyr: 3.0857e16 km ÷ (3.15576e16 s/Gyr). */
export const GYR_PER_KPC_KMS = 0.977792;
/**
 * One cm²/g in kpc²/M☉: 1 cm²/g × 1.98841e33 g/M☉ ÷ (3.0856776e21 cm/kpc)²
 *   = 1.98841e33 ÷ 9.5214062e42 = 2.08836e-10.
 */
export const CM2G = 2.08836e-10;
/** The constant C in Eq. (5), fixed by Yang et al. (2024) and Essig et al. (2019). */
export const C_COLLAPSE = 0.75;
/** Sharpness of the core in Eq. (4). */
export const BETA = 4;
/** Largest τ the model is used at: the end of the fitted range. */
export const TAU_MAX = 1;
/** Lens redshift of the page; the halo is defined with the critical density at this redshift. */
export const Z_LENS = 0.348;
/** Line-of-sight truncation of the projected integrals, in units of r_200. */
export const TRUNCATION = 2;

/** Parameters of the profile in Eq. (4); rC = 0 is an NFW halo. */
export interface Profile {
  /** Scale density, M☉/kpc³. */
  rhoS: number;
  /** Scale radius, kpc. */
  rS: number;
  /** Core radius, kpc. */
  rC: number;
}

/** The initial NFW halo. */
export interface Halo extends Profile {
  m200: number;
  c: number;
  /** Radius enclosing 200 times the critical density, kpc. */
  r200: number;
}

// --- cosmology and the initial halo ----------------------------------------------

/** Critical density at redshift z, M☉/kpc³: 3 H(z)² / (8π G), with H in km/s/kpc. */
export function criticalDensity(z: number): number {
  const h = (H0 / 1000) * Math.sqrt(OMEGA_M * (1 + z) ** 3 + 1 - OMEGA_M);
  return (3 * h * h) / (8 * Math.PI * G);
}

/**
 * Age of the Universe at redshift z, Gyr, for the flat ΛCDM of cosmology.ts (which has
 * no age function): t(z) = (1/H0) ∫_z^∞ dz′ / [(1 + z′) E(z′)] = (1/H0) ∫_0^a da / (a E(a)),
 * a = 1/(1 + z). With a = x² the integrand becomes 2x² / √(Ω_m + Ω_Λ x⁶), smooth at x = 0,
 * integrated by Simpson's rule. 1/H0 = 977.792 / H0 Gyr for H0 in km/s/Mpc.
 */
export function ageAt(z: number): number {
  const xMax = 1 / Math.sqrt(1 + z);
  const f = (x: number) => (2 * x * x) / Math.sqrt(OMEGA_M + (1 - OMEGA_M) * x ** 6);
  return ((GYR_PER_KPC_KMS * 1000) / H0) * simpson(f, 0, xMax, 512);
}

/** Age of the Universe at the lens redshift, Gyr: the longest a halo seen there can have evolved. */
export const T_LENS = ageAt(Z_LENS);

/** Length of the time axis, Gyr: up to 1.25 t_c so short collapse times are resolved, never past T_LENS. */
export function timeSpan(tc: number): number {
  return Number.isFinite(tc) ? Math.min(T_LENS, 1.25 * tc) : T_LENS;
}

/**
 * Where a time t (Gyr) goes when the span changes. A time that still fits is kept; one past
 * the new end of the axis goes to t_c, the end of the model, rather than to the end of the
 * axis (which lies past the collapse). Without a collapse time it goes to the end of the axis.
 */
export function clampTime(t: number, tc: number, span: number): number {
  if (!(t > 0)) return 0;
  return t > span ? Math.min(tc, span) : t;
}

// --- the cross-section slider ----------------------------------------------------------

/** Range of σ/m offered, cm²/g. */
export const SIGMA_MIN = 0.1;
export const SIGMA_MAX = 300;
/** The σ/m slider is an integer 0…SIGMA_STEPS, logarithmic in σ/m. */
export const SIGMA_STEPS = 1000;

/** σ/m (cm²/g) at slider position i; the two ends are exact. */
export function sigmaAtStep(i: number): number {
  if (i <= 0) return SIGMA_MIN;
  if (i >= SIGMA_STEPS) return SIGMA_MAX;
  return SIGMA_MIN * (SIGMA_MAX / SIGMA_MIN) ** (i / SIGMA_STEPS);
}

/** The nearest slider position to σ/m (cm²/g), clamped to the slider. */
export function stepOfSigma(sigma: number): number {
  const f = Math.log(sigma / SIGMA_MIN) / Math.log(SIGMA_MAX / SIGMA_MIN);
  return Math.min(Math.max(Math.round(f * SIGMA_STEPS), 0), SIGMA_STEPS);
}

/** m(x) = ln(1 + x) − x / (1 + x): NFW mass inside x = r / r_s in units of 4π ρ_s r_s³. */
const nfwM = (x: number) => Math.log1p(x) - x / (1 + x);

/** The NFW halo of mass m200 (M☉) and concentration c, at redshift z. */
export function nfwHalo(m200: number, c: number, z = Z_LENS): Halo {
  const rhoC = criticalDensity(z);
  const r200 = Math.cbrt((3 * m200) / (4 * Math.PI * 200 * rhoC));
  const rS = r200 / c;
  const rhoS = ((200 / 3) * rhoC * c ** 3) / nfwM(c);
  return { m200, c, r200, rhoS, rS, rC: 0 };
}

/** Analytic NFW mass inside radius r, M☉. */
export function nfwMass(h: Profile, r: number): number {
  return 4 * Math.PI * h.rhoS * h.rS ** 3 * nfwM(r / h.rS);
}

/** Analytic NFW surface density at projected radius R, M☉/kpc² (Bartelmann 1996; untruncated). */
export function nfwSurfaceDensity(h: Profile, R: number): number {
  const x = R / h.rS;
  let f: number;
  if (Math.abs(x - 1) < 1e-6) f = 1 / 3;
  else if (x < 1) f = (1 - (2 / Math.sqrt(1 - x * x)) * Math.atanh(Math.sqrt((1 - x) / (1 + x)))) / (x * x - 1);
  else f = (1 - (2 / Math.sqrt(x * x - 1)) * Math.atan(Math.sqrt((x - 1) / (1 + x)))) / (x * x - 1);
  return 2 * h.rhoS * h.rS * f;
}

// --- the parametric model ----------------------------------------------------------

/** Core-collapse time, Gyr, Eq. (5). The product (σ/m) ρ_s r_s is dimensionless. */
export function collapseTime(h: Profile, sigmaM: number): number {
  const opacity = sigmaM * CM2G * h.rhoS * h.rS;
  const dynamical = 1 / Math.sqrt(4 * Math.PI * G * h.rhoS); // kpc / (km/s)
  return (150 / C_COLLAPSE / opacity) * dynamical * GYR_PER_KPC_KMS;
}

const LN_EPS = Math.log(0.001);

/** The three ratios of Eq. (6) at dimensionless time τ (clamped to 0…TAU_MAX). */
export function ratios(tau: number): { rhoS: number; rS: number; rC: number } {
  const t = Math.min(Math.max(tau, 0), TAU_MAX);
  const lg = Math.log(t + 0.001) / LN_EPS;
  const t2 = t * t;
  const t3 = t2 * t;
  const t5 = t3 * t2;
  return {
    rhoS: 2.033 + 0.7381 * t + 7.264 * t5 - 12.73 * t5 * t2 + 9.915 * t5 * t2 * t2 + (1 - 2.033) * lg,
    rS: 0.7178 - 0.1026 * t + 0.2474 * t2 - 0.4079 * t3 + (1 - 0.7178) * lg,
    rC: 2.555 * Math.sqrt(t) - 3.632 * t + 2.131 * t2 - 1.415 * t3 + 0.4683 * t2 * t2,
  };
}

/** The profile of halo h at dimensionless time τ. */
export function profileAt(h: Profile, tau: number): Profile {
  const q = ratios(tau);
  return { rhoS: h.rhoS * q.rhoS, rS: h.rS * q.rS, rC: h.rS * Math.max(q.rC, 0) };
}

/** Density at radius r, M☉/kpc³, Eq. (4) with β = 4. */
export function density(p: Profile, r: number): number {
  const rc4 = p.rC ** 4;
  const soft = rc4 > 0 ? Math.sqrt(Math.sqrt(r ** 4 + rc4)) : r;
  const y = 1 + r / p.rS;
  return p.rhoS / ((soft / p.rS) * y * y);
}

/** Central density ρ(0) = ρ_s r_s / r_c, M☉/kpc³; infinite for an NFW cusp. */
export function centralDensity(p: Profile): number {
  return p.rC > 0 ? (p.rhoS * p.rS) / p.rC : Infinity;
}

/** The τ at which the core radius peaks, and its size there in units of r_s,0 (found numerically). */
export const CORE_PEAK = (() => {
  let best = { tau: 0, rC: 0 };
  for (let i = 1; i <= 2000; i++) {
    const tau = (i / 2000) * TAU_MAX;
    const rC = ratios(tau).rC;
    if (rC > best.rC) best = { tau, rC };
  }
  return best;
})();

export type Phase = 'core forming' | 'core at its largest' | 'collapsing' | 'collapsed';

/** Evolutionary phase at τ: "largest" while r_c is within 10% of its peak. */
export function phase(tau: number): Phase {
  if (tau >= 1) return 'collapsed';
  if (ratios(tau).rC >= 0.9 * CORE_PEAK.rC) return 'core at its largest';
  return tau < CORE_PEAK.tau ? 'core forming' : 'collapsing';
}

// --- numerical helpers ---------------------------------------------------------------

/** Simpson's rule for f on [a, b] with n (even) intervals. */
function simpson(f: (x: number) => number, a: number, b: number, n: number): number {
  const h = (b - a) / n;
  let sum = f(a) + f(b);
  for (let i = 1; i < n; i++) sum += (i % 2 ? 4 : 2) * f(a + i * h);
  return (sum * h) / 3;
}

/** n log-spaced values from a to b inclusive. */
export function logspace(a: number, b: number, n: number): Float64Array {
  const out = new Float64Array(n);
  const la = Math.log(a);
  const step = (Math.log(b) - la) / (n - 1);
  for (let i = 0; i < n; i++) out[i] = Math.exp(la + i * step);
  return out;
}

// --- 3D mass and rotation curve ----------------------------------------------------------

/** Enclosed mass on a log-spaced radius grid, for interpolation. */
export interface MassTable {
  r: Float64Array;
  m: Float64Array;
}

const TABLE_MIN = 1e-6; // kpc
const TABLE_N = 800;

/**
 * M(<r) = ∫ 4π r² ρ dr, integrated as ∫ 4π r³ ρ d ln r by the trapezoid rule on
 * a log grid from 1e-6 kpc to rMax (step ≈ 0.03 in ln r, error well below 0.1%).
 * The mass inside the first point is taken as 2π ρ r³, exact for a 1/r cusp.
 */
export function massTable(p: Profile, rMax: number): MassTable {
  const r = logspace(TABLE_MIN, rMax, TABLE_N);
  const m = new Float64Array(TABLE_N);
  const dl = Math.log(r[1] / r[0]);
  let prev = 4 * Math.PI * r[0] ** 3 * density(p, r[0]);
  m[0] = prev / 2;
  for (let i = 1; i < TABLE_N; i++) {
    const cur = 4 * Math.PI * r[i] ** 3 * density(p, r[i]);
    m[i] = m[i - 1] + 0.5 * (prev + cur) * dl;
    prev = cur;
  }
  return { r, m };
}

/** M(<r) from a table, interpolated linearly in log M against log r. */
export function massAt(t: MassTable, r: number): number {
  const n = t.r.length;
  if (r <= t.r[0]) return t.m[0] * (r / t.r[0]) ** 2;
  if (r >= t.r[n - 1]) return t.m[n - 1];
  const f = (Math.log(r / t.r[0]) / Math.log(t.r[n - 1] / t.r[0])) * (n - 1);
  const i = Math.min(Math.floor(f), n - 2);
  const w = f - i;
  return Math.exp((1 - w) * Math.log(t.m[i]) + w * Math.log(t.m[i + 1]));
}

/** Circular velocity √(G M(<r) / r), km/s. */
export function circularVelocity(t: MassTable, r: number): number {
  return Math.sqrt((G * massAt(t, r)) / r);
}

/** Peak circular velocity (km/s) and its radius (kpc), refined with a parabola in ln r. */
export function vmax(t: MassTable): { vmax: number; rmax: number } {
  let best = 1;
  const v2 = (i: number) => t.m[i] / t.r[i];
  for (let i = 1; i < t.r.length - 1; i++) if (v2(i) > v2(best)) best = i;
  const i = Math.min(Math.max(best, 1), t.r.length - 2);
  const [a, b, c] = [v2(i - 1), v2(i), v2(i + 1)];
  const denom = a - 2 * b + c;
  const shift = denom < 0 ? Math.max(-1, Math.min(1, (0.5 * (a - c)) / denom)) : 0;
  const dl = Math.log(t.r[1] / t.r[0]);
  const rmax = t.r[i] * Math.exp(shift * dl);
  return { vmax: Math.sqrt(G * (b - 0.25 * (a - c) * shift)), rmax };
}

// --- projection and lensing --------------------------------------------------------------------

/** Projected quantities on a grid of projected radii R (kpc). */
export interface Projection {
  R: Float64Array;
  /** Surface density Σ(R), M☉/kpc². */
  sigma: Float64Array;
  /** Projected (cylindrical) mass inside R, M☉. */
  mass: Float64Array;
}

const PROJ_N = 256; // Simpson intervals per line of sight

/**
 * Σ(R) and the cylinder mass M_2D(<R) of a halo truncated at radius rt.
 * With r = R cosh u along the line of sight:
 *   Σ(R)     = 2 ∫₀^U ρ(R cosh u) R cosh u du,                       U = arccosh(rt / R)
 *   M_2D(<R) = M(<R) + ∫_R^rt 4π r² ρ (1 − √(1 − R²/r²)) dr
 *            = M(<R) + ∫₀^U 4π R³ ρ(R cosh u) cosh²u sinh u · 2 / (e^{2u} + 1) du,
 * the second form avoiding the cancellation in M(<rt) − (mass outside the cylinder).
 */
export function project(p: Profile, rt: number, R: Float64Array, table: MassTable): Projection {
  const n = R.length;
  const sigma = new Float64Array(n);
  const mass = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const Ri = R[i];
    if (Ri >= rt) {
      mass[i] = massAt(table, rt);
      continue;
    }
    const U = Math.acosh(rt / Ri);
    sigma[i] = 2 * Ri * simpson((u) => density(p, Ri * Math.cosh(u)) * Math.cosh(u), 0, U, PROJ_N);
    const shell = simpson(
      (u) => {
        const ch = Math.cosh(u);
        return density(p, Ri * ch) * ch * ch * Math.sinh(u) * (2 / (Math.exp(2 * u) + 1));
      },
      0,
      U,
      PROJ_N,
    );
    mass[i] = massAt(table, Ri) + 4 * Math.PI * Ri ** 3 * shell;
  }
  return { R, sigma, mass };
}

/** c² / (4π G) in M☉/kpc. */
const C2_4PIG = (C_KMS * C_KMS) / (4 * Math.PI * G);

/** Angular-diameter distance to z, kpc (flat universe). */
export function angularDistance(z: number): number {
  return ((C_KMS / H0) * 1000 * comoving(z)) / (1 + z);
}

/** Critical surface density Σ_cr = c² D_s / (4π G D_l D_ls), M☉/kpc². Infinite if zs ≤ zl. */
export function criticalSurfaceDensity(zl: number, zs: number): number {
  const eff = efficiency(zl, zs); // D_ls / D_s
  return eff > 0 ? C2_4PIG / (angularDistance(zl) * eff) : Infinity;
}

/** Convergence κ(R) and mean convergence κ̄(<R) = M_2D / (π R² Σ_cr). */
export function convergence(proj: Projection, sigmaCr: number): { kappa: Float64Array; kappaBar: Float64Array } {
  const n = proj.R.length;
  const kappa = new Float64Array(n);
  const kappaBar = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    kappa[i] = proj.sigma[i] / sigmaCr;
    kappaBar[i] = proj.mass[i] / (Math.PI * proj.R[i] ** 2 * sigmaCr);
  }
  return { kappa, kappaBar };
}

/**
 * Einstein radius (kpc) of a halo embedded in a uniform external convergence κ_c:
 * the outermost R with κ̄_halo(<R) = 1 − κ_c, interpolated in log-log on the grid.
 * Zero if κ̄ never reaches the threshold at or above the grid's smallest radius.
 */
export function einsteinRadius(R: ArrayLike<number>, kappaBar: ArrayLike<number>, kappaC: number): number {
  const target = 1 - kappaC;
  for (let i = R.length - 2; i >= 0; i--) {
    if (kappaBar[i] >= target && kappaBar[i + 1] < target) {
      const w = Math.log(kappaBar[i] / target) / Math.log(kappaBar[i] / kappaBar[i + 1]);
      return Math.exp((1 - w) * Math.log(R[i]) + w * Math.log(R[i + 1]));
    }
  }
  return 0;
}

// --- everything the page needs, at one time or over all times ---------------------------------

/** The projected-radius grid used for lensing, kpc. */
export const LENS_GRID = logspace(1e-3, 100, 150);

export interface State {
  tau: number;
  profile: Profile;
  table: MassTable;
  vmax: number;
  rmax: number;
  /** κ̄(<R) on LENS_GRID. */
  kappaBar: Float64Array;
  /** Einstein radius, kpc (0 if none). */
  thetaE: number;
}

/** The halo at dimensionless time τ, with its derived quantities. */
export function stateAt(h: Halo, tau: number, sigmaCr: number, kappaC: number): State {
  const profile = profileAt(h, tau);
  const rt = TRUNCATION * h.r200;
  const table = massTable(profile, rt);
  const { vmax: v, rmax } = vmax(table);
  const { kappaBar } = convergence(project(profile, rt, LENS_GRID, table), sigmaCr);
  return { tau, profile, table, vmax: v, rmax, kappaBar, thetaE: einsteinRadius(LENS_GRID, kappaBar, kappaC) };
}

/** Central density and core radius on a grid of times 0…span (Gyr); NaN past TAU_MAX, where the model ends. */
export function history(h: Halo, tc: number, span: number, n = 400): { t: Float64Array; rho0: Float64Array; rC: Float64Array } {
  const t = new Float64Array(n);
  const rho0 = new Float64Array(n);
  const rC = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    t[i] = (i / (n - 1)) * span;
    const tau = t[i] / tc;
    if (tau > TAU_MAX) {
      rho0[i] = rC[i] = NaN;
      continue;
    }
    const p = profileAt(h, tau);
    rho0[i] = centralDensity(p);
    rC[i] = p.rC;
  }
  return { t, rho0, rC };
}
