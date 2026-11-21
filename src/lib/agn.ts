// A toy AGN-variability Hubble diagram: simulated quasar light curves, a toy
// relation between variability and luminosity, a maximum-likelihood fit of
// each light curve, and the distances that follow. Pure TypeScript, no DOM, so
// it runs on the page, in a Web Worker (agn-worker.ts) and in node. Used by
// src/pages/projects/agn-hubble.astro.
//
// Real physics:
//  - Light curves are a damped random walk (an Ornstein–Uhlenbeck process; Kelly
//    et al. 2009; MacLeod et al. 2010) in magnitudes with long-term rms σ and
//    rest-frame damping time τ: covariance σ² exp(−|Δt_rest|/τ), with
//    Δt_rest = Δt_obs / (1 + z). σ is the asymptotic standard deviation; the
//    asymptotic structure function is SF∞ = √2 σ (MacLeod et al. 2010). Kelly
//    et al. (2009) write the same process with a driving amplitude σ̂, related
//    by σ̂ = σ √(2/τ), i.e. long-term variance σ² = σ̂² τ / 2.
//    It is sampled exactly at any times:
//      x(t + Δt) = x(t) e^{−Δt/τ} + N(0, σ² [1 − e^{−2Δt/τ}]).
//  - Distances come from the CPL cosmology in cpl.ts.
//
// Toy choices (stated on the page):
//  - log10 σ = a_σ + b_σ (M + 23) + scatter, with σ = 0.2/√2 = 0.141 mag at
//    M = −23 and b_σ = 0.090 dex/mag, the luminosity coefficient MacLeod et al.
//    (2010, Table 1) fit for SF∞ when black hole mass is left out of the model
//    (it is 0.131 at fixed black hole mass). σ ∝ SF∞, so the slope is the same.
//  - log10 τ = a_τ + b_τ (M + 23) + scatter, with τ = 200 days (rest frame) at
//    M = −23 and b_τ = −0.08 dex/mag by default. This slope is invented so that
//    τ carries distance information; MacLeod et al. find τ nearly independent
//    of luminosity, which b_τ = 0 reproduces.
//  - Redshifts uniform in 0.2 to 4; absolute magnitudes uniform between −28 and
//    the brighter of −22 and the flux limit m = 23 (so brighter at high z). The
//    flux limit is applied with distances from a fixed reference ΛCDM model
//    (REFERENCE), so the AGN themselves, and so their fits, do not change when
//    the generating cosmology does.
//  - Uniform cadence with no seasonal gaps; Gaussian photometric errors of one
//    known size; no K-corrections, extinction, host galaxy light or selection
//    effects beyond the flux limit.
//
// Measurement (fitDRW): the exact Gaussian likelihood of each light curve under
// a damped random walk plus white noise with an unknown mean, evaluated in O(n)
// by a Kalman filter, maximised over (ln σ, ln τ_obs); the curvature at the
// maximum gives the errors. Luminosity (luminosity()): a generalised least
// squares combination of the two fitted parameters through the two relations.
import { LCDM, Universe, comovingSorted, withLambda, type CPLParams } from './cpl';

// --- random numbers -------------------------------------------------------------

/** Mixes integers into one 32-bit seed, so each object gets its own reproducible stream. */
export function hashSeed(...values: number[]): number {
  let h = 0x9e3779b9;
  for (const v of values) {
    h = Math.imul(h ^ (v | 0), 0x85ebca6b);
    h ^= h >>> 13;
    h = Math.imul(h, 0xc2b2ae35);
    h ^= h >>> 16;
  }
  return h >>> 0;
}

/** mulberry32 uniform deviates and Box–Muller normal deviates from a 32-bit seed. */
export class Random {
  private state: number;
  private spare = NaN;

  constructor(seed: number) {
    this.state = seed >>> 0;
  }

  /** Uniform on [0, 1). */
  uniform(): number {
    let t = (this.state = (this.state + 0x6d2b79f5) | 0);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  /** Standard normal. */
  normal(): number {
    if (!Number.isNaN(this.spare)) {
      const s = this.spare;
      this.spare = NaN;
      return s;
    }
    const u = 1 - this.uniform(); // (0, 1], so the log is finite
    const r = Math.sqrt(-2 * Math.log(u));
    const phi = 2 * Math.PI * this.uniform();
    this.spare = r * Math.sin(phi);
    return r * Math.cos(phi);
  }
}

// --- the damped random walk --------------------------------------------------------

/**
 * A damped random walk at the given observed times (days, ascending), in
 * magnitudes about zero, with long-term rms σ (mag) and rest-frame damping
 * time τ (days) for a source at redshift z. Starts from the stationary
 * distribution, so every epoch has variance σ².
 */
export function dampedRandomWalk(times: ArrayLike<number>, sigma: number, tauRest: number, z: number, rng: Random): Float64Array {
  const tauObs = tauRest * (1 + z); // time dilation: the observed clock runs (1 + z) times slower
  const x = new Float64Array(times.length);
  if (!times.length) return x;
  x[0] = sigma * rng.normal();
  for (let k = 1; k < times.length; k++) {
    const decay = Math.exp(-(times[k] - times[k - 1]) / tauObs);
    x[k] = x[k - 1] * decay + sigma * Math.sqrt(1 - decay * decay) * rng.normal();
  }
  return x;
}

// --- the likelihood ------------------------------------------------------------------

export interface Likelihood {
  /** ln L with the mean marginalised over a flat prior (see drwLikelihood). */
  lnL: number;
  /** Best estimate of the light curve's mean, mag, and its variance. */
  mean: number;
  meanVar: number;
}

const LN_2PI = Math.log(2 * Math.PI);

/**
 * Log-likelihood of magnitudes y at observed times t (days, ascending) for
 *   y_i = μ + x(t_i) + ε_i,  x a damped random walk with long-term rms σ and
 *   observed damping time τ,  ε_i ~ N(0, noise2_i) independent,
 * with the unknown mean μ treated as a diffuse parameter: integrated out
 * against a flat prior (the restricted, or REML, likelihood). With C the
 * covariance of x + ε and 1 the vector of ones,
 *
 *   ln L = −½ [ (n − 1) ln 2π + ln det C + ln(1ᵀC⁻¹1) + yᵀC⁻¹y − (1ᵀC⁻¹y)² / 1ᵀC⁻¹1 ].
 *
 * Integrating the mean out instead of maximising over it removes most of the
 * low bias in σ that comes from the mean absorbing slow wander.
 *
 * It is computed in O(n) by a Kalman filter on the state x: predict
 * x_{i+1} = a x_i + w with a = exp(−Δt/τ) and Var w = σ²(1 − a²), starting
 * from the stationary N(0, σ²); update with y_i. The filter is linear and its
 * gains do not depend on the data, so running it on y and on 1 at once gives
 * the innovations of y − μ1 for every μ: v_y − μ v_1, both with variance S_i.
 * Then ln det C = Σ ln S_i, 1ᵀC⁻¹1 = Σ v_1²/S, 1ᵀC⁻¹y = Σ v_1 v_y/S, yᵀC⁻¹y = Σ v_y²/S.
 */
export function drwLikelihood(t: ArrayLike<number>, y: ArrayLike<number>, noise2: number | ArrayLike<number>, sigma: number, tauObs: number): Likelihood {
  const n = t.length;
  const s2 = sigma * sigma;
  const scalarNoise = typeof noise2 === 'number';
  let mY = 0; // predicted state mean, run on y
  let m1 = 0; // and on the vector of ones
  let P = s2; // predicted state variance
  let logS = 0;
  let prodS = 1; // Σ ln S is accumulated as a product, with a log only when it nears the float range
  let yy = 0;
  let y1 = 0;
  let oo = 0;
  let lastDt = NaN;
  let a = 0;
  let q = 0;
  for (let i = 0; i < n; i++) {
    const e2 = scalarNoise ? noise2 : noise2[i];
    const S = P + e2;
    const vy = y[i] - mY;
    const v1 = 1 - m1;
    const inv = 1 / S;
    prodS *= S;
    if (prodS < 1e-200 || prodS > 1e200) {
      logS += Math.log(prodS);
      prodS = 1;
    }
    yy += vy * vy * inv;
    y1 += vy * v1 * inv;
    oo += v1 * v1 * inv;
    if (i === n - 1) break;
    const K = P * inv;
    const Pf = P * e2 * inv; // = P (1 − K), written so it is exactly zero with no noise
    const dt = t[i + 1] - t[i];
    if (dt !== lastDt) {
      // Uniform cadence reuses the same decay; only a new interval needs an exp.
      lastDt = dt;
      a = Math.exp(-dt / tauObs);
      q = -s2 * Math.expm1((-2 * dt) / tauObs); // σ² (1 − a²), accurate for Δt ≪ τ
    }
    mY = a * (mY + K * vy);
    m1 = a * (m1 + K * v1);
    P = a * a * Pf + q;
  }
  logS += Math.log(prodS);
  const lnL = -0.5 * ((n - 1) * LN_2PI + logS + Math.log(oo) + yy - (y1 * y1) / oo);
  return { lnL, mean: y1 / oo, meanVar: 1 / oo };
}

// --- the fit -------------------------------------------------------------------------

export interface DRWFit {
  /** Maximum-likelihood ln σ (σ in mag) and ln τ_obs (τ in days, observed frame). */
  lnSigma: number;
  lnTau: number;
  /** Covariance of (ln σ, ln τ_obs): [var σ, cov, var τ]. See fitDRW for the unconstrained case. */
  cov: [number, number, number];
  lnL: number;
  /** Fitted mean of the light curve (mag) and its variance. */
  mean: number;
  meanVar: number;
  /** −1 or +1 if the maximum is on the lower or upper bound of τ_obs, else 0. */
  tauBound: -1 | 0 | 1;
  /** True if the maximum is on a bound or the curvature there is not negative definite. */
  unconstrained: boolean;
  /** Likelihood evaluations used. */
  evaluations: number;
}

export interface FitOptions {
  /** Allowed τ_obs, days. */
  tauRange: [number, number];
  /** Allowed σ, mag. */
  sigmaRange?: [number, number];
}

/** The default τ_obs range for a light curve: half the shortest interval to ten times the baseline. */
export function defaultTauRange(t: ArrayLike<number>): [number, number] {
  let minDt = Infinity;
  for (let i = 1; i < t.length; i++) minDt = Math.min(minDt, t[i] - t[i - 1]);
  return [minDt / 2, 10 * (t[t.length - 1] - t[0])];
}

const STEP = 0.01; // finite-difference step in ln σ and ln τ
const GRID_TAU = 12;
const GOLD = 0.3819660112501051; // (3 − √5) / 2

/**
 * Brent's method (golden section with parabolic steps) for the minimum of fn on
 * [a, b], started from x with fn(x) = fx known. fn also returns a ln σ, kept
 * for the best point, so the profile can be warm-started.
 */
function brent(fn: (x: number) => { ls: number; value: number }, a: number, b: number, x: number, fx: number, lsx: number, tol: number) {
  let w = x;
  let v = x;
  let fw = fx;
  let fv = fx;
  let ls = lsx;
  let d = 0;
  let e = 0;
  for (let iter = 0; iter < 60; iter++) {
    const m = 0.5 * (a + b);
    const tol1 = tol + 1e-10 * Math.abs(x);
    if (Math.abs(x - m) <= 2 * tol1 - 0.5 * (b - a)) break;
    let useGolden = true;
    if (Math.abs(e) > tol1) {
      // Parabola through x, w, v.
      const r = (x - w) * (fx - fv);
      let q = (x - v) * (fx - fw);
      let p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) p = -p;
      q = Math.abs(q);
      if (Math.abs(p) < Math.abs(0.5 * q * e) && p > q * (a - x) && p < q * (b - x)) {
        e = d;
        d = p / q;
        const u = x + d;
        if (u - a < 2 * tol1 || b - u < 2 * tol1) d = x < m ? tol1 : -tol1;
        useGolden = false;
      }
    }
    if (useGolden) {
      e = (x < m ? b : a) - x;
      d = GOLD * e;
    }
    const u = Math.abs(d) >= tol1 ? x + d : x + (d > 0 ? tol1 : -tol1);
    const r = fn(u);
    const fu = r.value;
    if (fu <= fx) {
      if (u < x) b = x;
      else a = x;
      [v, fv, w, fw, x, fx, ls] = [w, fw, x, fx, u, fu, r.ls];
    } else {
      if (u < x) a = u;
      else b = u;
      if (fu <= fw || w === x) {
        [v, fv, w, fw] = [w, fw, u, fu];
      } else if (fu <= fv || v === x || v === w) {
        [v, fv] = [u, fu];
      }
    }
  }
  return { x, value: fx, ls };
}

/**
 * Maximum-likelihood σ and τ_obs for one light curve (drwLikelihood), over
 * p = (ln σ, ln τ_obs) inside the given bounds:
 *
 *  1. the profile likelihood max_σ L(σ, τ) on a grid of 12 values of ln τ
 *     evenly over its range (each a short Newton search in ln σ, started on the
 *     σ that would give the light curve's excess variance s² − e² for that τ),
 *     then Brent's method in ln τ on the bracket around the best grid value;
 *  2. Newton steps in both coordinates from there, with the gradient and Hessian from
 *     central differences (step 0.01), each step capped at 1 in either
 *     coordinate and halved until the likelihood rises; where the Hessian is not
 *     negative definite, minus the Hessian is shifted until it is positive
 *     definite (Levenberg–Marquardt). On a τ bound with the gradient pointing
 *     out, only σ moves.
 *  3. The covariance is the inverse of minus the Hessian at the maximum.
 *
 * When the maximum is on a τ bound, or the Hessian there is not negative
 * definite, τ is not measured and the fit is flagged `unconstrained`. Its
 * covariance is then made honest rather than small: ln τ gets the variance of a
 * uniform distribution over the whole allowed range, (range)² / 12, and ln σ
 * inherits the part of it that the likelihood ridge ties to σ. Along the ridge
 * ln σ moves by s = −H_στ / H_σσ per unit of ln τ (s ≈ ½ when τ is far longer
 * than the light curve, where only the driving amplitude σ̂² = 2σ²/τ is
 * measured), so var ln σ = −1/H_σσ + s² V_τ and cov = s V_τ.
 */
export function fitDRW(t: ArrayLike<number>, y: ArrayLike<number>, noise2: number | ArrayLike<number>, options: FitOptions): DRWFit {
  const n = t.length;
  const [tauMin, tauMax] = options.tauRange;
  const [sigMin, sigMax] = options.sigmaRange ?? [1e-4, 30];
  const lo = [Math.log(sigMin), Math.log(tauMin)];
  const hi = [Math.log(sigMax), Math.log(tauMax)];
  let evaluations = 0;
  const f = (ls: number, lt: number) => {
    evaluations++;
    return drwLikelihood(t, y, noise2, Math.exp(ls), Math.exp(lt)).lnL;
  };

  // Starting scale: the excess variance over the noise, s_x² = s² − e².
  let mean = 0;
  for (let i = 0; i < n; i++) mean += y[i];
  mean /= n;
  let s2 = 0;
  let e2 = 0;
  for (let i = 0; i < n; i++) {
    s2 += (y[i] - mean) ** 2;
    e2 += typeof noise2 === 'number' ? noise2 : noise2[i];
  }
  s2 /= Math.max(1, n - 1);
  e2 /= n;
  const excess = Math.max(s2 - e2, 0.09 * s2, 1e-12);
  const baseline = Math.max(t[n - 1] - t[0], 1e-9);

  // The likelihood is a long curved ridge in (ln σ, ln τ): a
  // light curve shorter than τ only shows a fraction of σ². So for each τ the σ
  // search starts on the σ that would give the observed excess variance,
  // E[s²] = σ² g(τ/T) with g(x) = 1 − 2x + 2x² (1 − e^{−1/x}) for a baseline T
  // (the expected mean-subtracted variance of a damped random walk; g → 1/(3x)
  // for x ≫ 1).
  const ridge = (lt: number) => {
    const x = Math.exp(lt) / baseline;
    const g = x < 1e-3 ? 1 - 2 * x : 1 - 2 * x - 2 * x * x * Math.expm1(-1 / x);
    return 0.5 * Math.log(excess / Math.max(g, 1e-12));
  };
  // The profile likelihood in τ: for each ln τ, maximise over ln σ by Newton
  // steps with central differences, starting from a given ln σ. In σ alone the
  // likelihood is smooth and single-peaked, so this is quick and safe.
  const profile = (lt: number, ls0: number): [number, number] => {
    let ls = Math.min(hi[0], Math.max(lo[0], ls0));
    let fl = f(ls, lt);
    for (let iter = 0; iter < 20; iter++) {
      const h = STEP;
      const f1 = f(ls + h, lt);
      const f0 = f(ls - h, lt);
      const g = (f1 - f0) / (2 * h);
      const H = (f1 - 2 * fl + f0) / (h * h);
      let d = H < 0 ? -g / H : Math.sign(g) * 0.5;
      d = Math.max(-1, Math.min(1, d));
      let moved = false;
      for (let k = 0; k < 8; k++, d /= 2) {
        const q = Math.min(hi[0], Math.max(lo[0], ls + d));
        const fq = f(q, lt);
        if (fq > fl) {
          moved = Math.abs(q - ls) > 1e-4;
          ls = q;
          fl = fq;
          break;
        }
      }
      if (!moved) break;
    }
    return [ls, fl];
  };

  // 1. The profile on a grid of 12 values of ln τ, then Brent's method on the
  // bracket around the best one.
  const grid: { lt: number; ls: number; fl: number }[] = [];
  for (let j = 0; j < GRID_TAU; j++) {
    const lt = lo[1] + ((j + 0.5) * (hi[1] - lo[1])) / GRID_TAU;
    const [ls, fl] = profile(lt, ridge(lt));
    grid.push({ lt, ls, fl });
  }
  let jb = 0;
  for (let j = 1; j < GRID_TAU; j++) if (grid[j].fl > grid[jb].fl) jb = j;
  const left = jb > 0 ? grid[jb - 1].lt : lo[1];
  const right = jb < GRID_TAU - 1 ? grid[jb + 1].lt : hi[1];
  const warm = grid[jb].ls;
  const negProfile = (lt: number) => {
    const [ls, fl] = profile(lt, warm + (ridge(lt) - ridge(grid[jb].lt)));
    return { ls, value: -fl };
  };
  const found = brent(negProfile, left, right, grid[jb].lt, -grid[jb].fl, grid[jb].ls, 2e-3);
  const best = [found.ls, found.x];
  const fBest = -found.value;

  // 2. Newton refinement.
  const clamp = (v: number, k: number) => Math.min(hi[k], Math.max(lo[k], v));
  const stencil = (p: number[], fp: number) => {
    const [a, b] = p;
    const h = STEP;
    const fa1 = f(a + h, b);
    const fa0 = f(a - h, b);
    const fb1 = f(a, b + h);
    const fb0 = f(a, b - h);
    const fpp = f(a + h, b + h);
    const fpm = f(a + h, b - h);
    const fmp = f(a - h, b + h);
    const fmm = f(a - h, b - h);
    return {
      g: [(fa1 - fa0) / (2 * h), (fb1 - fb0) / (2 * h)],
      H: [(fa1 - 2 * fp + fa0) / (h * h), (fpp - fpm - fmp + fmm) / (4 * h * h), (fb1 - 2 * fp + fb0) / (h * h)],
    };
  };
  const outward = (p: number[], g: number[]): -1 | 0 | 1 => {
    if (p[1] >= hi[1] - 1e-6 && g[1] > 0) return 1;
    if (p[1] <= lo[1] + 1e-6 && g[1] < 0) return -1;
    return 0;
  };

  let p = best;
  let fp = fBest;
  let d2 = stencil(p, fp);
  for (let iter = 0; iter < 40; iter++) {
    const { g, H } = d2;
    const bound = outward(p, g);
    let d: number[];
    const det = H[0] * H[2] - H[1] * H[1];
    if (bound !== 0) d = [H[0] < 0 ? -g[0] / H[0] : g[0] / Math.max(1, Math.abs(H[0])), 0];
    else if (H[0] < 0 && det > 0) d = [(-H[2] * g[0] + H[1] * g[1]) / det, (H[1] * g[0] - H[0] * g[1]) / det];
    else {
      // Not a maximum yet: shift minus the Hessian until it is positive definite
      // (a Levenberg–Marquardt step), which follows a curved ridge far better
      // than a plain gradient step.
      const [a00, a01, a11] = [-H[0], -H[1], -H[2]];
      const lmin = 0.5 * (a00 + a11 - Math.hypot(a00 - a11, 2 * a01));
      const shift = -lmin + 1e-3 * Math.max(Math.abs(a00), Math.abs(a11), 1);
      const [b00, b11] = [a00 + shift, a11 + shift];
      const dt = b00 * b11 - a01 * a01;
      d = [(b11 * g[0] - a01 * g[1]) / dt, (b00 * g[1] - a01 * g[0]) / dt];
    }
    d = d.map((v) => Math.max(-1, Math.min(1, v)));
    let step = 1;
    let moved = false;
    for (let k = 0; k < 12; k++, step /= 2) {
      const q = [clamp(p[0] + step * d[0], 0), clamp(p[1] + step * d[1], 1)];
      const fq = f(q[0], q[1]);
      if (fq > fp) {
        const shift = Math.max(Math.abs(q[0] - p[0]), Math.abs(q[1] - p[1]));
        const gain = fq - fp;
        p = q;
        fp = fq;
        moved = shift > 1e-6 && gain > 1e-10;
        break;
      }
    }
    d2 = stencil(p, fp);
    if (!moved) break;
  }

  // 3. Errors from the curvature.
  const { g, H } = d2;
  const tauBound = outward(p, g);
  const det = H[0] * H[2] - H[1] * H[1];
  const negDef = H[0] < 0 && det > 0;
  const atSigmaBound = p[0] <= lo[0] + 1e-6 || p[0] >= hi[0] - 1e-6;
  const unconstrained = tauBound !== 0 || !negDef || atSigmaBound;
  let cov: [number, number, number];
  if (!unconstrained) {
    cov = [-H[2] / det, H[1] / det, -H[0] / det];
  } else {
    const vTau = (hi[1] - lo[1]) ** 2 / 12;
    const condSigma = H[0] < 0 ? -1 / H[0] : 1;
    const slope = H[0] < 0 ? Math.max(-1, Math.min(1, -H[1] / H[0])) : tauBound > 0 ? 0.5 : 0;
    cov = [condSigma + slope * slope * vTau, slope * vTau, vTau];
  }
  const L = drwLikelihood(t, y, noise2, Math.exp(p[0]), Math.exp(p[1]));
  return { lnSigma: p[0], lnTau: p[1], cov, lnL: L.lnL, mean: L.mean, meanVar: L.meanVar, tauBound, unconstrained, evaluations };
}

// --- the toy relations -------------------------------------------------------------

export interface Relations {
  /** log10 σ (mag) at the pivot, and its slope d log10 σ / dM, dex per mag. */
  aSigma: number;
  bSigma: number;
  /** log10 τ_rest (days) at the pivot, and its slope, dex per mag. */
  aTau: number;
  bTau: number;
  /** Intrinsic Gaussian scatter about each relation, dex. */
  scatterSigma: number;
  scatterTau: number;
  pivot: number;
}

export const PIVOT = -23;
/** σ = 0.2/√2 mag (SF∞ = 0.2 mag) at M = −23. */
export const A_SIGMA = Math.log10(0.2 / Math.SQRT2);
/** MacLeod et al. (2010), Table 1: SF∞ against M_i with black hole mass left out. */
export const B_SIGMA = 0.09;
/** τ = 200 days, rest frame, at M = −23. */
export const A_TAU = Math.log10(200);
/** The toy default for the τ slope (more luminous AGN flicker more slowly). */
export const B_TAU_DEFAULT = -0.08;

export function relations(p: { bTau: number; scatterSigma: number; scatterTau: number }): Relations {
  return { aSigma: A_SIGMA, bSigma: B_SIGMA, aTau: A_TAU, bTau: p.bTau, scatterSigma: p.scatterSigma, scatterTau: p.scatterTau, pivot: PIVOT };
}

/** σ (mag) and rest-frame τ (days) on the relations, without scatter. */
export const relationSigma = (M: number, rel: Relations): number => 10 ** (rel.aSigma + rel.bSigma * (M - rel.pivot));
export const relationTau = (M: number, rel: Relations): number => 10 ** (rel.aTau + rel.bTau * (M - rel.pivot));

export interface Luminosity {
  /** Combined generalised-least-squares estimate of M and its variance, mag². */
  M: number;
  varM: number;
  /** From σ alone and from τ alone (NaN when b_τ = 0). */
  MSigma: number;
  varSigma: number;
  MTau: number;
  varTau: number;
}

/**
 * Absolute magnitude from the fitted log10 σ and log10 τ_rest with covariance
 * [var σ, cov, var τ] (dex²). With d = (log10 σ̂ − a_σ, log10 τ̂ − a_τ),
 * b = (b_σ, b_τ) and W = (C_meas + S_intrinsic)⁻¹, S = diag(scatter_σ², scatter_τ²),
 * the model is d = b (M − pivot) + noise of covariance W⁻¹, whose generalised
 * least squares solution is
 *   M̂ = pivot + bᵀWd / bᵀWb,   Var M̂ = 1 / bᵀWb.
 * The single-indicator estimates are d_σ / b_σ and d_τ / b_τ with variances
 * (C + S)_σσ / b_σ² and (C + S)_ττ / b_τ².
 */
export function luminosity(logSigma: number, logTau: number, cov: ArrayLike<number>, rel: Relations): Luminosity {
  const d0 = logSigma - rel.aSigma;
  const d1 = logTau - rel.aTau;
  const A00 = cov[0] + rel.scatterSigma ** 2;
  const A01 = cov[1];
  const A11 = cov[2] + rel.scatterTau ** 2;
  const det = A00 * A11 - A01 * A01;
  const [W00, W01, W11] = [A11 / det, -A01 / det, A00 / det];
  const b0 = rel.bSigma;
  const b1 = rel.bTau;
  const bWb = b0 * (W00 * b0 + W01 * b1) + b1 * (W01 * b0 + W11 * b1);
  const bWd = b0 * (W00 * d0 + W01 * d1) + b1 * (W01 * d0 + W11 * d1);
  const tauUseful = Math.abs(b1) > 1e-9;
  return {
    M: rel.pivot + bWd / bWb,
    varM: 1 / bWb,
    MSigma: rel.pivot + d0 / b0,
    varSigma: A00 / (b0 * b0),
    MTau: tauUseful ? rel.pivot + d1 / b1 : NaN,
    varTau: tauUseful ? A11 / (b1 * b1) : NaN,
  };
}

// --- the mock survey ---------------------------------------------------------------

export const Z_MIN = 0.2;
export const Z_MAX = 4;
export const M_BRIGHT = -28;
export const M_FAINT = -22;
/** Flux limit, apparent magnitude. */
export const M_LIMIT = 23;
/** The fixed cosmology that applies the flux limit (see the header). */
export const REFERENCE: CPLParams = LCDM;

/** Everything that changes the light curves or their fits. The cosmology is not here: changing it never refits. */
export interface SurveyParams {
  count: number;
  years: number;
  /** Days between epochs, observed frame. */
  cadence: number;
  /** Photometric error per epoch, mag. */
  photErr: number;
  scatterSigma: number;
  scatterTau: number;
  bTau: number;
  seed: number;
}

/** The simulated AGN, before any fit. Independent of the generating cosmology. */
export interface Population {
  params: SurveyParams;
  rel: Relations;
  /** Observation times, days. */
  times: Float64Array;
  z: Float64Array;
  M: Float64Array;
  /** True log10 σ (mag) and log10 τ_rest (days). */
  logSigma: Float64Array;
  logTau: Float64Array;
  /** True observed-frame damping time, days. */
  tauObs: Float64Array;
}

/** Observation times in days for a survey: 0, cadence, 2 cadence, … up to the length. */
export function surveyTimes(params: Pick<SurveyParams, 'years' | 'cadence'>): Float64Array {
  const n = Math.floor((params.years * 365.25) / params.cadence) + 1;
  return Float64Array.from({ length: n }, (_, k) => k * params.cadence);
}

const populationRandom = (seed: number, i: number) => new Random(hashSeed(seed, i, 1));
const curveRandom = (seed: number, i: number) => new Random(hashSeed(seed, i, 2));

/** Draws redshifts, luminosities and the true σ and τ of every AGN. Cheap. */
export function makePopulation(params: SurveyParams): Population {
  const { count, seed } = params;
  const rel = relations(params);
  const reference = new Universe(REFERENCE, Z_MAX + 0.1);
  const pop: Population = {
    params: { ...params },
    rel,
    times: surveyTimes(params),
    z: new Float64Array(count),
    M: new Float64Array(count),
    logSigma: new Float64Array(count),
    logTau: new Float64Array(count),
    tauObs: new Float64Array(count),
  };
  for (let i = 0; i < count; i++) {
    const r = populationRandom(seed, i);
    const z = Z_MIN + (Z_MAX - Z_MIN) * r.uniform();
    // Uniform between −28 and the brighter of −22 and the flux limit.
    const faint = Math.max(M_BRIGHT + 0.5, Math.min(M_FAINT, M_LIMIT - reference.distanceModulus(z)));
    const M = M_BRIGHT + r.uniform() * (faint - M_BRIGHT);
    pop.z[i] = z;
    pop.M[i] = M;
    pop.logSigma[i] = rel.aSigma + rel.bSigma * (M - rel.pivot) + rel.scatterSigma * r.normal();
    pop.logTau[i] = rel.aTau + rel.bTau * (M - rel.pivot) + rel.scatterTau * r.normal();
    pop.tauObs[i] = 10 ** pop.logTau[i] * (1 + z);
  }
  return pop;
}

/**
 * Object i's light curve as observed, in magnitudes about its true mean
 * apparent magnitude: σ u + e v, with u a unit damped random walk and v unit
 * white noise. Exact for the uniform cadence; reproducible from the seed.
 */
export function lightCurve(pop: Population, i: number): Float64Array {
  const t = pop.times;
  const rng = curveRandom(pop.params.seed, i);
  const sigma = 10 ** pop.logSigma[i];
  const e = pop.params.photErr;
  const y = new Float64Array(t.length);
  if (!t.length) return y;
  const rho = Math.exp(-pop.params.cadence / pop.tauObs[i]);
  const kick = Math.sqrt(1 - rho * rho);
  let u = rng.normal();
  y[0] = sigma * u + e * rng.normal();
  for (let k = 1; k < t.length; k++) {
    u = rho * u + kick * rng.normal();
    y[k] = sigma * u + e * rng.normal();
  }
  return y;
}

// --- fitting a whole survey -----------------------------------------------------------

/** Columns of one fitted object, in the order they are packed (FIELDS numbers per object). */
export const FIELDS = 9;
export const F = {
  /** Fitted log10 σ (mag) and log10 τ_rest (days). */
  LogSigma: 0,
  LogTau: 1,
  /** Their covariance, dex². */
  VarSigma: 2,
  Cov: 3,
  VarTau: 4,
  /** Fitted light-curve mean minus the true mean apparent magnitude, mag, and its variance. */
  Mean: 5,
  MeanVar: 6,
  /** 1 if τ is unconstrained. */
  Unconstrained: 7,
  /** −1, 0, +1: on the lower or upper τ bound. */
  TauBound: 8,
} as const;

/** Fits object i and writes its FIELDS numbers into out at offset. */
export function fitObject(pop: Population, i: number, out: Float64Array, offset = 0): DRWFit {
  const y = lightCurve(pop, i);
  const e = Math.max(pop.params.photErr, 1e-6);
  const fit = fitDRW(pop.times, y, e * e, { tauRange: defaultTauRange(pop.times) });
  const ln10 = Math.LN10;
  out[offset + F.LogSigma] = fit.lnSigma / ln10;
  out[offset + F.LogTau] = (fit.lnTau - Math.log(1 + pop.z[i])) / ln10;
  out[offset + F.VarSigma] = fit.cov[0] / (ln10 * ln10);
  out[offset + F.Cov] = fit.cov[1] / (ln10 * ln10);
  out[offset + F.VarTau] = fit.cov[2] / (ln10 * ln10);
  out[offset + F.Mean] = fit.mean;
  out[offset + F.MeanVar] = fit.meanVar;
  out[offset + F.Unconstrained] = fit.unconstrained ? 1 : 0;
  out[offset + F.TauBound] = fit.tauBound;
  return fit;
}

/**
 * Fits every object, a slice at a time: each slice runs until `sliceMs` has
 * passed, hands its rows to `emit`, then yields so the caller can cancel
 * (`stale()` returning true stops the job). Used by the Web Worker and by the
 * page's fallback on the main thread.
 */
export async function fitPopulation(
  pop: Population,
  emit: (from: number, rows: Float64Array) => void,
  stale: () => boolean,
  sliceMs = 30,
): Promise<boolean> {
  const count = pop.z.length;
  let i = 0;
  while (i < count) {
    if (stale()) return false;
    const start = performance.now();
    const from = i;
    const rows: number[] = [];
    const row = new Float64Array(FIELDS);
    do {
      fitObject(pop, i, row);
      for (let k = 0; k < FIELDS; k++) rows.push(row[k]);
      i++;
    } while (i < count && performance.now() - start < sliceMs);
    emit(from, Float64Array.from(rows));
    await new Promise((resolve) => setTimeout(resolve, 0));
  }
  return !stale();
}

// --- distances -------------------------------------------------------------------------

export interface Measurement {
  /** Objects with a fit so far. */
  fitted: number;
  /** 1 where object i has been fitted. */
  done: Uint8Array;
  muTrue: Float64Array;
  /** ΛCDM with the same Ωm and H0, for residuals. */
  muLCDM: Float64Array;
  /** Measured mean apparent magnitude, mag. */
  mObs: Float64Array;
  /** Combined estimate: absolute magnitude, distance modulus and its error (mag). */
  MEst: Float64Array;
  muEst: Float64Array;
  muErr: Float64Array;
  /** σ only and τ only. */
  muSigma: Float64Array;
  muErrSigma: Float64Array;
  muTau: Float64Array;
  muErrTau: Float64Array;
  unconstrained: number;
  /** Mean and standard deviation of μ̂ − μ_true, combined and single-indicator, mag. */
  bias: number;
  spread: number;
  spreadSigma: number;
  spreadTau: number;
  /** Median of survey length / true τ_obs over fitted objects. */
  medianBaseline: number;
}

/**
 * Distances for the fitted objects in a given generating cosmology:
 * m̂ = M + μ_true + (fitted light-curve mean offset), M̂ from luminosity(),
 * μ̂ = m̂ − M̂, Var μ̂ = Var M̂ + Var m̂. Nothing here refits.
 */
export function measure(pop: Population, fits: Float64Array, done: Uint8Array, cosmology: CPLParams): Measurement {
  const count = pop.z.length;
  const truth = new Universe(cosmology, Z_MAX + 0.1);
  const lambda = new Universe(withLambda(cosmology), Z_MAX + 0.1);
  const arr = () => new Float64Array(count).fill(NaN);
  const m: Measurement = {
    fitted: 0,
    done,
    muTrue: arr(),
    muLCDM: arr(),
    mObs: arr(),
    MEst: arr(),
    muEst: arr(),
    muErr: arr(),
    muSigma: arr(),
    muErrSigma: arr(),
    muTau: arr(),
    muErrTau: arr(),
    unconstrained: 0,
    bias: NaN,
    spread: NaN,
    spreadSigma: NaN,
    spreadTau: NaN,
    medianBaseline: NaN,
  };
  const length = pop.times[pop.times.length - 1] - pop.times[0];
  const ratios: number[] = [];
  const stats = [0, 0, 0, 0, 0, 0, 0];
  for (let i = 0; i < count; i++) {
    m.muTrue[i] = truth.distanceModulus(pop.z[i]);
    m.muLCDM[i] = lambda.distanceModulus(pop.z[i]);
    if (!done[i]) continue;
    const o = i * FIELDS;
    const lum = luminosity(fits[o + F.LogSigma], fits[o + F.LogTau], fits.subarray(o + F.VarSigma, o + F.VarTau + 1), pop.rel);
    const mObs = pop.M[i] + m.muTrue[i] + fits[o + F.Mean];
    const varm = fits[o + F.MeanVar];
    m.mObs[i] = mObs;
    m.MEst[i] = lum.M;
    m.muEst[i] = mObs - lum.M;
    m.muErr[i] = Math.sqrt(lum.varM + varm);
    m.muSigma[i] = mObs - lum.MSigma;
    m.muErrSigma[i] = Math.sqrt(lum.varSigma + varm);
    m.muTau[i] = mObs - lum.MTau;
    m.muErrTau[i] = Math.sqrt(lum.varTau + varm);
    m.fitted++;
    if (fits[o + F.Unconstrained]) m.unconstrained++;
    ratios.push(length / pop.tauObs[i]);
    const d = m.muEst[i] - m.muTrue[i];
    const ds = m.muSigma[i] - m.muTrue[i];
    const dt = m.muTau[i] - m.muTrue[i];
    stats[0] += d;
    stats[1] += d * d;
    stats[2] += ds;
    stats[3] += ds * ds;
    stats[4] += dt;
    stats[5] += dt * dt;
  }
  const n = m.fitted;
  if (n > 0) {
    const sd = (s: number, s2: number) => Math.sqrt(Math.max(0, s2 / n - (s / n) ** 2));
    m.bias = stats[0] / n;
    m.spread = sd(stats[0], stats[1]);
    m.spreadSigma = sd(stats[2], stats[3]);
    m.spreadTau = sd(stats[4], stats[5]);
    ratios.sort((a, b) => a - b);
    m.medianBaseline = n % 2 ? ratios[(n - 1) / 2] : 0.5 * (ratios[n / 2 - 1] + ratios[n / 2]);
  }
  return m;
}

// --- summaries -------------------------------------------------------------------------

export interface Bins {
  z: number[];
  mean: number[];
  /** Standard error of the weighted mean, from the weighted scatter inside the bin. */
  err: number[];
  count: number[];
}

/**
 * Inverse-variance weighted means of y in equal-width redshift bins over
 * [Z_MIN, Z_MAX], using objects with mask[i] set and finite y and err. The
 * error comes from the scatter inside the bin, not from the quoted errors, so it
 * stays honest if those are mis-sized. Bins with fewer than 3 objects are left out.
 */
export function binMeans(z: ArrayLike<number>, y: ArrayLike<number>, err: ArrayLike<number>, mask: ArrayLike<number>, nbins = 7): Bins {
  const out: Bins = { z: [], mean: [], err: [], count: [] };
  const width = (Z_MAX - Z_MIN) / nbins;
  const members: number[][] = Array.from({ length: nbins }, () => []);
  for (let i = 0; i < z.length; i++) {
    if (!mask[i] || !Number.isFinite(y[i]) || !(err[i] > 0)) continue;
    members[Math.min(nbins - 1, Math.max(0, Math.floor((z[i] - Z_MIN) / width)))].push(i);
  }
  for (const idx of members) {
    const n = idx.length;
    if (n < 3) continue;
    let sw = 0;
    let swy = 0;
    let swz = 0;
    for (const i of idx) {
      const w = 1 / err[i] ** 2;
      sw += w;
      swy += w * y[i];
      swz += w * z[i];
    }
    const mean = swy / sw;
    let sr = 0;
    for (const i of idx) sr += ((y[i] - mean) ** 2) / err[i] ** 2;
    out.z.push(swz / sw);
    out.mean.push(mean);
    out.err.push(Math.sqrt(sr / ((n - 1) * sw)));
    out.count.push(n);
  }
  return out;
}

export interface WFit {
  /** Best-fitting constant w, and the range where χ² is within 1 of its minimum. */
  best: number;
  lo: number;
  hi: number;
  chi2: number;
  dof: number;
  /** Which of the best value and the range ends sit on the edge of the scanned interval. */
  bestAtEdge: boolean;
  loAtEdge: boolean;
  hiAtEdge: boolean;
}

export const W_SCAN: [number, number] = [-3, 1];
const COARSE = 0.05;
const FINE = 0.005;

/**
 * Fits a constant w (wa = 0) to distance moduli mu with errors err (objects
 * with mask[i] set) by a χ² scan, with Ωm and H0 held at the generating values
 * and no free magnitude offset: the zero-points of the toy relations are taken
 * as known. (In a real analysis that absolute calibration is degenerate with
 * H0.) A coarse scan over −3 ≤ w ≤ 1 in steps of 0.05, then a fine one in
 * steps of 0.005 within 0.3 of its minimum.
 */
export function fitConstantW(z: ArrayLike<number>, mu: ArrayLike<number>, err: ArrayLike<number>, mask: ArrayLike<number>, cosmology: CPLParams): WFit | null {
  const idx: number[] = [];
  for (let i = 0; i < z.length; i++) if (mask[i] && Number.isFinite(mu[i]) && err[i] > 0) idx.push(i);
  if (idx.length < 3) return null;
  idx.sort((a, b) => z[a] - z[b]);
  const zs = Float64Array.from(idx, (i) => z[i]);
  const ys = Float64Array.from(idx, (i) => mu[i]);
  const inv = Float64Array.from(idx, (i) => 1 / err[i] ** 2);
  const offset = Float64Array.from(zs, (v) => 5 * Math.log10(1 + v));
  const chi = new Float64Array(zs.length);
  const dh = 5 * Math.log10(299792.458 / cosmology.h0) + 25;

  const chi2 = (w: number): number => {
    comovingSorted({ ...cosmology, w0: w, wa: 0 }, zs, chi);
    let sum = 0;
    for (let k = 0; k < zs.length; k++) {
      const d = ys[k] - (dh + offset[k] + 5 * Math.log10(chi[k]));
      sum += d * d * inv[k];
    }
    return sum;
  };
  const scan = (from: number, to: number, step: number) => {
    const ws: number[] = [];
    const cs: number[] = [];
    for (let w = from; w <= to + 1e-9; w += step) {
      ws.push(w);
      cs.push(chi2(w));
    }
    return { ws, cs };
  };

  const coarse = scan(W_SCAN[0], W_SCAN[1], COARSE);
  const centre = coarse.ws[coarse.cs.indexOf(Math.min(...coarse.cs))];
  const fine = scan(Math.max(W_SCAN[0], centre - 0.3), Math.min(W_SCAN[1], centre + 0.3), FINE);
  const min = Math.min(...fine.cs);
  const k = fine.cs.indexOf(min);
  // Parabola through the three lowest points for the minimum between grid nodes.
  let best = fine.ws[k];
  let chi2min = min;
  if (k > 0 && k < fine.ws.length - 1) {
    const [a, b, c] = [fine.cs[k - 1], fine.cs[k], fine.cs[k + 1]];
    const denom = a - 2 * b + c;
    if (denom > 0) {
      best += (FINE * (a - c)) / (2 * denom);
      chi2min = b - (a - c) ** 2 / (8 * denom);
    }
  }
  // Where χ² − χ²min crosses 1 on each side, by linear interpolation: on the fine
  // grid if the crossing is inside it, otherwise on the coarse one.
  const cross = (grid: { ws: number[]; cs: number[] }, from: number, dir: 1 | -1): number | null => {
    for (let j = from; j + dir >= 0 && j + dir < grid.ws.length; j += dir) {
      const c0 = grid.cs[j] - chi2min;
      const c1 = grid.cs[j + dir] - chi2min;
      if (c1 >= 1 && c0 < 1) return grid.ws[j] + ((1 - c0) / (c1 - c0)) * (grid.ws[j + dir] - grid.ws[j]);
    }
    return null;
  };
  const coarseAt = (w: number) => Math.max(0, Math.min(coarse.ws.length - 1, Math.round((w - W_SCAN[0]) / COARSE)));
  const side = (dir: 1 | -1): [number, boolean] => {
    const f = cross(fine, k, dir);
    if (f !== null) return [f, false];
    const c = cross(coarse, coarseAt(best), dir);
    return c !== null ? [c, false] : [W_SCAN[dir > 0 ? 1 : 0], true];
  };
  const [lo, loAtEdge] = side(-1);
  const [hi, hiAtEdge] = side(1);
  const bestAtEdge = best <= W_SCAN[0] + FINE || best >= W_SCAN[1] - FINE;
  return { best, lo, hi, chi2: chi2min, dof: idx.length - 1, bestAtEdge, loAtEdge, hiAtEdge };
}
