// Distances and expansion history of a flat universe containing matter and a
// dark energy whose equation of state evolves with the scale factor a:
//
//   w(a) = w0 + wa (1 − a),   i.e.   w(z) = w0 + wa z / (1 + z)
//
// (Chevallier & Polarski 2001; Linder 2003). w0 is the value today and wa sets
// how it changes; w0 = −1, wa = 0 is the cosmological constant. Radiation is
// ignored: with Ω_r ≈ 9 × 10⁻⁵ it would change distances by less than 0.05 per
// cent out to z = 10 (and the expansion rate itself by 0.17 per cent at z = 10).
//
// The Hubble rate in units of H0 is
//
//   E(z)² = Ωm (1 + z)³ + (1 − Ωm) f(z),
//   f(z)  = ρ_DE(z) / ρ_DE,0 = (1 + z)^{3(1 + w0 + wa)} exp(−3 wa z / (1 + z)),
//
// and everything else follows from it:
//
//   comoving distance        D_C(z) = (c/H0) χ(z),   χ(z) = ∫₀ᶻ dz′ / E(z′)
//   luminosity distance      D_L = (1 + z) D_C
//   angular diameter         D_A(z1, z2) = [D_C(z2) − D_C(z1)] / (1 + z2)   (flat)
//   distance modulus         μ = 5 log10(D_L / 10 pc)
//   lensing efficiency       D_ls / D_s = 1 − χ(z_l) / χ(z_s)               (flat)
//   time-delay distance      D_Δt = (1 + z_l) D_l D_s / D_ls                ∝ 1/H0
//   critical density         Σ_cr = c² D_s / (4π G D_l D_ls)                ∝ H0
//   deceleration parameter   q(z) = −1 + (1 + z) E′/E = ½ [1 + 3 w(z) Ω_DE(z)]
//   age                      t0 = (1/H0) ∫₀^∞ dz / [(1 + z) E(z)]
//
// χ(z) is integrated with five-point Gauss–Legendre on segments of width 0.1 in
// z: the cumulative integral is stored at the segment ends and the last, partial
// segment is integrated on demand, so a distance costs five evaluations of E and
// is accurate to about 1e-12 relative. No DOM; this module runs anywhere.

/** Speed of light, km/s. */
export const C_KMS = 299792.458;
/** 1/H0 in Gyr for H0 in km/s/Mpc: (1 Mpc / 1 km) / (1 Gyr in s). */
const HUBBLE_TIME_GYR = 977.7922216807891;
// The constants of src/lib/lens/cosmology.ts, for the lensing quantities below.
const G_SI = 6.6743e-11; // m³ kg⁻¹ s⁻²
const MPC_M = 3.0857e22;
const MSUN_KG = 1.989e30;
/** c² / (4π G) in M☉/kpc: Σ_cr in M☉/kpc² is SIGMA_CR_UNIT · D_s / (D_l D_ls) with the distances in kpc. */
const SIGMA_CR_UNIT = ((C_KMS * 1000) ** 2 / (4 * Math.PI * G_SI)) * (MPC_M / 1000) / MSUN_KG;
/** One M☉/kpc² in g/cm². */
export const MSUN_PER_KPC2_IN_G_PER_CM2 = (MSUN_KG * 1000) / (MPC_M * 0.1) ** 2;

export interface CPLParams {
  /** Hubble constant, km/s/Mpc. */
  h0: number;
  /** Matter density today in units of the critical density. Dark energy is 1 − Ωm (flat). */
  omegaM: number;
  /** Dark energy equation of state today. */
  w0: number;
  /** Its evolution: w(a) = w0 + wa (1 − a). */
  wa: number;
}

/** The ΛCDM model with the parameters of src/lib/lens/cosmology.ts. */
export const LCDM: CPLParams = { h0: 70, omegaM: 0.3, w0: -1, wa: 0 };

/** The same Ωm and H0 with a cosmological constant: the reference the pages compare against. */
export const withLambda = (p: CPLParams): CPLParams => ({ ...p, w0: -1, wa: 0 });

// Gauss–Legendre nodes and weights on [−1, 1].
const GL5_X = [-0.906179845938664, -0.5384693101056831, 0, 0.5384693101056831, 0.906179845938664];
const GL5_W = [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891];
const GL3_X = [-0.7745966692414834, 0, 0.7745966692414834];
const GL3_W = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];

/** w(z) = w0 + wa z / (1 + z). */
export const equationOfState = (p: CPLParams, z: number): number => p.w0 + (p.wa * z) / (1 + z);

/** Dark energy density relative to today, ρ_DE(z) / ρ_DE,0. */
export function darkEnergyRatio(p: CPLParams, z: number): number {
  return (1 + z) ** (3 * (1 + p.w0 + p.wa)) * Math.exp((-3 * p.wa * z) / (1 + z));
}

/** E(z)² = H(z)² / H0². */
export function hubble2(p: CPLParams, z: number): number {
  return p.omegaM * (1 + z) ** 3 + (1 - p.omegaM) * darkEnergyRatio(p, z);
}

/** Deceleration parameter q(z) = ½ [1 + 3 w(z) Ω_DE(z)] (matter has w = 0). */
export function deceleration(p: CPLParams, z: number): number {
  const omegaDE = ((1 - p.omegaM) * darkEnergyRatio(p, z)) / hubble2(p, z);
  return 0.5 * (1 + 3 * equationOfState(p, z) * omegaDE);
}

/** ∫ dz / E over [z0, z1] by n-point Gauss–Legendre (n = 3 or 5). */
function segment(p: CPLParams, z0: number, z1: number, five = true): number {
  const mid = 0.5 * (z0 + z1);
  const half = 0.5 * (z1 - z0);
  const xs = five ? GL5_X : GL3_X;
  const ws = five ? GL5_W : GL3_W;
  let sum = 0;
  for (let k = 0; k < xs.length; k++) sum += ws[k] / Math.sqrt(hubble2(p, mid + half * xs[k]));
  return sum * half;
}

/**
 * χ(z) = ∫₀ᶻ dz′/E in units of c/H0 for redshifts sorted in ascending order,
 * written into `out`. Integrates from one redshift to the next, so N distances
 * cost about 3N evaluations of E; used for fitting many objects at once.
 */
export function comovingSorted(p: CPLParams, zs: ArrayLike<number>, out: Float64Array, maxStep = 0.05): Float64Array {
  let z = 0;
  let chi = 0;
  for (let i = 0; i < zs.length; i++) {
    const target = zs[i];
    while (target - z > maxStep) {
      chi += segment(p, z, z + maxStep, false);
      z += maxStep;
    }
    if (target > z) chi += segment(p, z, target, false);
    z = Math.max(z, target);
    out[i] = chi;
  }
  return out;
}

/** A sign change of q(z): where the expansion starts (or stops) accelerating, going forward in time. */
export interface Transition {
  z: number;
  /** True if the expansion decelerates before (higher z) and accelerates after (lower z). */
  begins: boolean;
}

/**
 * One cosmology with its distances precomputed. Build one per parameter set
 * (it costs about 500 evaluations of E) and query it as often as needed.
 */
export class Universe {
  readonly params: CPLParams;
  /** False if E² is not positive and finite somewhere on 0 ≤ z ≤ zmax; every method then returns NaN. */
  readonly ok: boolean;
  /** Hubble distance c/H0, Mpc. */
  readonly hubbleDistance: number;

  private readonly step: number;
  private readonly zmax: number;
  /** χ at z = i · step. */
  private readonly nodes: Float64Array;

  constructor(params: CPLParams, zmax = 10, step = 0.1) {
    this.params = { ...params };
    this.hubbleDistance = C_KMS / params.h0;
    this.step = step;
    this.zmax = zmax;
    const n = Math.ceil(zmax / step);
    this.nodes = new Float64Array(n + 1);

    // For 0 < Ωm < 1 both terms of E² are positive, so this guard only trips on
    // overflow or on parameters outside the physical range (Ωm ≥ 1 with a
    // growing dark energy term, NaN inputs). Charts check `ok` and say so.
    let ok = Number.isFinite(params.h0) && params.h0 > 0;
    for (let i = 0; i <= 4 * n && ok; i++) {
      const e2 = hubble2(params, (i * step) / 4);
      if (!(e2 > 0) || !Number.isFinite(e2)) ok = false;
    }
    this.ok = ok;
    if (!ok) {
      this.nodes.fill(NaN);
      return;
    }
    for (let i = 1; i <= n; i++) this.nodes[i] = this.nodes[i - 1] + segment(params, (i - 1) * step, i * step);
  }

  /** E(z) = H(z) / H0. */
  E(z: number): number {
    return Math.sqrt(hubble2(this.params, z));
  }
  w(z: number): number {
    return equationOfState(this.params, z);
  }
  darkEnergyRatio(z: number): number {
    return darkEnergyRatio(this.params, z);
  }
  q(z: number): number {
    return this.ok ? deceleration(this.params, z) : NaN;
  }

  /** χ(z) = ∫₀ᶻ dz′/E, dimensionless (units of c/H0). */
  chi(z: number): number {
    if (!this.ok) return NaN;
    if (z <= 0) return 0;
    if (z > this.zmax) {
      // Beyond the table: continue from its end in steps of the table's width.
      let chi = this.nodes[this.nodes.length - 1];
      let a = (this.nodes.length - 1) * this.step;
      while (z - a > this.step) {
        chi += segment(this.params, a, a + this.step);
        a += this.step;
      }
      return chi + segment(this.params, a, z);
    }
    const i = Math.floor(z / this.step);
    const a = i * this.step;
    return this.nodes[i] + (z > a ? segment(this.params, a, z) : 0);
  }

  /** Comoving (line-of-sight) distance, Mpc. */
  comoving(z: number): number {
    return this.hubbleDistance * this.chi(z);
  }
  /** Luminosity distance, Mpc. */
  luminosity(z: number): number {
    return (1 + z) * this.comoving(z);
  }
  /** Angular-diameter distance to z, Mpc. */
  angular(z: number): number {
    return this.comoving(z) / (1 + z);
  }
  /** Angular-diameter distance from z1 to z2 > z1 (flat universe), Mpc. */
  angularBetween(z1: number, z2: number): number {
    return (this.comoving(z2) - this.comoving(z1)) / (1 + z2);
  }
  /** Distance modulus μ = 5 log10(D_L / 10 pc) = 5 log10(D_L / Mpc) + 25, mag. */
  distanceModulus(z: number): number {
    return 5 * Math.log10(this.luminosity(z)) + 25;
  }
  /** Lensing efficiency D_ls / D_s for a lens at zl and a source at zs; zero for a source in front of the lens. */
  lensingRatio(zl: number, zs: number): number {
    if (!this.ok) return NaN;
    if (zs <= zl) return 0;
    return 1 - this.chi(zl) / this.chi(zs);
  }
  /**
   * Time-delay distance D_Δt = (1 + zl) D_l D_s / D_ls, Mpc (Refsdal 1964): a
   * time delay between two images is D_Δt / c times a difference of Fermat
   * potentials set by the lens model, so for a fixed lens model delays scale with
   * D_Δt, and D_Δt ∝ 1/H0. NaN for a source in front of the lens.
   */
  timeDelayDistance(zl: number, zs: number): number {
    if (!this.ok || !(zs > zl)) return NaN;
    return ((1 + zl) * this.angular(zl) * this.angular(zs)) / this.angularBetween(zl, zs);
  }
  /**
   * Critical surface density Σ_cr = c² D_s / (4π G D_l D_ls), M☉/kpc² (multiply
   * by MSUN_PER_KPC2_IN_G_PER_CM2 for g/cm²). Σ_cr ∝ H0. The mass inside an
   * Einstein radius θ_E (radians) is π (θ_E D_l)² Σ_cr = c² θ_E² D_l D_s / (4 G D_ls).
   * NaN for a source in front of the lens.
   */
  criticalSurfaceDensity(zl: number, zs: number): number {
    if (!this.ok || !(zs > zl)) return NaN;
    const kpc = 1000; // kpc per Mpc
    return (SIGMA_CR_UNIT * this.angular(zs)) / (this.angular(zl) * this.angularBetween(zl, zs) * kpc);
  }

  /**
   * Age of the Universe, Gyr: t0 = (1/H0) ∫₀¹ da / (a E). With a = u² this is
   * ∫₀¹ 2 du / (u E), whose integrand vanishes smoothly at u = 0 in the matter
   * era (E ∝ u⁻³), so composite Gauss–Legendre converges quickly.
   */
  age(): number {
    if (!this.ok) return NaN;
    const segments = 64;
    let sum = 0;
    for (let s = 0; s < segments; s++) {
      const mid = (s + 0.5) / segments;
      const half = 0.5 / segments;
      for (let k = 0; k < 5; k++) {
        const u = mid + half * GL5_X[k];
        const z = 1 / (u * u) - 1;
        sum += (GL5_W[k] * half * 2) / (u * this.E(z));
      }
    }
    return (HUBBLE_TIME_GYR / this.params.h0) * sum;
  }

  /**
   * Every redshift in [0, zmax] where q(z) changes sign, lowest first. For
   * ΛCDM there is one, where acceleration begins; an evolving w can give none
   * (always decelerating, or accelerating since before zmax) or several.
   */
  accelerationTransitions(zmax = 20, step = 0.01): Transition[] {
    if (!this.ok) return [];
    const out: Transition[] = [];
    let z0 = 0;
    let q0 = this.q(0);
    for (let i = 1; i * step <= zmax + 1e-9; i++) {
      const z1 = i * step;
      const q1 = this.q(z1);
      if (q0 === 0 && i === 1) out.push({ z: 0, begins: q1 > 0 });
      else if (q0 * q1 < 0) {
        // Bisection to 1e-10 in z.
        let lo = z0;
        let hi = z1;
        let qlo = q0;
        while (hi - lo > 1e-10) {
          const mid = 0.5 * (lo + hi);
          const qm = this.q(mid);
          if (qm * qlo <= 0) hi = mid;
          else {
            lo = mid;
            qlo = qm;
          }
        }
        // q < 0 just below the root (later times) and > 0 above it: acceleration begins here.
        out.push({ z: 0.5 * (lo + hi), begins: q0 < 0 });
      }
      z0 = z1;
      q0 = q1;
    }
    return out;
  }
}
