// Distances in a flat ΛCDM universe, and the lensing quantities built on them.

export const H0 = 70; // km/s/Mpc
export const OMEGA_M = 0.3;
export const C_KMS = 299792.458;

const ARCSEC_PER_RAD = 206264.806;
const MPC_M = 3.0857e22;
const MSUN_KG = 1.989e30;
const G_SI = 6.6743e-11;

const invE = (z: number) => 1 / Math.sqrt(OMEGA_M * (1 + z) ** 3 + 1 - OMEGA_M);

/** Comoving distance to redshift z, in units of c/H0 (Simpson's rule). */
export function comoving(z: number): number {
  const n = 256;
  const h = z / n;
  let sum = invE(0) + invE(z);
  for (let i = 1; i < n; i++) sum += (i % 2 ? 4 : 2) * invE(i * h);
  return (sum * h) / 3;
}

/** Lensing efficiency D_ls / D_s for a lens at zl and a source at zs. Zero for a source in front of the lens. */
export function efficiency(zl: number, zs: number): number {
  if (zs <= zl) return 0;
  const dl = comoving(zl);
  const ds = comoving(zs);
  return (ds - dl) / ds;
}

/** Einstein radius of a singular isothermal sphere with velocity dispersion sigma (km/s), for D_ls/D_s = 1, in arcsec. */
export function sisEinsteinRadius(sigma: number): number {
  return 4 * Math.PI * (sigma / C_KMS) ** 2 * ARCSEC_PER_RAD;
}

/** Projected mass inside an Einstein radius theta (arcsec), in solar masses: M = c² θ² D_l D_s / (4 G D_ls). */
export function einsteinMass(theta: number, zl: number, zs: number): number {
  const hubble = (C_KMS / H0) * MPC_M; // c/H0 in metres
  const chiL = comoving(zl);
  const chiS = comoving(zs);
  if (chiS <= chiL) return 0;
  const dl = (hubble * chiL) / (1 + zl);
  // D_s / D_ls for a flat universe: the (1 + zs) factors cancel.
  const ratio = chiS / (chiS - chiL);
  const rad = theta / ARCSEC_PER_RAD;
  const c = C_KMS * 1000;
  return (c * c * rad * rad * dl * ratio) / (4 * G_SI) / MSUN_KG;
}

/** Proper kiloparsecs per arcsecond at redshift z. */
export function kpcPerArcsec(z: number): number {
  return (((C_KMS / H0) * comoving(z)) / (1 + z) / ARCSEC_PER_RAD) * 1000;
}
