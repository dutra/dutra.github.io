// Galaxy-galaxy strong lensing (GGSL) by one cluster member, for the page at
// /projects/ggsl/. Pure functions and typed arrays, no DOM.
//
// The model, in the few arcseconds around the member:
//
//  - The member is a dual pseudo-isothermal sphere (dPIE; Elíasdóttir et al.
//    2007) with velocity dispersion sigma, core radius rc and truncation radius
//    rt. Its deflection is radial, with magnitude
//        alpha(r) = A [sqrt(rc² + r²) - rc - sqrt(rt² + r²) + rt] / r,
//        A = b rt / (rt - rc),  b = 4 pi (sigma/c)² D_ls/D_s,
//    which tends to the singular isothermal sphere (alpha = b) for rc -> 0 and
//    rt -> infinity. This is the profile src/lib/lens/model.ts uses for members.
//    Those limits are allowed exactly: rc = 0 and rt = Infinity give the cored
//    isothermal sphere, alpha = b [sqrt(rc² + r²) - rc] / r (rt = Infinity), and the
//    singular isothermal sphere, alpha = b (both). makeLens takes A = b when rt is
//    infinite instead of forming rt / (rt - rc), and the rationalised forms below
//    then hold as written (1/Infinity = 0, and sqrt(r²) + 0 = r).
//  - The cluster's smooth mass is a constant convergence kc and shear gc, with
//    the shear along x. Both are the values for the chosen source redshift.
//  - Lens equation, member at the origin:
//        beta_x = (1 - kc - gc) theta_x - alpha_x,
//        beta_y = (1 - kc + gc) theta_y - alpha_y.
//
// Angles are in arcsec, x to the right and y up. Everything that depends on the
// source redshift enters through b (the member) and through kc, gc (the cluster).

/** Gravitational constant in kpc (km/s)² per solar mass. */
const G_KPC = 4.30091727e-6;

export interface LensParams {
  /** Lensing strength b = sisEinsteinRadius(sigma) × D_ls/D_s, arcsec. */
  b: number;
  /** dPIE core radius, arcsec; 0 for a singular profile. */
  core: number;
  /** dPIE truncation radius, arcsec; larger than the core, or Infinity for no truncation. */
  cut: number;
  /** Cluster convergence and shear at the member, for this source redshift. */
  kappa: number;
  gamma: number;
}

/** A lens with its constants worked out. */
export interface Lens extends LensParams {
  /** Amplitude A = b rt / (rt - rc), or b when rt is infinite. */
  amp: number;
  /** Diagonal of the cluster's part of the Jacobian. */
  m11: number;
  m22: number;
}

export function makeLens(p: LensParams): Lens {
  // A core of exactly 0 is the singular limit; a tiny positive one is kept off rounding noise.
  const core = p.core > 0 ? Math.max(p.core, 1e-6) : 0;
  const cut = Number.isFinite(p.cut) ? Math.max(p.cut, core * 1.001, 1e-6) : Infinity;
  return {
    ...p,
    core,
    cut,
    amp: Number.isFinite(cut) ? (p.b * cut) / (cut - core) : p.b,
    m11: 1 - p.kappa - p.gamma,
    m22: 1 - p.kappa + p.gamma,
  };
}

/** True for a profile with no core (rc = 0), whose centre is a singularity. */
export function isSingular(lens: Lens): boolean {
  return lens.core === 0;
}

/**
 * Mean convergence of the member inside radius r, kbar = alpha(r) / r.
 * Written as A [1/(sqrt(rc² + r²) + rc) - 1/(sqrt(rt² + r²) + rt)], which is the
 * same expression rationalised, so it stays accurate at r = 0 (where it equals
 * the central convergence; infinite for a singular profile).
 */
export function meanConvergence(lens: Lens, r: number): number {
  const r2 = r * r;
  return lens.amp * (1 / (Math.sqrt(lens.core * lens.core + r2) + lens.core) - 1 / (Math.sqrt(lens.cut * lens.cut + r2) + lens.cut));
}

/** Convergence of the member at radius r: kappa = (A/2) [1/sqrt(rc² + r²) - 1/sqrt(rt² + r²)]. */
export function convergence(lens: Lens, r: number): number {
  const r2 = r * r;
  return 0.5 * lens.amp * (1 / Math.sqrt(lens.core * lens.core + r2) - 1 / Math.sqrt(lens.cut * lens.cut + r2));
}

/**
 * The member's largest deflection, with a small safety factor. alpha(r) rises
 * from 0, peaks near the truncation radius and falls as b rt / r, so a scan in
 * log r finds the peak. (The analytic bounds alpha < A and alpha < b rt / r are
 * loose when rt is close to rc, which would waste the mesh on empty sky.)
 */
function maxDeflection(lens: Lens): number {
  // Untruncated (cored or singular isothermal): alpha rises monotonically to b.
  if (!Number.isFinite(lens.cut)) return 1.05 * lens.b;
  let peak = 0;
  const r0 = 0.01 * (lens.core > 0 ? lens.core : 1e-3 * lens.cut);
  const r1 = 100 * lens.cut;
  for (let k = 0; k <= 400; k++) {
    const r = r0 * (r1 / r0) ** (k / 400);
    peak = Math.max(peak, r * meanConvergence(lens, r));
  }
  return 1.05 * peak;
}

/** Source-plane position of the ray through (x, y), written into out. */
export function toSource(lens: Lens, x: number, y: number, out: Float64Array | number[]): void {
  const kbar = meanConvergence(lens, Math.hypot(x, y));
  out[0] = (lens.m11 - kbar) * x;
  out[1] = (lens.m22 - kbar) * y;
}

/**
 * Determinant of the lensing Jacobian at (x, y), analytically. For a radial
 * deflection alpha = kbar(r) theta, d alpha_i / d theta_j = kbar delta_ij +
 * 2 (kappa - kbar) theta_i theta_j / r², so
 *   A = diag(m11, m22) - kbar I - 2 (kappa - kbar) n nᵀ,  n = theta / r.
 */
export function jacobianDet(lens: Lens, x: number, y: number): number {
  const r2 = x * x + y * y;
  const r = Math.sqrt(r2);
  const kbar = meanConvergence(lens, r);
  const d = r2 > 0 ? (2 * (convergence(lens, r) - kbar)) / r2 : 0;
  const a11 = lens.m11 - kbar - d * x * x;
  const a22 = lens.m22 - kbar - d * y * y;
  const a12 = -d * x * y;
  return a11 * a22 - a12 * a12;
}

/**
 * Einstein radius, arcsec: the radius where the circularly averaged mean
 * convergence inside it, kc + kbar(r), equals 1. (The shear averages to zero
 * around a circle.) For gc = 0 this is the solution of alpha(r) = (1 - kc) r.
 * Returns 0 when the member is not dense enough at its centre to reach it.
 */
export function einsteinRadius(lens: Lens): number {
  const target = 1 - lens.kappa;
  if (meanConvergence(lens, 0) <= target) return 0;
  // kbar falls monotonically from its central value to 0, so bracket and bisect.
  let lo = 0;
  let hi = Math.max(lens.b, lens.core);
  while (meanConvergence(lens, hi) > target) hi *= 2;
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi);
    if (meanConvergence(lens, mid) > target) lo = mid;
    else hi = mid;
  }
  return 0.5 * (lo + hi);
}

/**
 * Total mass of the dPIE, in solar masses, for sigma in km/s and the
 * truncation radius in kpc: M = pi sigma² rt / G.
 *
 * Derived from the profile rather than quoted: the projected mass inside R is
 * M(<R) = pi R² Sigma_cr kbar(R) = pi Sigma_cr R alpha(R) (physical units), and
 * R alpha(R) -> A (rt - rc) = b rt as R -> infinity, so M = pi Sigma_cr b rt.
 * With b = 4 pi sigma² D_ls / (c² D_s) and Sigma_cr = c² D_s / (4 pi G D_l D_ls)
 * the distances cancel, leaving M = pi sigma² rt / G, independent of the core.
 *
 * Here sigma is the equivalent singular-isothermal dispersion, defined by
 * b = 4 pi (sigma/c)² D_ls/D_s. Elíasdóttir et al. (2007, arXiv:0710.5636), and
 * Lenstool, use a fiducial sigma_dPIE with E0 = 6 pi (sigma_dPIE/c)² D_ls/D_s and
 * sigma_dPIE² = (4/3) G Sigma0 rc rt² / (rt² - rc²), M = 2 pi Sigma0 rc rt. The
 * two are related by sigma² = (3/2) sigma_dPIE² (1 - rc²/rt²), so in their form
 * M = 3 pi sigma_dPIE² rt (1 - rc²/rt²) / (2 G).
 */
export function dpieMass(sigma: number, cutKpc: number): number {
  return (Math.PI * sigma * sigma * cutKpc) / G_KPC;
}

// --- eigenvalues of the Jacobian ------------------------------------------------

/**
 * The two eigenvalues of the lensing Jacobian at (x, y), written into out:
 *   out[0] = lambda_t = 1 - kappa_tot - |gamma_tot|  (tangential),
 *   out[1] = lambda_r = 1 - kappa_tot + |gamma_tot|  (radial),
 * where kappa_tot and gamma_tot include the cluster's kc, gc and the member's
 * own convergence and shear. With A as in jacobianDet, 1 - kappa_tot is half
 * the trace and |gamma_tot| = sqrt(((a11 - a22)/2)² + a12²); det A = lambda_t lambda_r.
 * The tangential critical curve is lambda_t = 0 and the radial one lambda_r = 0.
 */
export function eigenvalues(lens: Lens, x: number, y: number, out: Float64Array | number[]): void {
  const r2 = x * x + y * y;
  const r = Math.sqrt(r2);
  const kbar = meanConvergence(lens, r);
  const d = r2 > 0 ? (2 * (convergence(lens, r) - kbar)) / r2 : 0;
  const a11 = lens.m11 - kbar - d * x * x;
  const a22 = lens.m22 - kbar - d * y * y;
  const a12 = -d * x * y;
  const half = 0.5 * (a11 + a22);
  const g = Math.hypot(0.5 * (a11 - a22), a12);
  out[0] = half - g;
  out[1] = half + g;
}

/**
 * Distance from the member, along x (axis 0) or y (axis 1), of the outermost
 * image of a source at distance beta along that axis: the largest root of
 * (m - kbar(r)) r = beta, with m = 1 - kc - gc for x and 1 - kc + gc for y.
 *
 * It bounds every image of every source in the box |beta_x|, |beta_y| <= beta:
 * from the lens equation |theta_x| (m11 - kbar(r)) = |beta_x| with r >= |theta_x|,
 * and kbar falls with r, so either kbar(r) >= m11 (the image lies inside the
 * root of kbar = m11, which is inside this radius) or
 * (m11 - kbar(|theta_x|)) |theta_x| <= beta, which puts |theta_x| below the largest
 * root. The same holds along y.
 */
export function outerImageRadius(lens: Lens, beta: number, axis: 0 | 1): number {
  const m = axis === 0 ? lens.m11 : lens.m22;
  const f = (r: number) => (m - meanConvergence(lens, r)) * r - beta;
  let hi = Math.max(lens.b, lens.core, beta, 1e-3);
  while (f(hi) < 0) hi *= 2;
  // Scan down from hi for the last sign change, then bisect inside it.
  let lo = hi;
  for (let k = 1; k <= 400; k++) {
    const r = hi * 1e-4 ** (k / 400);
    if (f(r) < 0) {
      lo = r;
      break;
    }
    hi = r;
  }
  if (lo === hi) return hi;
  for (let i = 0; i < 60; i++) {
    const mid = 0.5 * (lo + hi);
    if (f(mid) < 0) lo = mid;
    else hi = mid;
  }
  return hi;
}

// --- critical curves, caustics and the cross-section --------------------------

export interface Resolution {
  /** Rings and spokes of the polar image-plane mesh. */
  rings: number;
  spokes: number;
  /** Cells of the source-plane multiplicity grid along its longer side. */
  grid: number;
  /** Scan lines used to measure the area inside the tangential caustic. */
  areaRows: number;
}

/** Good enough for the live figures. */
export const FINE: Resolution = { rings: 220, spokes: 288, grid: 288, areaRows: 2048 };
/** For the chart, which needs many models quickly. */
export const COARSE: Resolution = { rings: 110, spokes: 144, grid: 144, areaRows: 768 };

/** A polyline [x0, y0, x1, y1, ...] in arcsec; closed if its last point joins its first. */
export interface Curve {
  points: Float64Array;
  closed: boolean;
}

export interface Multiplicity {
  /** Number of images of a point source at each cell centre, row by row from the bottom (y up). */
  counts: Uint8Array;
  nx: number;
  ny: number;
  /** Lower-left corner of the grid and cell size, arcsec. */
  x0: number;
  y0: number;
  cell: number;
}

export interface Analysis {
  /** Tangential (lambda_t = 0) and radial (lambda_r = 0) critical curves, and their caustics. */
  tangentialCritical: Curve[];
  radialCritical: Curve[];
  tangentialCaustic: Curve[];
  radialCaustic: Curve[];
  /**
   * For a singular profile (no core), which has no radial critical curve: the
   * cut, the source-plane image of the galaxy's centre, a circle of radius
   * radius lim_{r->0} alpha(r) = A around the member. Crossing it adds or removes one image
   * (the one absorbed by the central singularity). Empty for a cored profile.
   */
  cut: Curve[];
  /** Largest |x| and |y| reached by all the critical curves and by all the caustics and the cut, arcsec (0 if there are none). */
  criticalExtent: [number, number];
  causticExtent: [number, number];
  /** Number of images over a box around the caustics, or null if there are no caustics. */
  multiplicity: Multiplicity | null;
  /**
   * GGSL cross-section, arcsec²: the source-plane area enclosed by the member's
   * tangential caustic, following Meneghetti et al. (2020).
   */
  crossSection: number;
  /** Source-plane area in which a point source has 3 or more images, arcsec². */
  multiplyImagedArea: number;
}

/**
 * Traces the critical curves and caustics and measures both areas.
 *
 * Curves. The image plane is covered by a polar mesh centred on the member:
 * rings spaced logarithmically from well inside the core out to a radius r_max,
 * and evenly spaced spokes. The two eigenvalues of the Jacobian are evaluated
 * at every vertex, and each zero contour is found by marching squares: on every
 * mesh edge where the eigenvalue changes sign, the zero is located by root
 * finding on the analytic eigenvalue (so the curve points are exact, not linear
 * interpolations), and the crossings are joined into polylines. Each polyline is
 * oriented with lambda < 0 on its left and mapped point by point to the source
 * plane, which gives the caustics. The log-spaced rings resolve the radial
 * critical curve even for a core of a few hundredths of an arcsecond.
 *
 * GGSL cross-section (area enclosed by the tangential caustic). The caustic is
 * a closed curve with cusps; it need not be simple (it can cross itself) and
 * there may be more than one. The area is that of the points with a nonzero
 * winding number about the caustic polylines, which is right for any closed
 * curve and does not depend on how each fold is oriented. It is integrated by
 * scan lines: on each of areaRows horizontal lines across the caustic, every
 * caustic segment that crosses the line adds +1 or -1 to the winding number
 * (by its direction), the crossings are sorted, and the lengths of the
 * intervals with nonzero winding are summed exactly. The area is the sum of
 * those lengths times the line spacing (a midpoint rule in y only).
 *
 * With no external shear (gc = 0) the model is circularly symmetric: every
 * point of the tangential critical curve, the Einstein ring, maps to beta = 0,
 * so the tangential caustic is a single point and the area is exactly 0. That
 * case returns 0 directly rather than measuring a polygon of rounding errors.
 *
 * Area with 3 or more images. Each mesh quad is split into two triangles and
 * every vertex is mapped to the source plane. This makes the lens mapping
 * piecewise linear and continuous, so the number of mapped triangles covering
 * a source-plane point is the number of images of that point (for the
 * piecewise-linear map, exactly: there are no gaps or overlaps between
 * neighbouring triangles, only folds along the critical curves). Counting
 * triangles over the cell centres of a grid around the caustics gives the image
 * multiplicity N(beta), and the area of the cells with N >= 3 is the second
 * number. It includes the sources inside the radial caustic as well as those
 * inside the tangential one.
 *
 * For a singular profile the centre maps not to beta = 0 but to the cut, the
 * circle |beta| = A, and the image a cored profile would have there is absorbed
 * by the singularity. So the central fan of triangles is left out and the mesh's
 * inner ring (at 0.1% of b) traces the cut; across it the count changes by one.
 * The area with 3 or more images is then the area inside the tangential caustic.
 *
 * r_max follows from the lens equation: |beta| >= (1 - kc - gc) |theta| - alpha_max,
 * so a ray with |theta| > (beta_max + alpha_max) / (1 - kc - gc) cannot reach a
 * source-plane point closer than beta_max to the member.
 */
export function analyse(lens: Lens, res: Resolution = FINE): Analysis {
  return run(lens, res, true);
}

/** The GGSL cross-section alone (area inside the tangential caustic), arcsec². Skips the multiplicity grid. */
export function crossSection(lens: Lens, res: Resolution = FINE): number {
  return run(lens, res, false).crossSection;
}

function run(lens: Lens, res: Resolution, withMultiplicity: boolean): Analysis {
  const { rings, spokes } = res;
  const out = [0, 0];
  const eig = [0, 0];
  const lamT = (x: number, y: number) => {
    eigenvalues(lens, x, y, eig);
    return eig[0];
  };
  const lamR = (x: number, y: number) => {
    eigenvalues(lens, x, y, eig);
    return eig[1];
  };

  // Inner and outer radii of the mesh. The radial critical curve of a cored
  // profile lies at a radius of order the core, so start well inside it.
  const singular = isSingular(lens);
  const rMin = singular ? 1e-3 * Math.max(lens.b, 1e-3) : 0.02 * Math.min(lens.core, Math.max(lens.b, 1e-3));
  const lambda = Math.max(lens.m11, 1e-3);
  // A first outer radius for finding the curves: the tangential critical curve
  // lies where the member's mean convergence is about lambda, inside b / lambda.
  let rMax = (2.5 * maxDeflection(lens)) / lambda + 4 * lens.core;

  const nv = rings * spokes;
  const mesh: Mesh = {
    rings,
    spokes,
    vx: new Float64Array(nv),
    vy: new Float64Array(nv),
  };
  const bx = new Float64Array(nv);
  const by = new Float64Array(nv);
  const lt = new Float64Array(nv);
  const lr = new Float64Array(nv);
  const cosS = new Float64Array(spokes);
  const sinS = new Float64Array(spokes);
  // Spokes are offset by half a step so no vertex sits on an axis, where the
  // symmetric model would put mapped vertices exactly on grid-cell centres.
  for (let j = 0; j < spokes; j++) {
    const phi = ((j + 0.5) * 2 * Math.PI) / spokes;
    cosS[j] = Math.cos(phi);
    sinS[j] = Math.sin(phi);
  }

  const fillMesh = () => {
    const step = Math.log(rMax / rMin) / (rings - 1);
    for (let i = 0; i < rings; i++) {
      const r = rMin * Math.exp(i * step);
      for (let j = 0; j < spokes; j++) {
        const v = i * spokes + j;
        const x = r * cosS[j];
        const y = r * sinS[j];
        mesh.vx[v] = x;
        mesh.vy[v] = y;
        toSource(lens, x, y, out);
        bx[v] = out[0];
        by[v] = out[1];
        eigenvalues(lens, x, y, eig);
        lt[v] = eig[0];
        lr[v] = eig[1];
      }
    }
  };

  const traceBoth = () => {
    const t = traceContour(mesh, lt, lamT);
    const r = traceContour(mesh, lr, lamR);
    return { t, r, extent: Math.max(t.extent, r.extent) };
  };

  fillMesh();
  let traced = traceBoth();
  // The tangential curve can reach further than the first guess when kc + gc is
  // close to 1; widen the mesh until it sits comfortably inside.
  for (let k = 0; k < 4 && traced.extent > 0.6 * rMax; k++) {
    rMax *= 2;
    fillMesh();
    traced = traceBoth();
  }

  const toCaustic = (c: Curve): Curve => {
    const p = c.points;
    const q = new Float64Array(p.length);
    for (let s = 0; s < p.length; s += 2) {
      toSource(lens, p[s], p[s + 1], out);
      q[s] = out[0];
      q[s + 1] = out[1];
    }
    return { points: q, closed: c.closed };
  };
  const tangentialCritical = traced.t.curves;
  const radialCritical = traced.r.curves;
  const tangentialCaustic = tangentialCritical.map(toCaustic);
  const radialCaustic = radialCritical.map(toCaustic);
  const cut: Curve[] = [];
  if (singular) {
    const n = 256;
    const points = new Float64Array(2 * n);
    for (let k = 0; k < n; k++) {
      const phi = (2 * Math.PI * k) / n;
      points[2 * k] = lens.amp * Math.cos(phi);
      points[2 * k + 1] = lens.amp * Math.sin(phi);
    }
    cut.push({ points, closed: true });
  }
  const criticalExtent = extentOf([...tangentialCritical, ...radialCritical]);
  const causticExtent = extentOf([...tangentialCaustic, ...radialCaustic, ...cut]);

  // The cross-section: area inside the tangential caustic (see above for gc = 0).
  let cross = 0;
  if (lens.gamma !== 0 && tangentialCaustic.length > 0) {
    const ce = extentOf(tangentialCaustic);
    const te = extentOf(tangentialCritical);
    // A caustic many orders of magnitude smaller than its critical curve is a point.
    if (Math.max(ce[0], ce[1]) > 1e-9 * Math.max(te[0], te[1])) cross = windingArea(tangentialCaustic, res.areaRows);
  }

  const base = {
    tangentialCritical,
    radialCritical,
    tangentialCaustic,
    radialCaustic,
    cut,
    criticalExtent,
    causticExtent,
    crossSection: cross,
  };
  if (!withMultiplicity || Math.max(causticExtent[0], causticExtent[1]) === 0) {
    return { ...base, multiplicity: null, multiplyImagedArea: 0 };
  }

  // Source-plane box around the caustics (multiplicity only changes across
  // them), with a margin, and square cells.
  const hx = causticExtent[0] * 1.06 + 1e-6;
  const hy = causticExtent[1] * 1.06 + 1e-6;
  const cell = (2 * Math.max(hx, hy)) / res.grid;
  // Even cell counts keep cell centres off the axes.
  const nx = 2 * Math.ceil(hx / cell);
  const ny = 2 * Math.ceil(hy / cell);
  const x0 = (-nx / 2) * cell;
  const y0 = (-ny / 2) * cell;

  // Make sure every ray that can land in the box is inside the mesh.
  const betaMax = Math.hypot(nx, ny) * 0.5 * cell;
  const rNeed = (1.1 * (betaMax + maxDeflection(lens))) / lambda;
  if (rNeed > rMax) {
    rMax = rNeed;
    fillMesh();
  }

  const counts = new Uint8Array(nx * ny);
  const raster = (ax: number, ay: number, bxv: number, byv: number, cx: number, cy: number) =>
    rasterTriangle(counts, nx, ny, x0, y0, cell, ax, ay, bxv, byv, cx, cy);

  // Central fan: the origin maps to beta = 0. (Not for a singular profile; see above.)
  for (let j = 0; j < spokes && !singular; j++) {
    const j2 = (j + 1) % spokes;
    raster(0, 0, bx[j], by[j], bx[j2], by[j2]);
  }
  // Ring quads, two triangles each.
  for (let i = 0; i < rings - 1; i++) {
    for (let j = 0; j < spokes; j++) {
      const j2 = (j + 1) % spokes;
      const a = i * spokes + j;
      const b = i * spokes + j2;
      const c = (i + 1) * spokes + j2;
      const d = (i + 1) * spokes + j;
      raster(bx[a], by[a], bx[b], by[b], bx[c], by[c]);
      raster(bx[a], by[a], bx[c], by[c], bx[d], by[d]);
    }
  }

  // A round singular profile with no shear has at most two images (the Einstein
  // ring maps to the single point beta = 0); the folded mesh triangles would
  // otherwise count a cell or two there.
  if (singular && lens.gamma === 0) for (let p = 0; p < counts.length; p++) counts[p] = Math.min(counts[p], 2);
  let multiple = 0;
  for (let p = 0; p < counts.length; p++) if (counts[p] >= 3) multiple++;

  return {
    ...base,
    multiplicity: { counts, nx, ny, x0, y0, cell },
    multiplyImagedArea: multiple * cell * cell,
  };
}

/** Largest |x| and |y| over a set of curves. */
function extentOf(curves: Curve[]): [number, number] {
  let ex = 0;
  let ey = 0;
  for (const c of curves) {
    const p = c.points;
    for (let s = 0; s < p.length; s += 2) {
      ex = Math.max(ex, Math.abs(p[s]));
      ey = Math.max(ey, Math.abs(p[s + 1]));
    }
  }
  return [ex, ey];
}

interface Mesh {
  rings: number;
  spokes: number;
  /** Vertex positions, ring by ring, arcsec. */
  vx: Float64Array;
  vy: Float64Array;
}

/**
 * Zero contour of a field sampled at the mesh vertices (values in val, the
 * analytic function in field), as oriented polylines.
 *
 * Marching squares: each mesh edge whose end values differ in sign holds one
 * crossing, located by the Illinois variant of regula falsi on the analytic
 * field along the edge, and stored once per edge so that neighbouring quads
 * share it exactly. A quad with two crossings gives one segment; one with four
 * (a saddle, rare at this resolution) is joined in a fixed order. Segments that
 * share an edge are then chained into polylines, closed where the chain comes
 * back to its start. Each polyline is turned so the field is negative on its
 * left (counter-clockwise around a region where it is negative).
 *
 * Edge ids: 2v for the radial edge from vertex v = i spokes + j outwards to
 * v + spokes, and 2v + 1 for the edge from v to the next spoke on the same ring.
 */
function traceContour(mesh: Mesh, val: Float64Array, field: (x: number, y: number) => number): { curves: Curve[]; extent: number } {
  const { rings, spokes, vx, vy } = mesh;
  const px: number[] = [];
  const py: number[] = [];
  const index = new Map<number, number>();

  /** Index of the crossing on edge `id` from vertex u to w, or -1 if there is none. */
  const cross = (id: number, u: number, w: number): number => {
    const fu = val[u];
    const fw = val[w];
    if (fu < 0 === fw < 0) return -1;
    const known = index.get(id);
    if (known !== undefined) return known;
    const ux = vx[u];
    const uy = vy[u];
    const dx = vx[w] - ux;
    const dy = vy[w] - uy;
    let ta = 0;
    let fa = fu;
    let tb = 1;
    let fb = fw;
    let t = fa / (fa - fb);
    let side = 0;
    for (let it = 0; it < 40; it++) {
      const prev = t;
      t = (ta * fb - tb * fa) / (fb - fa);
      const f = field(ux + t * dx, uy + t * dy);
      if (f === 0 || (it > 0 && Math.abs(t - prev) < 1e-14)) break;
      if (f < 0 === fa < 0) {
        ta = t;
        fa = f;
        if (side === -1) fb *= 0.5;
        side = -1;
      } else {
        tb = t;
        fb = f;
        if (side === 1) fa *= 0.5;
        side = 1;
      }
      if (tb - ta < 1e-13) break;
    }
    const k = px.length;
    px.push(ux + t * dx);
    py.push(uy + t * dy);
    index.set(id, k);
    return k;
  };

  // Segments as pairs of crossing indices.
  const segA: number[] = [];
  const segB: number[] = [];
  const found: number[] = [];
  for (let i = 0; i < rings - 1; i++) {
    for (let j = 0; j < spokes; j++) {
      const j2 = (j + 1) % spokes;
      const a = i * spokes + j;
      const b = i * spokes + j2;
      const c = (i + 1) * spokes + j2;
      const d = (i + 1) * spokes + j;
      found.length = 0;
      // Round the quad: a-b along the ring, b-c out along a spoke, c-d back along
      // the next ring, d-a in along the first spoke.
      for (const k of [cross(2 * a + 1, a, b), cross(2 * b, b, c), cross(2 * d + 1, d, c), cross(2 * a, a, d)]) {
        if (k >= 0) found.push(k);
      }
      for (let k = 0; k + 1 < found.length; k += 2) {
        segA.push(found[k]);
        segB.push(found[k + 1]);
      }
    }
  }

  // Which segments meet at each crossing (at most two).
  const n = px.length;
  const s1 = new Int32Array(n).fill(-1);
  const s2 = new Int32Array(n).fill(-1);
  for (let s = 0; s < segA.length; s++) {
    for (const p of [segA[s], segB[s]]) {
      if (s1[p] < 0) s1[p] = s;
      else s2[p] = s;
    }
  }
  const otherSeg = (p: number, s: number) => (s1[p] === s ? s2[p] : s1[p]);
  const otherEnd = (s: number, p: number) => (segA[s] === p ? segB[s] : segA[s]);

  const used = new Uint8Array(segA.length);
  const curves: Curve[] = [];
  let extent = 0;
  for (let s0 = 0; s0 < segA.length; s0++) {
    if (used[s0]) continue;
    used[s0] = 1;
    const chain = [segA[s0], segB[s0]];
    let closed = false;
    // Forwards from the B end.
    let p = segB[s0];
    let s = s0;
    for (;;) {
      const next = otherSeg(p, s);
      if (next < 0) break;
      if (next === s0) {
        closed = true;
        break;
      }
      if (used[next]) break;
      used[next] = 1;
      p = otherEnd(next, p);
      s = next;
      chain.push(p);
    }
    if (closed) {
      // The last crossing reached is the chain's first; do not repeat it.
      if (chain.length > 2 && chain[chain.length - 1] === chain[0]) chain.pop();
    } else {
      // Backwards from the A end.
      const back: number[] = [];
      p = segA[s0];
      s = s0;
      for (;;) {
        const next = otherSeg(p, s);
        if (next < 0 || used[next]) break;
        used[next] = 1;
        p = otherEnd(next, p);
        s = next;
        back.push(p);
      }
      chain.unshift(...back.reverse());
    }

    // Orientation: probe the field on either side of the longest segment.
    let best = 0;
    let bestLen = -1;
    for (let k = 0; k + 1 < chain.length; k++) {
      const len = Math.hypot(px[chain[k + 1]] - px[chain[k]], py[chain[k + 1]] - py[chain[k]]);
      if (len > bestLen) {
        bestLen = len;
        best = k;
      }
    }
    if (bestLen > 0) {
      const ax = px[chain[best]];
      const ay = py[chain[best]];
      const ex = px[chain[best + 1]] - ax;
      const ey = py[chain[best + 1]] - ay;
      const mx = ax + 0.5 * ex;
      const my = ay + 0.5 * ey;
      // Left normal, a quarter of the segment long.
      const nx = -0.25 * ey;
      const ny = 0.25 * ex;
      if (field(mx + nx, my + ny) > field(mx - nx, my - ny)) chain.reverse();
    }

    const points = new Float64Array(2 * chain.length);
    chain.forEach((k, m) => {
      points[2 * m] = px[k];
      points[2 * m + 1] = py[k];
      extent = Math.max(extent, Math.hypot(px[k], py[k]));
    });
    curves.push({ points, closed });
  }
  return { curves, extent };
}

/**
 * Area of the points with a nonzero winding number about a set of polylines
 * (each treated as closed), by scan lines. See analyse() for the method.
 */
function windingArea(curves: Curve[], rows: number): number {
  let ymin = Infinity;
  let ymax = -Infinity;
  for (const c of curves) {
    const p = c.points;
    for (let s = 1; s < p.length; s += 2) {
      ymin = Math.min(ymin, p[s]);
      ymax = Math.max(ymax, p[s]);
    }
  }
  if (!(ymax > ymin)) return 0;
  const dy = (ymax - ymin) / rows;
  // Scan line k is at y = ymin + (k + 1/2) dy. A segment from y1 to y2 crosses
  // the lines with min(y1, y2) <= y < max(y1, y2) (half open, so a vertex shared
  // by two segments is counted once).
  const firstRow = (y: number) => Math.max(0, Math.ceil((y - ymin) / dy - 0.5));
  const forEachEdge = (fn: (x1: number, y1: number, x2: number, y2: number) => void) => {
    for (const c of curves) {
      const p = c.points;
      const m = p.length / 2;
      for (let k = 0; k < m; k++) {
        const k2 = (k + 1) % m;
        fn(p[2 * k], p[2 * k + 1], p[2 * k2], p[2 * k2 + 1]);
      }
    }
  };

  // Count crossings per line, then fill (compressed rows).
  const start = new Int32Array(rows + 1);
  forEachEdge((_x1, y1, _x2, y2) => {
    if (y1 === y2) return;
    const k0 = firstRow(Math.min(y1, y2));
    const k1 = Math.min(rows, firstRow(Math.max(y1, y2)));
    for (let k = k0; k < k1; k++) start[k + 1]++;
  });
  for (let k = 0; k < rows; k++) start[k + 1] += start[k];
  const xs = new Float64Array(start[rows]);
  const ws = new Int8Array(start[rows]);
  const fill = start.slice(0, rows);
  forEachEdge((x1, y1, x2, y2) => {
    if (y1 === y2) return;
    const k0 = firstRow(Math.min(y1, y2));
    const k1 = Math.min(rows, firstRow(Math.max(y1, y2)));
    const slope = (x2 - x1) / (y2 - y1);
    const w = y2 > y1 ? 1 : -1;
    for (let k = k0; k < k1; k++) {
      const y = ymin + (k + 0.5) * dy;
      const q = fill[k]++;
      xs[q] = x1 + (y - y1) * slope;
      ws[q] = w;
    }
  });

  let area = 0;
  for (let k = 0; k < rows; k++) {
    const a = start[k];
    const b = start[k + 1];
    // Insertion sort by x; a line crosses only a handful of segments.
    for (let i = a + 1; i < b; i++) {
      const x = xs[i];
      const w = ws[i];
      let j = i - 1;
      while (j >= a && xs[j] > x) {
        xs[j + 1] = xs[j];
        ws[j + 1] = ws[j];
        j--;
      }
      xs[j + 1] = x;
      ws[j + 1] = w;
    }
    let wind = 0;
    for (let i = a; i < b; i++) {
      if (wind !== 0) area += xs[i] - xs[i - 1];
      wind += ws[i];
    }
  }
  return area * dy;
}

/**
 * Adds 1 to every grid cell whose centre lies inside the triangle (a, b, c),
 * whatever its orientation. Degenerate triangles cover nothing.
 */
function rasterTriangle(
  counts: Uint8Array,
  nx: number,
  ny: number,
  x0: number,
  y0: number,
  cell: number,
  ax: number,
  ay: number,
  bx: number,
  by: number,
  cx: number,
  cy: number,
): void {
  // Cell index range of the triangle's bounding box (cell centres at x0 + (i + 1/2) cell).
  const i0 = Math.max(0, Math.ceil((Math.min(ax, bx, cx) - x0) / cell - 0.5));
  const i1 = Math.min(nx - 1, Math.floor((Math.max(ax, bx, cx) - x0) / cell - 0.5));
  if (i0 > i1) return;
  const j0 = Math.max(0, Math.ceil((Math.min(ay, by, cy) - y0) / cell - 0.5));
  const j1 = Math.min(ny - 1, Math.floor((Math.max(ay, by, cy) - y0) / cell - 0.5));
  if (j0 > j1) return;
  const area = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
  if (area === 0) return;
  const s = area > 0 ? 1 : -1;
  for (let j = j0; j <= j1; j++) {
    const py = y0 + (j + 0.5) * cell;
    for (let i = i0; i <= i1; i++) {
      const px = x0 + (i + 0.5) * cell;
      // Edge functions, all non-negative (after orienting) inside the triangle.
      const w0 = s * ((bx - ax) * (py - ay) - (by - ay) * (px - ax));
      if (w0 < 0) continue;
      const w1 = s * ((cx - bx) * (py - by) - (cy - by) * (px - bx));
      if (w1 < 0) continue;
      const w2 = s * ((ax - cx) * (py - cy) - (ay - cy) * (px - cx));
      if (w2 < 0) continue;
      const p = j * nx + i;
      if (counts[p] < 255) counts[p]++;
    }
  }
}

// --- views ------------------------------------------------------------------------

/**
 * Field-of-view half-widths, arcsec. The page's views snap to these so the
 * curves are seen to grow and shrink, rather than the frame following them.
 */
export const VIEW_STEPS = [0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80];
const snapView = (v: number) => VIEW_STEPS.find((s) => s >= v) ?? VIEW_STEPS[VIEW_STEPS.length - 1];

export interface Views {
  /** Half-widths of the source-plane and image-plane views, arcsec. */
  srcHalf: number;
  imgHalf: number;
  /** How far the source may go from the member along x and y, arcsec, so that all its images stay in the image-plane view. */
  limX: number;
  limY: number;
}

/**
 * The two views for a lens, and the box the source is kept in.
 *
 * Source view: the caustics (and cut) with a margin. Image view: large enough to
 * hold every image of every source in the source view (the outermost image along
 * each axis, with room for the source's width, srcWidth its standard deviation),
 * but no more than four times the critical curves, so they stay visible. Where
 * the view is smaller than that, the source is kept to the box whose images all
 * fit: a source with |beta_x| <= (m11 - kbar(R)) R - margin has every image at
 * |theta_x| <= R (see outerImageRadius), and likewise in y.
 *
 * When the galaxy is only just supercritical that box can shrink to a line, so
 * while either side is below a quarter of the source view the image view is
 * widened by one step.
 */
export function viewsFor(lens: Lens, an: Analysis, srcWidth: number): Views {
  const crit = Math.max(an.criticalExtent[0], an.criticalExtent[1]);
  const caus = Math.max(an.causticExtent[0], an.causticExtent[1]);
  const srcHalf = snapView(Math.max(1.4 * caus, 0.6 * lens.b, 0.25));
  const margin = 2 * srcWidth;
  const reach = srcHalf + margin;
  const need = 1.02 * Math.max(outerImageRadius(lens, reach, 0), outerImageRadius(lens, reach, 1));
  let imgHalf = snapView(Math.max(crit > 0 ? Math.min(need, 4 * crit) : need, 1.3 * crit, 0.25));
  const limits = () => {
    const R = 0.98 * imgHalf;
    const kbarR = meanConvergence(lens, R);
    return [
      Math.max(0, Math.min(srcHalf, (lens.m11 - kbarR) * R - margin)),
      Math.max(0, Math.min(srcHalf, (lens.m22 - kbarR) * R - margin)),
    ];
  };
  let [limX, limY] = limits();
  const largest = VIEW_STEPS[VIEW_STEPS.length - 1];
  while (Math.min(limX, limY) < 0.25 * srcHalf && imgHalf < largest) {
    imgHalf = VIEW_STEPS[VIEW_STEPS.indexOf(imgHalf) + 1];
    [limX, limY] = limits();
  }
  return { srcHalf, imgHalf, limX, limY };
}
