// The lens as the hero sees it: the default cluster model (src/lib/lens/model.ts)
// acting on a background galaxy at a fixed redshift (other redshifts are a
// rescaling, see efficiencyRatio). The simulator page
// uses the same model with everything adjustable.
//
// Angles are in units of the image width, measured from the brightest cluster
// galaxy, x to the right and y down.
import { efficiency } from './lens/cosmology';
import { DEFAULT_MODEL, Z_LENS, makeDeflector } from './lens/model';

export { BCG, IMAGE_ASPECT } from './lens/model';

/** Half-width of the square the lens is drawn over, centred on the BCG. 0.6 covers the whole photograph. */
export const HALF = 0.6;

/** Redshift of the hero's background galaxy. */
export const Z_SOURCE = 2;

const eff = efficiency(Z_LENS, Z_SOURCE);
const model = makeDeflector(DEFAULT_MODEL);

/** Deflection angle at (x, y), written into `out` as [ax, ay]. */
export function deflect(x: number, y: number, out: Float64Array | number[]): void {
  model(x, y, out);
  out[0] *= eff;
  out[1] *= eff;
}

/** Factor that turns `deflect()` into the deflection for a source at redshift z instead of Z_SOURCE. */
export function efficiencyRatio(z: number): number {
  return efficiency(Z_LENS, z) / eff;
}
