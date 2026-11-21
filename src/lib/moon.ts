/** Moon phase from the date alone: days since a known new moon, modulo the synodic month. */
const SYNODIC = 29.530588853; // days
const REF_NEW_MOON = Date.UTC(2000, 0, 6, 18, 14); // 2000-01-06 18:14 UTC

const NAMES = [
  'new moon', 'waxing crescent', 'first quarter', 'waxing gibbous',
  'full moon', 'waning gibbous', 'last quarter', 'waning crescent',
];

export function moonPhase(date: Date = new Date()) {
  const days = (date.getTime() - REF_NEW_MOON) / 86_400_000;
  const phase = (((days / SYNODIC) % 1) + 1) % 1; // 0 = new, 0.5 = full
  const lit = Math.round(((1 - Math.cos(2 * Math.PI * phase)) / 2) * 100);
  const name = NAMES[Math.floor(((phase + 1 / 16) % 1) * 8)];
  const glyph = phase < 0.5 ? '☽' : '☾';
  return { phase, lit, name, label: `${glyph} ${name} · ${lit}% lit` };
}
