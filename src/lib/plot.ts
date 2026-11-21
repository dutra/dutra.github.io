// A small canvas chart used by the project pages: axes with ticks, then lines,
// points, bands and labels in data coordinates. It takes its colours and font
// from the page's CSS variables, so charts follow the light and dark themes.
//
// Typical use:
//   const plot = new Plot(canvas, { x: { min: 0, max: 3, label: 'Redshift' }, y: { min: -1, max: 1, label: 'Δμ' } });
//   plot.observe(() => {            // runs now, and again on resize or theme change
//     plot.begin();                 // clears and draws the axes
//     plot.line(xs, ys, { color: SERIES[0] });
//   });
//   ...and call plot.redraw() whenever the data changes.

export interface Axis {
  min: number;
  max: number;
  label?: string;
  /** Logarithmic scale; min and max must be positive. */
  log?: boolean;
  /** Tick positions. Defaults to round numbers (decades on a log axis). */
  ticks?: number[];
  format?: (value: number) => string;
}

export interface PlotOptions {
  x: Axis;
  y: Axis;
  margin?: Partial<Record<'left' | 'right' | 'top' | 'bottom', number>>;
}

export interface Stroke {
  color: string;
  width?: number;
  dash?: number[];
  /** 0 to 1. */
  alpha?: number;
}

/** Series colours that read on both the paper and the dark background. Use them in this order. */
export const SERIES = ['#d98a2b', '#3d85d6', '#d9534f', '#2aa198', '#9575cd'] as const;

type Numbers = ArrayLike<number>;

function niceTicks(lo: number, hi: number, target = 6): number[] {
  // An axis may run backwards (min above max, as for magnitudes); the ticks are the same.
  const min = Math.min(lo, hi);
  const max = Math.max(lo, hi);
  const span = max - min;
  if (!(span > 0)) return [min];
  const raw = span / target;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw)!;
  const ticks: number[] = [];
  for (let v = Math.ceil(min / step - 1e-9) * step; v <= max + step * 1e-9; v += step) {
    ticks.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  }
  return ticks;
}

function decadeTicks(lo: number, hi: number): number[] {
  const min = Math.min(lo, hi);
  const max = Math.max(lo, hi);
  const ticks: number[] = [];
  for (let e = Math.ceil(Math.log10(min) - 1e-9); e <= Math.log10(max) + 1e-9; e++) ticks.push(10 ** e);
  return ticks;
}

const SUPERSCRIPT = '⁰¹²³⁴⁵⁶⁷⁸⁹';
const superscript = (n: number) =>
  (n < 0 ? '⁻' : '') + String(Math.abs(n)).replace(/\d/g, (d) => SUPERSCRIPT[Number(d)]);

function defaultFormat(axis: Axis): (v: number) => string {
  if (axis.log) {
    return (v) => {
      const e = Math.round(Math.log10(v));
      return e >= -2 && e <= 3 ? String(Number(v.toPrecision(3))) : `10${superscript(e)}`;
    };
  }
  return (v) => String(Number(v.toPrecision(6))).replace('-', '−');
}

export class Plot {
  x: Axis;
  y: Axis;
  /** Size in CSS pixels, set by begin(). */
  width = 0;
  height = 0;
  /** Colours read from the page's CSS variables, set by begin(). */
  theme = { ink: '#000', soft: '#333', muted: '#777', line: 'rgba(0,0,0,0.14)', strong: 'rgba(0,0,0,0.3)', bg: '#fff', font: 'sans-serif' };

  private ctx: CanvasRenderingContext2D;
  private margin: Record<'left' | 'right' | 'top' | 'bottom', number>;
  private draw: (() => void) | null = null;

  constructor(
    readonly canvas: HTMLCanvasElement,
    options: PlotOptions,
  ) {
    this.ctx = canvas.getContext('2d')!;
    this.x = options.x;
    this.y = options.y;
    this.margin = { left: 58, right: 16, top: 14, bottom: 44, ...options.margin };
  }

  /** Registers the drawing function, runs it, and reruns it when the canvas is resized or the theme changes. */
  observe(draw: () => void): void {
    this.draw = draw;
    new ResizeObserver(() => this.redraw()).observe(this.canvas);
    window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => this.redraw());
    this.redraw();
  }

  /** Runs the drawing function given to observe(). Call it when the data changes. */
  redraw(): void {
    this.draw?.();
  }

  // --- coordinates ---------------------------------------------------------------
  private fraction(axis: Axis, v: number): number {
    return axis.log
      ? (Math.log(v) - Math.log(axis.min)) / (Math.log(axis.max) - Math.log(axis.min))
      : (v - axis.min) / (axis.max - axis.min);
  }
  /** Data to CSS pixels. */
  px(x: number): number {
    return this.margin.left + this.fraction(this.x, x) * (this.width - this.margin.left - this.margin.right);
  }
  py(y: number): number {
    return this.height - this.margin.bottom - this.fraction(this.y, y) * (this.height - this.margin.top - this.margin.bottom);
  }
  /** CSS pixels (relative to the canvas) to data, for pointer handling. */
  toData(px: number, py: number): [number, number] {
    const fx = (px - this.margin.left) / (this.width - this.margin.left - this.margin.right);
    const fy = (this.height - this.margin.bottom - py) / (this.height - this.margin.top - this.margin.bottom);
    const un = (axis: Axis, f: number) =>
      axis.log ? Math.exp(Math.log(axis.min) + f * (Math.log(axis.max) - Math.log(axis.min))) : axis.min + f * (axis.max - axis.min);
    return [un(this.x, fx), un(this.y, fy)];
  }

  // --- frame -----------------------------------------------------------------------
  /** Sizes the canvas to its box, clears it and draws the grid, ticks and axis labels. */
  begin(): void {
    const { canvas, ctx } = this;
    const rect = canvas.getBoundingClientRect();
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    this.width = rect.width;
    this.height = rect.height;
    canvas.width = Math.max(1, Math.round(rect.width * dpr));
    canvas.height = Math.max(1, Math.round(rect.height * dpr));
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, rect.width, rect.height);

    const css = getComputedStyle(canvas);
    const v = (name: string, fallback: string) => css.getPropertyValue(name).trim() || fallback;
    this.theme = {
      ink: v('--ink', '#111'),
      soft: v('--ink-soft', '#333'),
      muted: v('--muted', '#777'),
      line: v('--line', 'rgba(0,0,0,0.14)'),
      strong: v('--line-strong', 'rgba(0,0,0,0.3)'),
      bg: v('--bg-raised', '#fff'),
      font: v('--font-ui', 'system-ui, sans-serif'),
    };

    const x0 = this.margin.left;
    const x1 = this.width - this.margin.right;
    const y0 = this.margin.top;
    const y1 = this.height - this.margin.bottom;
    const ticks = (axis: Axis) =>
      (axis.ticks ?? (axis.log ? decadeTicks(axis.min, axis.max) : niceTicks(axis.min, axis.max))).filter(
        (t) => t >= Math.min(axis.min, axis.max) && t <= Math.max(axis.min, axis.max),
      );

    ctx.lineWidth = 1;
    ctx.font = `500 11.5px ${this.theme.font}`;
    ctx.fillStyle = this.theme.muted;
    ctx.strokeStyle = this.theme.line;
    ctx.setLineDash([]);
    ctx.globalAlpha = 1;

    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    const fx = this.x.format ?? defaultFormat(this.x);
    for (const t of ticks(this.x)) {
      const p = Math.round(this.px(t)) + 0.5;
      ctx.beginPath();
      ctx.moveTo(p, y0);
      ctx.lineTo(p, y1);
      ctx.stroke();
      ctx.fillText(fx(t), p, y1 + 7);
    }
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    const fy = this.y.format ?? defaultFormat(this.y);
    for (const t of ticks(this.y)) {
      const p = Math.round(this.py(t)) + 0.5;
      ctx.beginPath();
      ctx.moveTo(x0, p);
      ctx.lineTo(x1, p);
      ctx.stroke();
      ctx.fillText(fy(t), x0 - 8, p);
    }

    ctx.strokeStyle = this.theme.strong;
    ctx.strokeRect(Math.round(x0) + 0.5, Math.round(y0) + 0.5, Math.round(x1 - x0), Math.round(y1 - y0));

    ctx.fillStyle = this.theme.soft;
    ctx.font = `600 12px ${this.theme.font}`;
    if (this.x.label) {
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';
      ctx.fillText(this.x.label, (x0 + x1) / 2, this.height - 6);
    }
    if (this.y.label) {
      ctx.save();
      ctx.translate(13, (y0 + y1) / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText(this.y.label, 0, 0);
      ctx.restore();
    }
  }

  /** Runs `body` with drawing clipped to the plot area. */
  private clipped(body: () => void): void {
    const { ctx } = this;
    ctx.save();
    ctx.beginPath();
    ctx.rect(this.margin.left, this.margin.top, this.width - this.margin.left - this.margin.right, this.height - this.margin.top - this.margin.bottom);
    ctx.clip();
    body();
    ctx.restore();
  }

  private usable(x: number, y: number): boolean {
    return Number.isFinite(x) && Number.isFinite(y) && (!this.x.log || x > 0) && (!this.y.log || y > 0);
  }

  // --- marks ---------------------------------------------------------------------
  /** A polyline through (xs[i], ys[i]). Points that are not finite (or not positive on a log axis) break the line. */
  line(xs: Numbers, ys: Numbers, stroke: Stroke): void {
    this.clipped(() => {
      const { ctx } = this;
      ctx.strokeStyle = stroke.color;
      ctx.lineWidth = stroke.width ?? 2;
      ctx.globalAlpha = stroke.alpha ?? 1;
      ctx.setLineDash(stroke.dash ?? []);
      ctx.lineJoin = 'round';
      ctx.beginPath();
      let pen = false;
      for (let i = 0; i < xs.length; i++) {
        if (!this.usable(xs[i], ys[i])) {
          pen = false;
          continue;
        }
        if (pen) ctx.lineTo(this.px(xs[i]), this.py(ys[i]));
        else ctx.moveTo(this.px(xs[i]), this.py(ys[i]));
        pen = true;
      }
      ctx.stroke();
    });
  }

  /** Dots at (xs[i], ys[i]), with optional vertical error bars of half-length yerr[i]. */
  points(xs: Numbers, ys: Numbers, style: { color: string; radius?: number; alpha?: number; yerr?: Numbers }): void {
    this.clipped(() => {
      const { ctx } = this;
      const r = style.radius ?? 3;
      ctx.fillStyle = ctx.strokeStyle = style.color;
      ctx.globalAlpha = style.alpha ?? 1;
      ctx.lineWidth = 1;
      ctx.setLineDash([]);
      for (let i = 0; i < xs.length; i++) {
        if (!this.usable(xs[i], ys[i])) continue;
        const px = this.px(xs[i]);
        if (style.yerr) {
          ctx.beginPath();
          ctx.moveTo(px, this.py(ys[i] - style.yerr[i]));
          ctx.lineTo(px, this.py(ys[i] + style.yerr[i]));
          ctx.stroke();
        }
        ctx.beginPath();
        ctx.arc(px, this.py(ys[i]), r, 0, 2 * Math.PI);
        ctx.fill();
      }
    });
  }

  /** A filled band between lo[i] and hi[i]. */
  band(xs: Numbers, lo: Numbers, hi: Numbers, style: { color: string; alpha?: number }): void {
    this.clipped(() => {
      const { ctx } = this;
      ctx.fillStyle = style.color;
      ctx.globalAlpha = style.alpha ?? 0.18;
      ctx.beginPath();
      for (let i = 0; i < xs.length; i++) ctx.lineTo(this.px(xs[i]), this.py(hi[i]));
      for (let i = xs.length - 1; i >= 0; i--) ctx.lineTo(this.px(xs[i]), this.py(lo[i]));
      ctx.closePath();
      ctx.fill();
    });
  }

  /** A vertical or horizontal rule across the plot area at a data value. */
  rule(direction: 'x' | 'y', value: number, stroke: Stroke): void {
    if (direction === 'x') this.line([value, value], [this.y.min, this.y.max], { width: 1.2, ...stroke });
    else this.line([this.x.min, this.x.max], [value, value], { width: 1.2, ...stroke });
  }

  /** Text anchored at a data position, offset by (dx, dy) CSS pixels. */
  text(x: number, y: number, label: string, style: { color?: string; align?: CanvasTextAlign; dx?: number; dy?: number; bold?: boolean } = {}): void {
    const { ctx } = this;
    ctx.globalAlpha = 1;
    ctx.fillStyle = style.color ?? this.theme.soft;
    ctx.font = `${style.bold ? 600 : 500} 12px ${this.theme.font}`;
    ctx.textAlign = style.align ?? 'left';
    ctx.textBaseline = 'middle';
    ctx.fillText(label, this.px(x) + (style.dx ?? 0), this.py(y) + (style.dy ?? 0));
  }

  /** A legend in a corner of the plot area: one coloured swatch and label per entry. */
  legend(entries: { label: string; color: string; dash?: number[] }[], corner: 'top-left' | 'top-right' | 'bottom-left' | 'bottom-right' = 'top-right'): void {
    const { ctx } = this;
    ctx.font = `500 12px ${this.theme.font}`;
    ctx.textBaseline = 'middle';
    ctx.textAlign = 'left';
    const row = 18;
    const widest = Math.max(...entries.map((e) => ctx.measureText(e.label).width));
    const left = corner.endsWith('left') ? this.margin.left + 12 : this.width - this.margin.right - 12 - widest - 28;
    const top = corner.startsWith('top') ? this.margin.top + 14 : this.height - this.margin.bottom - 14 - row * (entries.length - 1);
    // A translucent backdrop in the page colour keeps the labels readable over data.
    ctx.globalAlpha = 0.78;
    ctx.fillStyle = this.theme.bg;
    ctx.fillRect(left - 8, top - row / 2 - 3, widest + 44, row * entries.length + 6);
    entries.forEach((e, i) => {
      const y = top + i * row;
      ctx.globalAlpha = 1;
      ctx.strokeStyle = e.color;
      ctx.lineWidth = 2;
      ctx.setLineDash(e.dash ?? []);
      ctx.beginPath();
      ctx.moveTo(left, y);
      ctx.lineTo(left + 20, y);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = this.theme.soft;
      ctx.fillText(e.label, left + 28, y);
    });
  }
}
