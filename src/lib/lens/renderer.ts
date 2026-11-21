// Draws the lensed sky with WebGL 2. For every pixel the fragment shader looks
// up the deflection there, follows the ray back to each source plane and reads
// the brightness of the galaxy it lands on: surface brightness is conserved by
// lensing, so that is the brightness of the image. No lensing physics lives in
// the shader beyond the lens equation; the deflection comes from field.ts.
import type { LensField } from './field';

export const MAX_SOURCES = 4;

export interface SourceView {
  /** Position in the source plane, width units from the BCG. */
  x: number;
  y: number;
  /** Lensing efficiency D_ls/D_s at the source's redshift. */
  eff: number;
  /** Scale radius, width units. */
  size: number;
  /** Axis ratio and orientation (radians) of the disc. */
  q: number;
  angle: number;
  /** Picks the arrangement of spiral arms and star-forming knots. */
  seed: number;
  disc: [number, number, number];
  knot: [number, number, number];
}

export interface Frame {
  /** Canvas size in CSS pixels. */
  width: number;
  height: number;
  /** Where the BCG sits on the canvas and how many CSS pixels one image width spans. */
  originX: number;
  originY: number;
  scale: number;
  sources: SourceView[];
  /**
   * 0: lensed galaxies only, as light on black for the page to blend with `mix-blend-mode: screen`.
   * 1: convergence map. 2: magnification map; both translucent, to be laid over the page normally.
   */
  map: 0 | 1 | 2;
  /** Lensing efficiency the map is drawn for. */
  mapEff: number;
}

const VERTEX = `#version 300 es
in vec2 aPos;
uniform vec2 uCanvas;
uniform vec3 uView;
out vec2 vTheta;
void main() {
  vec2 px = vec2(aPos.x * 0.5 + 0.5, 0.5 - aPos.y * 0.5) * uCanvas;
  vTheta = (px - uView.xy) / uView.z;
  gl_Position = vec4(aPos, 0.0, 1.0);
}`;

const FRAGMENT = `#version 300 es
precision highp float;
precision highp sampler2D;

uniform sampler2D uField;
uniform vec2 uGridOrigin;
uniform vec2 uGridSize;
uniform float uGridStep;

uniform int uCount;
uniform vec4 uSrcA[${MAX_SOURCES}]; // x, y, efficiency, size
uniform vec4 uSrcB[${MAX_SOURCES}]; // axis ratio, angle, seed, unused
uniform vec3 uDisc[${MAX_SOURCES}];
uniform vec3 uKnot[${MAX_SOURCES}];

uniform int uMap;
uniform float uMapEff;

in vec2 vTheta;
out vec4 frag;

// Deflection for D_ls/D_s = 1, interpolated between grid nodes.
vec2 alphaAt(vec2 theta) {
  vec2 g = clamp((theta - uGridOrigin) / uGridStep, vec2(0.0), uGridSize - 1.001);
  ivec2 i = ivec2(g);
  vec2 f = g - vec2(i);
  vec2 a00 = texelFetch(uField, i, 0).rg;
  vec2 a10 = texelFetch(uField, i + ivec2(1, 0), 0).rg;
  vec2 a01 = texelFetch(uField, i + ivec2(0, 1), 0).rg;
  vec2 a11 = texelFetch(uField, i + ivec2(1, 1), 0).rg;
  return mix(mix(a00, a10, f.x), mix(a01, a11, f.x), f.y);
}

float hash(float n) {
  return fract(sin(n) * 43758.5453);
}

// A late-type galaxy: exponential disc with two spiral arms, a warm bulge and
// a handful of star-forming knots. d is the offset from its centre.
vec3 galaxy(vec2 d, vec4 a, vec4 b, vec3 disc, vec3 knot) {
  float c = cos(b.y);
  float s = sin(b.y);
  vec2 p = vec2(c * d.x + s * d.y, (c * d.y - s * d.x) / b.x) / a.w;
  float r2 = dot(p, p);
  if (r2 > 40.0) return vec3(0.0);
  float r = sqrt(r2);
  float arms = 1.0 + 0.6 * cos(2.0 * atan(p.y, p.x) - 4.5 * log(r + 0.1) + b.z);
  vec3 col = disc * (1.25 * exp(-1.7 * r) * arms);
  col += mix(disc, vec3(1.0, 0.93, 0.8), 0.7) * (0.9 * exp(-r2 / 0.08));
  for (int k = 0; k < 5; k++) {
    float fk = float(k) + b.z * 7.31;
    float kr = 0.35 + 1.3 * hash(fk);
    float ka = 6.2832 * hash(fk + 17.0);
    vec2 e = p - kr * vec2(cos(ka), sin(ka));
    col += knot * (exp(-dot(e, e) / 0.03) * (0.8 + hash(fk + 31.0)) * exp(-0.5 * kr));
  }
  return col;
}

vec3 ramp(float t) {
  const vec3 c0 = vec3(0.0, 0.0, 0.03);
  const vec3 c1 = vec3(0.23, 0.06, 0.44);
  const vec3 c2 = vec3(0.72, 0.21, 0.47);
  const vec3 c3 = vec3(0.99, 0.55, 0.38);
  const vec3 c4 = vec3(0.99, 0.99, 0.75);
  t *= 4.0;
  if (t < 1.0) return mix(c0, c1, t);
  if (t < 2.0) return mix(c1, c2, t - 1.0);
  if (t < 3.0) return mix(c2, c3, t - 2.0);
  return mix(c3, c4, t - 3.0);
}

void main() {
  vec2 alpha = alphaAt(vTheta);
  vec3 light = vec3(0.0);
  for (int k = 0; k < ${MAX_SOURCES}; k++) {
    if (k >= uCount) break;
    // The lens equation: beta = theta - (D_ls/D_s) alpha.
    vec2 beta = vTheta - uSrcA[k].z * alpha;
    light += galaxy(beta - uSrcA[k].xy, uSrcA[k], uSrcB[k], uDisc[k], uKnot[k]);
  }

  if (uMap == 0) {
    // Light on black: the page blends this canvas over the photograph with "screen".
    frag = vec4(min(light, vec3(1.0)), 1.0);
    return;
  }

  // Jacobian of the lens mapping from differences of the deflection.
  float h = uGridStep;
  vec2 dx = (alphaAt(vTheta + vec2(h, 0.0)) - alphaAt(vTheta - vec2(h, 0.0))) * (uMapEff / (2.0 * h));
  vec2 dy = (alphaAt(vTheta + vec2(0.0, h)) - alphaAt(vTheta - vec2(0.0, h))) * (uMapEff / (2.0 * h));
  float kappa = 0.5 * (dx.x + dy.y);
  float det = (1.0 - dx.x) * (1.0 - dy.y) - dx.y * dy.x;
  float t = uMap == 1
    ? log2(max(kappa, 1e-4) / 0.125) / 4.0          // convergence from 1/8 to 2
    : log(1.0 / max(abs(det), 1e-4)) / log(100.0);  // magnification from 1 to 100
  const float cover = 0.8;
  // A map is laid over the photograph instead, mostly covering it (premultiplied alpha).
  frag = vec4(min(light + ramp(clamp(t, 0.0, 1.0)), vec3(1.0)) * cover, cover);
}`;

export class LensRenderer {
  private gl: WebGL2RenderingContext;
  private texture!: WebGLTexture;
  private uniforms: Record<string, WebGLUniformLocation | null> = {};

  constructor(
    private canvas: HTMLCanvasElement,
    private field: LensField,
  ) {
    const gl = canvas.getContext('webgl2', { alpha: true, premultipliedAlpha: true, antialias: false });
    if (!gl) throw new Error('WebGL 2 is not available');
    this.gl = gl;
    this.setup();
  }

  /** (Re)creates every GL object; also used after the browser restores a lost context. */
  setup(): void {
    const gl = this.gl;
    const compile = (type: number, source: string) => {
      const shader = gl.createShader(type)!;
      gl.shaderSource(shader, source);
      gl.compileShader(shader);
      if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        throw new Error(`Shader: ${gl.getShaderInfoLog(shader)}`);
      }
      return shader;
    };
    const program = gl.createProgram()!;
    gl.attachShader(program, compile(gl.VERTEX_SHADER, VERTEX));
    gl.attachShader(program, compile(gl.FRAGMENT_SHADER, FRAGMENT));
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) throw new Error(`Program: ${gl.getProgramInfoLog(program)}`);
    gl.useProgram(program);

    // One triangle that covers the canvas.
    gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    const loc = gl.getAttribLocation(program, 'aPos');
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);

    this.texture = gl.createTexture()!;
    gl.bindTexture(gl.TEXTURE_2D, this.texture);
    // Float textures can always be read texel by texel; the shader interpolates itself.
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

    for (const name of ['uField', 'uGridOrigin', 'uGridSize', 'uGridStep', 'uCount', 'uSrcA', 'uSrcB', 'uDisc', 'uKnot', 'uMap', 'uMapEff', 'uCanvas', 'uView']) {
      this.uniforms[name] = gl.getUniformLocation(program, name);
    }
    gl.uniform1i(this.uniforms.uField, 0);
    gl.uniform2f(this.uniforms.uGridOrigin, this.field.x0, this.field.y0);
    gl.uniform2f(this.uniforms.uGridSize, this.field.gw, this.field.gh);
    gl.uniform1f(this.uniforms.uGridStep, this.field.h);
    gl.clearColor(0, 0, 0, 0);
    this.uploadField();
  }

  /** Sends the current deflection field to the GPU. Call after LensField.update(). */
  uploadField(): void {
    const gl = this.gl;
    gl.bindTexture(gl.TEXTURE_2D, this.texture);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RG32F, this.field.gw, this.field.gh, 0, gl.RG, gl.FLOAT, this.field.a);
  }

  /** Sets the drawing-buffer size, in device pixels. */
  resize(width: number, height: number): void {
    if (this.canvas.width !== width || this.canvas.height !== height) {
      this.canvas.width = width;
      this.canvas.height = height;
    }
  }

  render(frame: Frame): void {
    const gl = this.gl;
    const u = this.uniforms;
    const n = Math.min(frame.sources.length, MAX_SOURCES);
    const a = new Float32Array(4 * MAX_SOURCES);
    const b = new Float32Array(4 * MAX_SOURCES);
    const disc = new Float32Array(3 * MAX_SOURCES);
    const knot = new Float32Array(3 * MAX_SOURCES);
    for (let k = 0; k < n; k++) {
      const s = frame.sources[k];
      a.set([s.x, s.y, s.eff, s.size], 4 * k);
      b.set([s.q, s.angle, s.seed, 0], 4 * k);
      disc.set(s.disc, 3 * k);
      knot.set(s.knot, 3 * k);
    }
    gl.viewport(0, 0, this.canvas.width, this.canvas.height);
    gl.uniform2f(u.uCanvas, frame.width, frame.height);
    gl.uniform3f(u.uView, frame.originX, frame.originY, frame.scale);
    gl.uniform1i(u.uCount, n);
    gl.uniform4fv(u.uSrcA, a);
    gl.uniform4fv(u.uSrcB, b);
    gl.uniform3fv(u.uDisc, disc);
    gl.uniform3fv(u.uKnot, knot);
    gl.uniform1i(u.uMap, frame.map);
    gl.uniform1f(u.uMapEff, frame.mapEff);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }
}
