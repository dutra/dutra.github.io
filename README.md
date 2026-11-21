# dutra.cc

Personal site of Isaque Dutra, built with [Astro](https://astro.build) as a fully static, single-page site.

## Commands

| Command           | Action                                           |
| :---------------- | :----------------------------------------------- |
| `npm install`     | Install dependencies                             |
| `npm run dev`     | Start the dev server at `localhost:4321`         |
| `npm run check`   | Type-check the project                           |
| `npm run build`   | Build the static site into `dist/`               |
| `npm run preview` | Serve the built site locally                     |

## Where things live

- `src/data/site.ts` — name, affiliation, email, links, navigation, and the About and Research text.
- `src/content/publications/` — one Markdown file per paper (front matter only). `first_author: true` puts a paper in the First-author group; `date` is the year.
- `src/pages/index.astro` — the page: hero, About, Research, Publications, Contact.
- `src/components/LensHero.astro` — the hero. A background galaxy follows the pointer and its lensed images are added to the Webb photograph by inverse ray shooting (adapted from the kappaccino landing page).
- `src/pages/projects.astro` — the Projects page at `/projects/`. Two card grids: Projects (research code on GitHub, the `software` list at the top of the file) and Experiments (the interactive pages on this site, the `projects` list).
- `src/pages/projects/*.astro` — one page per interactive tool, each with its physics in a DOM-free module under `src/lib/`:
  - `ggsl.astro` + `src/lib/ggsl.ts` — galaxy-galaxy strong lensing: one member galaxy (dPIE) in the cluster's convergence and shear.
  - `sidm.astro` + `src/lib/sidm.ts` — SIDM core formation and collapse with the parametric model of Yang et al. (arXiv:2305.16176), and the halo's lensing strength.
  - `dark-energy.astro` + `src/lib/cpl.ts` — expansion history, distances and lensing quantities for ΛCDM, wCDM and w0–wa models, with the Planck 2018 and DESI DR2 fits marked. The fits table at the top of the page script is the one place those numbers live.
  - `agn-hubble.astro` + `src/lib/agn.ts`, `src/lib/agn-worker.ts` (and `cpl.ts`) — a toy AGN-variability Hubble diagram: damped-random-walk light curves, a maximum-likelihood fit of σ and τ per AGN run in a Web Worker, and distances from toy relations. The relations are illustrative.
- `src/lib/plot.ts` — the small canvas chart all the tool pages use; it takes its colours from the CSS variables, so charts follow the theme.
- `src/styles/tools.css` — shared layout and control styles for the tool pages (`.tool-*` classes).
- `src/lib/controls.ts` — shared slider behaviour for the tool pages: `setFrozen` (a parameter the selected model fixes is disabled, muted and tagged "fixed") and `setMarks` (clickable reference ticks under a slider, such as the Planck and DESI fits).
- `src/assets/projects/` — SVG thumbnails for the cards.
- `src/pages/lens.astro` and `src/components/LensLab.astro` — the simulator at `/lens/`: several background galaxies at their own redshifts, an adjustable mass model, critical curves and caustics, convergence and magnification maps, and readouts. The component holds only the controls and pointer handling; each control is bound by its `data-*` attribute.
- `src/lib/lens/` — the physics, shared by the hero and the simulator, with no page code in it:
  - `model.ts` — the mass model in physical units (velocity dispersions in km/s, radii in arcsec): one or two cored isothermal ellipsoids for the cluster and a truncated isothermal sphere per member galaxy, scaled from its brightness. `DEFAULT_MODEL` is tuned by eye to the photograph, not fitted. The image scale (126″ across) comes from ESA's stated field of view.
  - `cosmology.ts` — flat ΛCDM distances (H0 = 70, Ωm = 0.3), the lensing efficiency D_ls/D_s, and the mass inside an Einstein radius.
  - `field.ts` — the deflection tabulated on a grid over the whole photograph, plus what is read from it: magnification, Einstein radius, total magnification of a source.
  - `trace.ts` — critical curves and caustics from that grid (marching squares on the sign of the Jacobian determinant).
  - `renderer.ts` — WebGL 2 drawing: per pixel, the lens equation for each source plane and the brightness of the galaxy the ray lands on.
- `src/lib/lens.ts` — the hero's view of that model: the default parameters and one source at z = 2.
- `src/lib/curves.ts` — traces the hero's critical curves and caustics on a fine grid at build time and returns them as SVG paths.
- `scripts/find-members.mjs` — measures the member galaxies on the image and writes `src/data/lens-members.json`. Rerun it (`node scripts/find-members.mjs`) if the image or the seed list changes.
- `src/assets/abell-s1063-webb.jpg` — the hero backdrop, a crop of ESA/Webb image [potm2505a](https://esawebb.org/images/potm2505a/) (CC BY 4.0; credit ESA/Webb, NASA & CSA, H. Atek, M. Zamani; acknowledgement R. Endsley). `public/og.jpg` is the link-preview crop of the same image.
- `src/assets/portrait.jpeg` — the portrait; Astro resizes and converts it at build time.
- `src/styles/global.css` — colours, type and shared styles.
- `src/lib/moon.ts` — the footer moon phase, computed in the browser from the date.

Deploy `dist/` to any static host.
