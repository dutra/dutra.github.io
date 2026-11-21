// Web Worker for src/pages/projects/agn-hubble.astro: fits every simulated
// light curve (agn.ts, fitPopulation) off the main thread, so the page stays
// responsive while it works.
//
// Messages in:  { job, params }   start fitting the survey `params`; a newer
//                                 job makes any older one stop at its next slice.
// Messages out: { type: 'rows', job, from, rows }  FIELDS numbers per object,
//                                                  objects from, from + 1, …
//               { type: 'done', job, ms }          all fitted, after ms milliseconds.
import { fitPopulation, makePopulation, type SurveyParams } from './agn';

export interface FitJob {
  job: number;
  params: SurveyParams;
}

export type FitMessage = { type: 'rows'; job: number; from: number; rows: Float64Array } | { type: 'done'; job: number; ms: number };

// The worker global, typed by hand: the project's TypeScript lib is the DOM one.
const scope = self as unknown as {
  onmessage: ((ev: MessageEvent<FitJob>) => void) | null;
  postMessage(message: FitMessage, transfer?: Transferable[]): void;
};

let current = -1;

scope.onmessage = (ev) => {
  const { job, params } = ev.data;
  current = job;
  const start = performance.now();
  const pop = makePopulation(params);
  void fitPopulation(
    pop,
    (from, rows) => scope.postMessage({ type: 'rows', job, from, rows }, [rows.buffer]),
    () => current !== job,
    30,
  ).then((finished) => {
    if (finished) scope.postMessage({ type: 'done', job, ms: performance.now() - start });
  });
};
