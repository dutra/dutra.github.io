// Two small behaviours shared by the sliders on the project pages (styles in
// src/styles/tools.css):
//
//  - Freezing: when the chosen model does not let a parameter vary, its slider
//    is disabled, drawn in a muted colour and tagged "fixed".
//  - Marks: labelled ticks under a slider at reference values (for example a
//    published measurement). Clicking a mark moves the slider there.
//
// Both work on a slider inside a `.tool-row` label:
//   <label class="tool-row"><span>Name <output></output></span><input type="range" …></label>

/** Freezes or releases a slider. A frozen slider keeps showing its value; set that value before freezing. */
export function setFrozen(input: HTMLInputElement, frozen: boolean, reason = 'Fixed in this model'): void {
  const row = input.closest<HTMLElement>('.tool-row');
  input.disabled = frozen;
  if (!row) return;
  row.toggleAttribute('data-frozen', frozen);
  row.title = frozen ? reason : '';
  let tag = row.querySelector<HTMLElement>('.tool-fixed');
  if (frozen && !tag) {
    tag = document.createElement('small');
    tag.className = 'tool-fixed';
    tag.textContent = 'fixed';
    // After the name, before the value.
    row.querySelector('span')?.insertBefore(tag, row.querySelector('output'));
  } else if (!frozen && tag) {
    tag.remove();
  }
}

export interface Mark {
  value: number;
  /** Short text under the tick, e.g. "Planck". */
  label: string;
  /** Longer description shown on hover: the source and the exact value. */
  title?: string;
}

/** Shows reference marks under a slider, replacing any it had. Pass an empty list to remove them. */
export function setMarks(input: HTMLInputElement, marks: Mark[]): void {
  const row = input.closest<HTMLElement>('.tool-row');
  if (!row) return;
  row.querySelector('.tool-marks')?.remove();
  const min = Number(input.min);
  const max = Number(input.max);
  const inside = marks.filter((m) => m.value >= min && m.value <= max);
  if (!inside.length) return;
  const strip = document.createElement('span');
  strip.className = 'tool-marks';
  for (const m of inside) {
    const b = document.createElement('button');
    b.type = 'button';
    b.className = 'tool-mark';
    b.style.setProperty('--at', String((m.value - min) / (max - min)));
    b.textContent = m.label;
    b.title = m.title ?? `${m.label}: ${m.value}`;
    b.addEventListener('click', (ev) => {
      // The row is a <label>: without this the click would also focus the slider.
      ev.preventDefault();
      if (input.disabled) return;
      input.value = String(m.value);
      input.dispatchEvent(new Event('input', { bubbles: true }));
    });
    strip.append(b);
  }
  row.append(strip);
}
