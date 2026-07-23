// ---------------------------------------------------------------------------
// ChIP Genome Viewer — front-end
//
// Design notes:
//   * The SVG is rendered ONCE per zoom/data change, spanning the full data
//     range. Its pixel width equals `container.clientWidth * zoom`, where
//     zoom = 1 means "fit all data into the viewport".
//   * Panning is native horizontal scrolling of the container. No re-render
//     happens on scroll, which eliminates the jitter caused by re-computing
//     the view on every scroll event.
//   * Zoom keeps the coordinate under the anchor point (cursor for wheel,
//     viewport center for buttons) fixed on screen.
// ---------------------------------------------------------------------------

const svg = document.getElementById('chart');
const status = document.getElementById('status');
const container = document.getElementById('chart-container');

const MIN_ZOOM = 1;               // 1 = fit entire data range in viewport
const MAX_TOTAL_WIDTH = 2_000_000;  // cap SVG width so browsers stay happy
const ZOOM_STEP = 1.6;
const LABEL_MIN_WIDTH_PX = 20;    // hide bar labels below this rendered width

// Signal (histogram) lane sizing. Signals may hold millions of points, so
// each render is viewport-scoped: we only build path vertices for the bins
// that fall within (visible x range + padding), at roughly one bin per pixel.
// This keeps the DOM small regardless of zoom or track size.
const SIGNAL_LANE_HEIGHT = 100;
const SIGNAL_LANE_GAP = 10;
const SIGNAL_BINS_PER_PX = 1;      // aggregation density inside the viewport
const SIGNAL_MAX_BINS_PER_LANE = 4000; // safety cap per lane per render
const SIGNAL_VIEWPORT_PAD_FRAC = 0.75; // extra viewport-widths rendered off-screen

const state = {
  genes: [],
  peaks: [],
  signals: [],
  dataMin: 0,
  dataMax: 1,
  zoom: 1,
};

// Saved geometry from the last render() so we can re-render only the signal
// lanes on scroll without recomputing the full layout.
const layoutState = {
  margin: null,
  plotWidth: 0,
  signalYs: [],
};

// ---------------------------------------------------------------------------
// Data loading
// ---------------------------------------------------------------------------

async function loadDataFromFiles() {
  const genesFile = document.getElementById('genes-file').files[0];
  const peaksFile = document.getElementById('peaks-file').files[0];
  const signalFiles = Array.from(
    document.getElementById('signals-file').files || [],
  );

  if (!genesFile || !peaksFile) {
    status.textContent = 'Please select both files first.';
    return;
  }

  try {
    const [genesText, peaksText, ...signalTexts] = await Promise.all([
      readFileSmart(genesFile),
      readFileSmart(peaksFile),
      ...signalFiles.map(readFileSmart),
    ]);

    state.genes = parseGenes(genesText);
    state.peaks = parsePeaks(peaksText);
    state.signals = signalTexts
      .map((text, i) => parseSignal(text, signalFiles[i].name))
      .filter(Boolean);

    if (state.genes.length === 0 && state.peaks.length === 0) {
      status.textContent = 'No features found in the provided files.';
      return;
    }

    const { min, max } = computeDataExtent(
      state.genes,
      state.peaks,
      state.signals,
    );
    // Pad by 1% each side so features never sit flush against the edges.
    const pad = Math.max(1, Math.round((max - min) * 0.01));
    state.dataMin = Math.max(0, min - pad);
    state.dataMax = max + pad;
    state.zoom = clampZoom(ZOOM_STEP ** 8);

    render();
    container.scrollLeft = 0;
    buildSignalControls();
    const signalMsg = state.signals.length
      ? `, ${state.signals.length} signal track${state.signals.length > 1 ? 's' : ''}`
      : '';
    status.textContent =
      `Loaded ${state.genes.length} genes and ${state.peaks.length} peaks${signalMsg}.`;
  } catch (err) {
    console.error(err);
    status.textContent = `Failed to load files: ${err.message}`;
  }
}

// Handles both UTF-8 and UTF-16 (LE/BE) GFF files. The bundled dummy peak.gff
// is UTF-16 LE with a BOM — plain File.text() would still decode it, but being
// explicit avoids surprises for exports from Excel/Notepad.
async function readFileSmart(file) {
  const buf = await file.arrayBuffer();
  const bytes = new Uint8Array(buf);
  let encoding = 'utf-8';
  if (bytes.length >= 2 && bytes[0] === 0xff && bytes[1] === 0xfe) {
    encoding = 'utf-16le';
  } else if (bytes.length >= 2 && bytes[0] === 0xfe && bytes[1] === 0xff) {
    encoding = 'utf-16be';
  }
  return new TextDecoder(encoding).decode(buf);
}

function parseGenes(text) {
  return text
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter((line) => line && !line.startsWith('#'))
    .map((line) => {
      const cols = line.split('\t');
      const start = Number(cols[3]);
      const end = Number(cols[4]);
      const direction = cols[6] === '-' ? '-' : '+';
      const name = extractName(cols[8] || '');
      return { start, end, direction, name };
    })
    .filter((row) => Number.isFinite(row.start) && Number.isFinite(row.end));
}

function parsePeaks(text) {
  return text
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter((line) => line && !line.startsWith('#'))
    .map((line) => {
      const cols = line.split('\t');
      const start = Number(cols[3]);
      const end = Number(cols[4]);
      const score = Number(cols[5]);
      return { start, end, score: Number.isFinite(score) ? score : null };
    })
    .filter((row) => Number.isFinite(row.start) && Number.isFinite(row.end));
}

// Signal GFF: each row is a single position (start == end) with an integer
// score. Files often list all positive scores first and then repeat the same
// coordinate range with negative scores, so we split them into two dense
// Int32Arrays indexed by (pos - dataMin). Both are rendered as mirrored
// halves of a single centered histogram lane.
function parseSignal(text, name) {
  const lines = text.split(/\r?\n/);
  const positions = [];
  const scores = [];
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    if (!line) continue;
    if (line.charCodeAt(0) === 35 /* # */) continue;
    const cols = line.split('\t');
    if (cols.length < 6) continue;
    const p = +cols[3];
    const s = +cols[5];
    if (p === p && s === s) {   // NaN-safe: NaN !== NaN
      positions.push(p);
      scores.push(s);
    }
  }
  if (positions.length === 0) return null;

  let min = positions[0];
  let max = positions[0];
  for (let i = 1; i < positions.length; i++) {
    const p = positions[i];
    if (p < min) min = p;
    else if (p > max) max = p;
  }
  const span = max - min + 1;
  const posData = new Int32Array(span);
  let negData = null;   // allocated lazily — most tracks have no negatives
  let posMax = 0;
  let negMax = 0;       // magnitude of the most-negative score

  for (let i = 0; i < positions.length; i++) {
    const idx = positions[i] - min;
    const s = scores[i] | 0;
    if (s >= 0) {
      if (s > posData[idx]) posData[idx] = s;
      if (s > posMax) posMax = s;
    } else {
      if (negData === null) negData = new Int32Array(span);
      const a = -s;
      if (a > negData[idx]) negData[idx] = a;
      if (a > negMax) negMax = a;
    }
  }
  return {
    name,
    dataMin: min,
    dataMax: max,
    posData,
    negData,
    posMax,
    negMax,
    // User-tunable display caps (null → auto). Set via the controls panel.
    viewPosMax: null,
    viewNegMax: null,
  };
}

function extractName(attr) {
  const match = attr.match(/name=([^;"]+)/i);
  return match ? match[1].trim() : 'unknown';
}

function computeDataExtent(genes, peaks, signals = []) {
  let min = Infinity;
  let max = -Infinity;
  for (const g of genes) {
    if (g.start < min) min = g.start;
    if (g.end > max) max = g.end;
  }
  for (const p of peaks) {
    if (p.start < min) min = p.start;
    if (p.end > max) max = p.end;
  }
  for (const s of signals) {
    if (s.dataMin < min) min = s.dataMin;
    if (s.dataMax > max) max = s.dataMax;
  }
  if (!Number.isFinite(min) || !Number.isFinite(max) || min === max) {
    return { min: 0, max: 1 };
  }
  return { min, max };
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

const SVG_NS = 'http://www.w3.org/2000/svg';

function svgEl(tag, attrs = {}) {
  const el = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) {
    el.setAttribute(k, String(v));
  }
  return el;
}

function getNiceTickStep(span, targetTicks = 10) {
  if (span <= 0) return 1;
  const roughStep = span / targetTicks;
  const magnitude = 10 ** Math.floor(Math.log10(roughStep));
  const normalized = roughStep / magnitude;
  let nice;
  if (normalized < 1.5) nice = 1;
  else if (normalized < 3) nice = 2;
  else if (normalized < 7) nice = 5;
  else nice = 10;
  return nice * magnitude;
}

function getTickValues(start, end, targetTicks) {
  const step = getNiceTickStep(end - start, targetTicks);
  const first = Math.ceil(start / step) * step;
  const ticks = [];
  for (let v = first; v <= end; v += step) {
    ticks.push(Math.round(v));
  }
  return ticks;
}

function totalWidthPx() {
  const viewport = Math.max(600, container.clientWidth || 1200);
  return Math.min(MAX_TOTAL_WIDTH, Math.round(viewport * state.zoom));
}

function render() {
  const dataSpan = Math.max(1, state.dataMax - state.dataMin);
  const width = totalWidthPx();
  const margin = { top: 24, right: 24, bottom: 46, left: 24 };

  const laneHeight = 22;
  const laneGap = 10;
  const peakHeight = 18;
  const basePlotHeight = 270; // preserves original no-signal layout exactly

  // Extra vertical space needed to fit the histogram lanes.
  const signalsBlock = state.signals.length
    * (SIGNAL_LANE_HEIGHT + SIGNAL_LANE_GAP);
  const plotHeight = basePlotHeight + signalsBlock;
  const height = plotHeight + margin.top + margin.bottom;
  const plotWidth = width - margin.left - margin.right;

  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.setAttribute('width', width);
  svg.setAttribute('height', height);
  svg.style.width = `${width}px`;
  svg.style.height = `${height}px`;
  svg.replaceChildren();

  const xOf = (coord) =>
    margin.left + ((coord - state.dataMin) / dataSpan) * plotWidth;

  // Gene lane geometry (unchanged from the original layout).
  const axisY = margin.top + basePlotHeight / 2;
  const forwardY = axisY - laneGap - laneHeight;
  const reverseY = axisY + laneGap;

  // Signal lanes stack directly below the reverse-strand gene lane. When
  // there are no signals, peakY falls back to its original position.
  const firstSignalY = reverseY + laneHeight + SIGNAL_LANE_GAP;
  const signalYs = state.signals.map(
    (_, i) => firstSignalY + i * (SIGNAL_LANE_HEIGHT + SIGNAL_LANE_GAP),
  );
  const peakY = state.signals.length
    ? signalYs[signalYs.length - 1] + SIGNAL_LANE_HEIGHT + SIGNAL_LANE_GAP
    : axisY + laneGap + laneHeight + 14;

  // ---- Grid & ticks ----
  const targetTicks = Math.max(6, Math.round(plotWidth / 110));
  const ticks = getTickValues(state.dataMin, state.dataMax, targetTicks);
  for (const t of ticks) {
    const xPos = xOf(t);
    svg.appendChild(svgEl('line', {
      x1: xPos, x2: xPos,
      y1: margin.top, y2: margin.top + plotHeight,
      class: 'grid-line',
    }));
    const label = svgEl('text', {
      x: xPos,
      y: margin.top + plotHeight + 20,
      'text-anchor': 'middle',
      class: 'tick-label',
    });
    label.textContent = t.toLocaleString();
    svg.appendChild(label);
  }

  // ---- Central axis ----
  svg.appendChild(svgEl('line', {
    x1: margin.left, x2: margin.left + plotWidth,
    y1: axisY, y2: axisY,
    class: 'axis-line',
  }));

  // ---- Genes ----
  for (const gene of state.genes) {
    const x1 = xOf(gene.start);
    const x2 = xOf(gene.end);
    const w = Math.max(2, x2 - x1);
    const y = gene.direction === '+' ? forwardY : reverseY;

    const group = svgEl('g', { class: 'feature gene' });
    group.appendChild(svgEl('rect', {
      x: Math.min(x1, x2), y,
      width: w, height: laneHeight,
      rx: 4, class: 'gene-bar',
    }));

    if (w >= LABEL_MIN_WIDTH_PX) {
      const label = svgEl('text', {
        x: Math.min(x1, x2) + 6,
        y: y + laneHeight / 2 + 4,
        class: 'gene-label',
      });
      label.textContent = gene.name;
      group.appendChild(label);
    }

    const title = svgEl('title');
    title.textContent =
      `${gene.name}\nStrand: ${gene.direction}\n` +
      `${gene.start.toLocaleString()}\u2013${gene.end.toLocaleString()} bp ` +
      `(${(gene.end - gene.start + 1).toLocaleString()} bp)`;
    group.appendChild(title);

    svg.appendChild(group);
  }

  // ---- Signal histograms ----
  layoutState.margin = margin;
  layoutState.plotWidth = plotWidth;
  layoutState.signalYs = signalYs;
  for (let i = 0; i < state.signals.length; i++) {
    renderSignalLane(state.signals[i], signalYs[i], plotWidth, xOf, margin);
  }

  // ---- Peaks ----
  for (const peak of state.peaks) {
    const x1 = xOf(peak.start);
    const x2 = xOf(peak.end);
    const w = Math.max(2, x2 - x1);

    const group = svgEl('g', { class: 'feature peak' });
    group.appendChild(svgEl('rect', {
      x: Math.min(x1, x2), y: peakY,
      width: w, height: peakHeight,
      rx: 3, class: 'peak-bar',
    }));

    if (w >= LABEL_MIN_WIDTH_PX && peak.score !== null) {
      const label = svgEl('text', {
        x: Math.min(x1, x2) + 6,
        y: peakY + peakHeight / 2 + 4,
        class: 'peak-label',
      });
      label.textContent = peak.score.toFixed(3);
      group.appendChild(label);
    }

    const title = svgEl('title');
    title.textContent =
      `Peak\n${peak.start.toLocaleString()}\u2013${peak.end.toLocaleString()} bp` +
      (peak.score !== null ? `\nScore: ${peak.score}` : '');
    group.appendChild(title);

    svg.appendChild(group);
  }
}

// Signal lane: centered horizontal axis, positive scores mirrored above,
// negative-score magnitudes mirrored below. Rendering is viewport-scoped —
// we only walk bins that overlap the visible region (plus some padding), at
// ~1 bin per CSS pixel. This makes both very fine detail at high zoom AND
// panning cheap. Custom viewPosMax/viewNegMax clamp the display range so
// tall outliers don't flatten the smaller peaks.
function renderSignalLane(signal, laneTop, plotWidth, xOf, margin) {
  const {
    name, dataMin: sMin, dataMax: sMax,
    posData, negData, posMax, negMax,
    viewPosMax, viewNegMax,
  } = signal;

  const centerY = laneTop + SIGNAL_LANE_HEIGHT / 2;
  const halfH = SIGNAL_LANE_HEIGHT / 2;

  const group = svgEl('g', { class: 'signal-lane' });

  // Full-range background so empty stretches read as "no data" rather than
  // as the parent chart backdrop.
  group.appendChild(svgEl('rect', {
    x: xOf(state.dataMin), y: laneTop,
    width: Math.max(0, xOf(state.dataMax) - xOf(state.dataMin)),
    height: SIGNAL_LANE_HEIGHT,
    class: 'signal-lane-bg',
  }));
  group.appendChild(svgEl('line', {
    x1: xOf(state.dataMin), x2: xOf(state.dataMax),
    y1: centerY, y2: centerY,
    class: 'signal-axis',
  }));

  const effPosMax = viewPosMax != null ? viewPosMax : posMax;
  const effNegMax = viewNegMax != null ? viewNegMax : negMax;

  // Clip work to the visible x range in data coordinates.
  const dataSpan = Math.max(1, state.dataMax - state.dataMin);
  const viewWidth = Math.max(1, container.clientWidth);
  const padPx = viewWidth * SIGNAL_VIEWPORT_PAD_FRAC;
  const visStartPx = Math.max(margin.left, container.scrollLeft - padPx);
  const visEndPx = Math.min(
    margin.left + plotWidth,
    container.scrollLeft + viewWidth + padPx,
  );
  const visDataMin =
    state.dataMin + ((visStartPx - margin.left) / plotWidth) * dataSpan;
  const visDataMax =
    state.dataMin + ((visEndPx - margin.left) / plotWidth) * dataSpan;

  const rangeMin = Math.max(sMin, Math.floor(visDataMin));
  const rangeMax = Math.min(sMax, Math.ceil(visDataMax));

  if (rangeMin <= rangeMax) {
    const visSpanBp = rangeMax - rangeMin + 1;
    const visSpanPx = Math.max(1, xOf(rangeMax) - xOf(rangeMin));
    const targetBins = Math.min(
      SIGNAL_MAX_BINS_PER_LANE,
      Math.max(50, Math.round(visSpanPx * SIGNAL_BINS_PER_PX)),
    );
    const binBp = Math.max(1, Math.ceil(visSpanBp / targetBins));

    if (effPosMax > 0) {
      appendStepPath(
        group, posData, sMin, rangeMin, rangeMax, binBp,
        xOf, centerY, -halfH / effPosMax, effPosMax,
        'signal-path signal-path-pos',
      );
    }
    if (negData && effNegMax > 0) {
      appendStepPath(
        group, negData, sMin, rangeMin, rangeMax, binBp,
        xOf, centerY, halfH / effNegMax, effNegMax,
        'signal-path signal-path-neg',
      );
    }
  }

  // Label follows the viewport left edge so the track name stays visible
  // while panning.
  const labelX = Math.min(
    xOf(state.dataMax) - 4,
    Math.max(xOf(state.dataMin) + 6, container.scrollLeft + margin.left + 6),
  );
  const label = svgEl('text', {
    x: labelX,
    y: laneTop + 12,
    class: 'signal-label',
  });
  const rangeText = negData
    ? `+${posMax} / -${negMax}`
    : `max ${posMax}`;
  label.textContent = `${name}  (${rangeText})`;
  group.appendChild(label);

  const title = svgEl('title');
  title.textContent =
    `${name}\n${sMin.toLocaleString()}\u2013${sMax.toLocaleString()} bp\n` +
    `positive max ${posMax}` +
    (negData ? `\nnegative max ${negMax}` : '');
  group.appendChild(title);

  svg.appendChild(group);
}

// Draws a step-histogram path for one half of a signal lane.
// `yPerUnit` maps a data value → pixel offset from `baselineY`. Use a
// negative value to draw upward from the baseline, positive to draw down.
// Values above `cap` are clipped so the caller controls the visual ceiling.
function appendStepPath(
  group, data, sMin, rangeMin, rangeMax, binBp,
  xOf, baselineY, yPerUnit, cap, className,
) {
  const startIdx = Math.max(0, rangeMin - sMin);
  const endIdx = Math.min(data.length - 1, rangeMax - sMin);
  if (startIdx > endIdx) return;

  // Snap bin boundaries to multiples of binBp so pos/neg halves align.
  const startBin = Math.floor(startIdx / binBp);
  const endBin = Math.floor(endIdx / binBp);

  const parts = [];
  let started = false;
  let lastX2 = 0;

  for (let bin = startBin; bin <= endBin; bin++) {
    const lo = bin * binBp;
    const hi = Math.min(lo + binBp - 1, data.length - 1);
    if (lo >= data.length) break;

    let m = 0;
    for (let i = lo; i <= hi; i++) {
      const v = data[i];
      if (v > m) m = v;
    }
    const clipped = m > cap ? cap : m;

    const x1 = xOf(sMin + lo);
    const x2 = xOf(sMin + hi + 1);
    const y = baselineY + clipped * yPerUnit;

    if (!started) {
      parts.push(`M${x1.toFixed(2)} ${baselineY.toFixed(2)}`);
      parts.push(`L${x1.toFixed(2)} ${y.toFixed(2)}`);
      started = true;
    } else {
      parts.push(`L${x1.toFixed(2)} ${y.toFixed(2)}`);
    }
    parts.push(`L${x2.toFixed(2)} ${y.toFixed(2)}`);
    lastX2 = x2;
  }

  if (started) {
    parts.push(`L${lastX2.toFixed(2)} ${baselineY.toFixed(2)} Z`);
    group.appendChild(svgEl('path', { d: parts.join(' '), class: className }));
  }
}

// Re-renders just the signal lanes using cached layout. Called from the
// scroll handler (rAF-throttled) and after the user tweaks Y-axis caps.
function renderSignalsOnly() {
  if (!layoutState.margin || state.signals.length === 0) return;
  svg.querySelectorAll('g.signal-lane').forEach((n) => n.remove());
  const { margin, plotWidth, signalYs } = layoutState;
  const dataSpan = Math.max(1, state.dataMax - state.dataMin);
  const xOf = (coord) =>
    margin.left + ((coord - state.dataMin) / dataSpan) * plotWidth;
  for (let i = 0; i < state.signals.length; i++) {
    renderSignalLane(state.signals[i], signalYs[i], plotWidth, xOf, margin);
  }
}

// Builds one row of Y-axis cap inputs per signal underneath the toolbar.
function buildSignalControls() {
  const host = document.getElementById('signal-controls');
  host.replaceChildren();
  if (state.signals.length === 0) {
    host.hidden = true;
    return;
  }
  host.hidden = false;

  const heading = document.createElement('div');
  heading.className = 'signal-controls-heading';
  heading.textContent = 'Signal Y-axis caps (blank = auto)';
  host.appendChild(heading);

  for (const signal of state.signals) {
    const row = document.createElement('div');
    row.className = 'signal-control-row';

    const nameEl = document.createElement('span');
    nameEl.className = 'signal-control-name';
    nameEl.textContent = signal.name;
    row.appendChild(nameEl);

    row.appendChild(makeCapInput(
      signal, 'viewPosMax', '+max', `auto ${signal.posMax}`,
    ));
    if (signal.negData) {
      row.appendChild(makeCapInput(
        signal, 'viewNegMax', '−max', `auto ${signal.negMax}`,
      ));
    }

    host.appendChild(row);
  }
}

function makeCapInput(signal, field, labelText, placeholder) {
  const wrap = document.createElement('label');
  wrap.className = 'signal-cap';
  const lbl = document.createElement('span');
  lbl.textContent = labelText;
  const input = document.createElement('input');
  input.type = 'number';
  input.min = '1';
  input.step = '1';
  input.placeholder = placeholder;
  input.addEventListener('input', () => {
    const raw = input.value.trim();
    const v = raw === '' ? null : Math.max(1, Number(raw));
    signal[field] = Number.isFinite(v) ? v : null;
    renderSignalsOnly();
  });
  wrap.appendChild(lbl);
  wrap.appendChild(input);
  return wrap;
}

// ---------------------------------------------------------------------------
// Zoom & pan
// ---------------------------------------------------------------------------

function clampZoom(z) {
  const viewport = Math.max(600, container.clientWidth || 1200);
  const maxZoom = MAX_TOTAL_WIDTH / viewport;
  return Math.max(MIN_ZOOM, Math.min(maxZoom, z));
}

function setZoom(nextZoom, anchorClientX) {
  if (state.dataMax <= state.dataMin) return;

  const rect = container.getBoundingClientRect();
  const anchor =
    anchorClientX === undefined
      ? container.clientWidth / 2
      : anchorClientX - rect.left;

  const oldWidth = totalWidthPx();
  const coordAtAnchor =
    state.dataMin +
    ((container.scrollLeft + anchor) / oldWidth) *
      (state.dataMax - state.dataMin);

  state.zoom = clampZoom(nextZoom);
  render();

  const newWidth = totalWidthPx();
  const newScrollLeft =
    ((coordAtAnchor - state.dataMin) / (state.dataMax - state.dataMin)) *
      newWidth -
    anchor;

  const maxScroll = Math.max(0, newWidth - container.clientWidth);
  container.scrollLeft = Math.max(0, Math.min(maxScroll, newScrollLeft));
}

function zoomBy(factor, anchorClientX) {
  setZoom(state.zoom * factor, anchorClientX);
}

function panByPixels(dx) {
  container.scrollBy({ left: dx, behavior: 'smooth' });
}

function resetView() {
  state.zoom = clampZoom(ZOOM_STEP ** 8);
  render();
  container.scrollLeft = 0;
}

// ---------------------------------------------------------------------------
// Events
// ---------------------------------------------------------------------------

document.getElementById('zoom-in').addEventListener('click', () => zoomBy(ZOOM_STEP));
document.getElementById('zoom-out').addEventListener('click', () => zoomBy(1 / ZOOM_STEP));
document.getElementById('reset').addEventListener('click', resetView);
document.getElementById('pan-left').addEventListener('click', () =>
  panByPixels(-container.clientWidth * 0.6),
);
document.getElementById('pan-right').addEventListener('click', () =>
  panByPixels(container.clientWidth * 0.6),
);
document.getElementById('load-files').addEventListener('click', loadDataFromFiles);

// Ctrl/Cmd + wheel → zoom toward cursor. Plain wheel → horizontal pan.
container.addEventListener(
  'wheel',
  (event) => {
    if (event.ctrlKey || event.metaKey) {
      event.preventDefault();
      const factor = event.deltaY < 0 ? ZOOM_STEP : 1 / ZOOM_STEP;
      zoomBy(factor, event.clientX);
    } else {
      const delta = event.deltaX !== 0 ? event.deltaX : event.deltaY;
      if (delta !== 0) {
        event.preventDefault();
        container.scrollLeft += delta;
      }
    }
  },
  { passive: false },
);

// Signal lanes render only the visible slice, so we redraw them on pan.
// rAF-throttled so a burst of scroll events collapses to at most one repaint.
let signalScrollRafId = 0;
container.addEventListener('scroll', () => {
  if (state.signals.length === 0) return;
  if (signalScrollRafId) return;
  signalScrollRafId = requestAnimationFrame(() => {
    signalScrollRafId = 0;
    renderSignalsOnly();
  });
});

// Keyboard: ←/→ pan, +/- zoom, 0 reset. Only fires when the chart is focused.
container.addEventListener('keydown', (event) => {
  switch (event.key) {
    case 'ArrowLeft':
      event.preventDefault();
      panByPixels(-container.clientWidth * 0.3);
      break;
    case 'ArrowRight':
      event.preventDefault();
      panByPixels(container.clientWidth * 0.3);
      break;
    case '+':
    case '=':
      event.preventDefault();
      zoomBy(ZOOM_STEP);
      break;
    case '-':
    case '_':
      event.preventDefault();
      zoomBy(1 / ZOOM_STEP);
      break;
    case '0':
      event.preventDefault();
      resetView();
      break;
  }
});

// Debounced re-render on resize so the SVG keeps filling the viewport.
let resizeTimer = null;
window.addEventListener('resize', () => {
  if (
    state.genes.length === 0 &&
    state.peaks.length === 0 &&
    state.signals.length === 0
  ) return;
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    const anchor =
      container.getBoundingClientRect().left + container.clientWidth / 2;
    setZoom(state.zoom, anchor);
  }, 120);
});

// ---------------------------------------------------------------------------
// Drag-to-zoom selection
// ---------------------------------------------------------------------------

const selectionOverlay = document.getElementById('selection-overlay');
let dragState = null;

container.addEventListener('mousedown', (event) => {
  // Only respond to primary button, ignore if Ctrl is held (that's for wheel zoom)
  if (event.button !== 0 || event.ctrlKey || event.metaKey) return;
  // Don't start drag on scrollbar area
  if (event.offsetY > container.clientHeight - 12) return;

  dragState = {
    startX: event.clientX,
    scrollLeftAtStart: container.scrollLeft,
  };

  // Prevent text selection while dragging
  event.preventDefault();
});

window.addEventListener('mousemove', (event) => {
  if (!dragState) return;

  const currentX = event.clientX;
  const rect = container.getBoundingClientRect();

  // Compute positions relative to the container's left edge
  const startRel = dragState.startX - rect.left;
  const currentRel = currentX - rect.left;

  const left = Math.max(0, Math.min(startRel, currentRel));
  const right = Math.min(container.clientWidth, Math.max(startRel, currentRel));
  const width = right - left;

  if (width > 4) {
    // Position overlay accounting for scroll — it's inside the scrollable container
    selectionOverlay.style.left = `${left + container.scrollLeft}px`;
    selectionOverlay.style.width = `${width}px`;
    selectionOverlay.hidden = false;
  } else {
    selectionOverlay.hidden = true;
  }
});

window.addEventListener('mouseup', (event) => {
  if (!dragState) return;

  selectionOverlay.hidden = true;

  const rect = container.getBoundingClientRect();
  const startRel = dragState.startX - rect.left;
  const endRel = event.clientX - rect.left;

  const leftPx = Math.max(0, Math.min(startRel, endRel));
  const rightPx = Math.min(container.clientWidth, Math.max(startRel, endRel));
  const widthPx = rightPx - leftPx;

  dragState = null;

  // Only zoom if the user dragged a meaningful distance (> 8px)
  if (widthPx < 8 || state.dataMax <= state.dataMin) return;

  // Convert pixel positions to data coordinates
  const svgWidth = totalWidthPx();
  const dataSpan = state.dataMax - state.dataMin;

  const coordLeft =
    state.dataMin +
    ((container.scrollLeft + leftPx) / svgWidth) * dataSpan;
  const coordRight =
    state.dataMin +
    ((container.scrollLeft + rightPx) / svgWidth) * dataSpan;

  // Compute new zoom so that [coordLeft, coordRight] fills the viewport
  const selectedFraction = (coordRight - coordLeft) / dataSpan;
  const newZoom = clampZoom(1 / selectedFraction);

  state.zoom = newZoom;
  render();

  // Scroll so coordLeft aligns with the left edge of the viewport
  const newSvgWidth = totalWidthPx();
  const newScrollLeft =
    ((coordLeft - state.dataMin) / dataSpan) * newSvgWidth;

  container.scrollLeft = Math.max(
    0,
    Math.min(newSvgWidth - container.clientWidth, newScrollLeft),
  );
});