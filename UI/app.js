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
const signalCanvas = document.getElementById('signal-canvas');
const signalCtx = signalCanvas.getContext('2d');
const rulerOverlay = document.getElementById('ruler-overlay');
const rulerLabel = document.getElementById('ruler-label');

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
  colors: {
    gene: '#2e9f57',
    peak: '#3b82f6',
  },
  dataMin: 0,
  dataMax: 1,
  zoom: 1,
  // Order of lanes from top to bottom. Each entry: { kind, index?, key }.
  // kind ∈ { 'gene+', 'gene-', 'signal', 'peak' }. `key` is unique per lane.
  // Reset to null on data load so defaultLaneOrder() rebuilds it.
  laneOrder: null,
};

// Per-lane vertical size in SVG pixels. Signal lanes are much taller since
// they contain a histogram, not just a bar. LANE_GAP is the vertical space
// inserted between adjacent lanes.
const LANE_HEIGHTS = {
  'gene+': 22,
  'gene-': 22,
  signal: SIGNAL_LANE_HEIGHT,
  peak: 18,
};
const LANE_GAP = 10;
// Gene strand lanes move as one group, keeping + above −.
const GENE_GROUP_KEY = 'gene-group';

// Saved geometry from the last render() so we can re-render only the signal
// lanes on scroll without recomputing the full layout.
const layoutState = {
  margin: null,
  plotWidth: 0,
  plotHeight: 0,
  laneRects: null, // Map<laneKey, { top, height, kind, index, key }>
  handleEls: null, // Map<laneKey, HTMLDivElement[]>
};

// ---------------------------------------------------------------------------
// Data loading
// ---------------------------------------------------------------------------

async function loadDataFromFiles() {
  const genesFile = document.getElementById('genes-file').files[0];
  const peaksFiles = Array.from(
    document.getElementById('peaks-file').files || [],
  );
  const signalFiles = Array.from(
    document.getElementById('signals-file').files || [],
  );

  if (!genesFile || peaksFiles.length === 0) {
    status.textContent = 'Please select both files first.';
    return;
  }

  try {
    const [genesText, ...fileTexts] = await Promise.all([
      readFileSmart(genesFile),
      ...peaksFiles.map(readFileSmart),
      ...signalFiles.map(readFileSmart),
    ]);
    const peakTexts = fileTexts.slice(0, peaksFiles.length);
    const signalTexts = fileTexts.slice(peaksFiles.length);

    state.genes = parseGenes(genesText);
    state.peaks = peakTexts.map((text, i) => ({
      name: peaksFiles[i].name,
      peaks: parsePeaks(text),
    }));
    state.signals = signalTexts
      .map((text, i) => parseSignal(text, signalFiles[i].name))
      .filter(Boolean);
    // Data changed — reset lane order so it's regenerated with the new lanes.
    state.laneOrder = null;

    const peakCount = state.peaks.reduce((total, lane) => total + lane.peaks.length, 0);
    if (state.genes.length === 0 && peakCount === 0) {
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
      `Loaded ${state.genes.length} genes and ${peakCount} peaks across ` +
      `${state.peaks.length} peak lane${state.peaks.length > 1 ? 's' : ''}${signalMsg}.`;
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
    color: '#f49e4c',
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
  for (const lane of peaks) {
    for (const p of lane.peaks) {
      if (p.start < min) min = p.start;
      if (p.end > max) max = p.end;
    }
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

function hexToRgba(hex, alpha) {
  const value = hex.replace('#', '');
  const red = parseInt(value.slice(0, 2), 16);
  const green = parseInt(value.slice(2, 4), 16);
  const blue = parseInt(value.slice(4, 6), 16);
  return `rgba(${red}, ${green}, ${blue}, ${alpha})`;
}

function contrastingTextColor(hex) {
  const value = hex.replace('#', '');
  const channels = [0, 2, 4].map((offset) => parseInt(value.slice(offset, offset + 2), 16) / 255);
  const luminance = channels.reduce((total, channel, index) => {
    const linear = channel <= 0.04045
      ? channel / 12.92
      : ((channel + 0.055) / 1.055) ** 2.4;
    return total + linear * [0.2126, 0.7152, 0.0722][index];
  }, 0);
  const contrastWithBlack = (luminance + 0.05) / 0.05;
  const contrastWithWhite = 1.05 / (luminance + 0.05);
  return contrastWithBlack >= contrastWithWhite ? '#000000' : '#ffffff';
}

function applyColors() {
  document.documentElement.style.setProperty('--gene', state.colors.gene);
  document.documentElement.style.setProperty('--gene-hover', state.colors.gene);
  document.documentElement.style.setProperty('--gene-label', contrastingTextColor(state.colors.gene));
  document.documentElement.style.setProperty('--peak', state.colors.peak);
  document.documentElement.style.setProperty('--peak-hover', state.colors.peak);
  document.documentElement.style.setProperty('--peak-label', contrastingTextColor(state.colors.peak));
  renderSignalsOnly();
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

function formatBp(bp) {
  if (bp >= 1_000_000) return `${(bp / 1_000_000).toFixed(2)} Mbp`;
  if (bp >= 1_000) return `${(bp / 1_000).toFixed(2)} kbp`;
  return `${Math.round(bp)} bp`;
}

function clearRuler() {
  rulerOverlay.hidden = true;
  rulerLabel.hidden = true;
  rulerState.active = false;
}

function render() {
  clearRuler();
  syncLaneOrder();

  const dataSpan = Math.max(1, state.dataMax - state.dataMin);
  const width = totalWidthPx();
  const margin = { top: 24, right: 24, bottom: 46, left: 24 };

  // Compute a Y rect per lane based on state.laneOrder. This is the single
  // source of truth for every renderer (SVG genes/peaks, canvas signals,
  // and lane-handle overlays).
  const { laneRects, blockHeight } = computeLaneLayout(margin.top);
  // Preserve a minimum plot height (matches original layout) so the empty
  // state and datasets without signals still look reasonable.
  const plotHeight = Math.max(120, blockHeight);
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

  // ---- Genes ----
  // Draw each gene into whichever strand lane exists in the current order.
  // Lanes may be absent (e.g. peaks-only view) — skip genes with no lane.
  const genePosRect = laneRects.get('gene+');
  const geneNegRect = laneRects.get('gene-');
  for (const gene of state.genes) {
    const rect = gene.direction === '+' ? genePosRect : geneNegRect;
    if (!rect) continue;
    const x1 = xOf(gene.start);
    const x2 = xOf(gene.end);
    const w = Math.max(2, x2 - x1);
    const y = rect.top;

    const group = svgEl('g', { class: 'feature gene' });
    group.appendChild(svgEl('rect', {
      x: Math.min(x1, x2), y,
      width: w, height: rect.height,
      rx: 4, class: 'gene-bar',
    }));

    if (w >= LABEL_MIN_WIDTH_PX) {
      const label = svgEl('text', {
        x: Math.min(x1, x2) + 6,
        y: y + rect.height / 2 + 4,
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

  // ---- Signal histograms (canvas-based) ----
  layoutState.margin = margin;
  layoutState.plotWidth = plotWidth;
  layoutState.plotHeight = plotHeight;
  layoutState.laneRects = laneRects;
  setupSignalCanvas();
  paintSignalCanvas();

  // ---- Peaks ----
  for (const [, rect] of laneRects) {
    if (rect.kind !== 'peak') continue;
    const lane = state.peaks[rect.index];
    if (!lane) continue;
    const peakY = rect.top;

    for (const peak of lane.peaks) {
      const x1 = xOf(peak.start);
      const x2 = xOf(peak.end);
      const w = Math.max(2, x2 - x1);

      const group = svgEl('g', { class: 'feature peak' });
      group.appendChild(svgEl('rect', {
        x: Math.min(x1, x2), y: peakY,
        width: w, height: rect.height,
        rx: 3, class: 'peak-bar',
      }));

      if (w >= LABEL_MIN_WIDTH_PX && peak.score !== null) {
        const label = svgEl('text', {
          x: Math.min(x1, x2) + 6,
          y: peakY + rect.height / 2 + 4,
          class: 'peak-label',
        });
        label.textContent = peak.score.toFixed(3);
        group.appendChild(label);
      }

      const title = svgEl('title');
      title.textContent =
        `${lane.name}\nPeak\n${peak.start.toLocaleString()}\u2013${peak.end.toLocaleString()} bp` +
        (peak.score !== null ? `\nScore: ${peak.score}` : '');
      group.appendChild(title);

      svg.appendChild(group);
    }
  }

  // ---- Draggable lane handles overlay ----
  buildLaneHandles(laneRects, margin);
}

// Build the default lane order for the current dataset. Skips gene lanes
// when no genes are loaded so peaks-only views don't have empty strand rows.
function defaultLaneOrder() {
  const order = [];
  if (state.genes.length > 0) {
    order.push({ kind: 'gene+', key: 'gene+' });
    order.push({ kind: 'gene-', key: 'gene-' });
  }
  state.signals.forEach((_, i) => order.push({ kind: 'signal', index: i, key: `signal:${i}` }));
  state.peaks.forEach((_, i) => order.push({ kind: 'peak', index: i, key: `peak:${i}` }));
  return order;
}

// Reconciles state.laneOrder with the currently-loaded data. Keeps any
// user reordering intact across re-renders, drops stale keys, and appends
// any newly-added lanes at the end. Gene strand lanes are always forced
// to the top of the order (+ above −), regardless of prior state.
function syncLaneOrder() {
  const defaults = defaultLaneOrder();
  const wanted = new Set(defaults.map((d) => d.key));
  const existing = (state.laneOrder || []).filter((l) => wanted.has(l.key));
  const have = new Set(existing.map((l) => l.key));
  const missing = defaults.filter((d) => !have.has(d.key));
  const merged = [...existing, ...missing];
    state.laneOrder = merged;
}

// Assigns a { top, height, kind, index } rect to each lane based on its
// position in state.laneOrder. Returns a Map keyed by lane key plus the
// total height of the lane block (for sizing the plot).
function computeLaneLayout(startY) {
  const laneRects = new Map();
  let y = startY;
  let first = true;
  for (const lane of state.laneOrder) {
    const h = LANE_HEIGHTS[lane.kind] ?? 20;
    if (!first) y += LANE_GAP;
    laneRects.set(lane.key, {
      top: y,
      height: h,
      kind: lane.kind,
      index: lane.index,
      key: lane.key,
    });
    y += h;
    first = false;
  }
  return { laneRects, blockHeight: first ? 0 : y - startY };
}

// ---------------------------------------------------------------------------
// Canvas-based signal rendering
// The canvas is viewport-sized and positioned over the signal lanes area.
// On each scroll/zoom we just clear + repaint.
// ---------------------------------------------------------------------------

function setupSignalCanvas() {
  if (state.signals.length === 0 && state.peaks.length === 0) {
    signalCanvas.style.display = 'none';
    return;
  }
  const dpr = window.devicePixelRatio || 1;
  const viewWidth = container.clientWidth;
  const canvasH = signalCanvasHeight();
  if (canvasH <= 0) {
    signalCanvas.style.display = 'none';
    return;
  }
  signalCanvas.style.display = 'block';
  signalCanvas.style.width = `${viewWidth}px`;
  signalCanvas.style.height = `${canvasH}px`;
  // The canvas now covers the entire plot area (from margin.top down to
  // the bottom of the lane block) so signal + peak lanes can appear
  // anywhere in the lane order. We pull it up over the SVG with a
  // negative margin-top and compensate with margin-bottom so its flow
  // height stays zero.
  const svgHeight = svg.clientHeight || svg.getBoundingClientRect().height;
  const overlapPx = svgHeight - layoutState.margin.top;
  signalCanvas.style.marginTop = `-${overlapPx}px`;
  signalCanvas.style.marginBottom = `${overlapPx - canvasH}px`;

  signalCanvas.width = Math.round(viewWidth * dpr);
  signalCanvas.height = Math.round(canvasH * dpr);
  signalCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
}

function paintSignalCanvas() {
  if ((state.signals.length === 0 && state.peaks.length === 0) || !layoutState.margin) return;
  const dpr = window.devicePixelRatio || 1;
  const viewWidth = container.clientWidth;
  const canvasH = signalCanvasHeight();
  if (canvasH <= 0) return;

  // Ensure canvas is sized for current viewport (handles resize).
  const neededW = Math.round(viewWidth * dpr);
  const neededH = Math.round(canvasH * dpr);
  if (signalCanvas.width !== neededW || signalCanvas.height !== neededH) {
    signalCanvas.width = neededW;
    signalCanvas.height = neededH;
    signalCanvas.style.width = `${viewWidth}px`;
    signalCanvas.style.height = `${canvasH}px`;
    signalCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  signalCtx.clearRect(0, 0, viewWidth, canvasH);

  const scrollLeft = container.scrollLeft;
  const { margin, plotWidth, laneRects } = layoutState;
  if (!laneRects) return;
  const dataSpan = Math.max(1, state.dataMax - state.dataMin);

  // xOf maps a data coord to SVG pixel space (absolute). We convert to
  // canvas-local by subtracting scrollLeft.
  const dataToCanvasX = (coord) =>
    margin.left + ((coord - state.dataMin) / dataSpan) * plotWidth - scrollLeft;

  // Signal lanes may now sit anywhere in the lane order — resolve each
  // lane's top position from the lane rects. Canvas-local Y = rect.top
  // minus margin.top (since the canvas top edge is at margin.top).
  for (const [, rect] of laneRects) {
    if (rect.kind !== 'signal') continue;
    const signal = state.signals[rect.index];
    if (!signal) continue;
    paintSignalLane(signal, rect.top - margin.top,
      viewWidth, scrollLeft, dataToCanvasX, margin, plotWidth, dataSpan);
  }

  // Sticky peak-lane labels, positioned at each peak lane's Y.
  signalCtx.font = '700 10px system-ui, sans-serif';
  signalCtx.lineWidth = 3;
  signalCtx.strokeStyle = 'rgba(9, 12, 32, 0.75)';
  signalCtx.fillStyle = 'rgba(230, 242, 255, 0.9)';
  const peakLabelX = Math.max(6, dataToCanvasX(state.dataMin) + 6);
  for (const [, rect] of laneRects) {
    if (rect.kind !== 'peak') continue;
    const lane = state.peaks[rect.index];
    if (!lane) continue;
    const laneTop = rect.top - margin.top;
    signalCtx.strokeText(lane.name, peakLabelX, laneTop + 12);
    signalCtx.fillText(lane.name, peakLabelX, laneTop + 12);
  }
}

function signalCanvasHeight() {
  return layoutState.plotHeight || 0;
}

function paintSignalLane(signal, laneLocalTop, viewWidth, scrollLeft, dataToCanvasX, margin, plotWidth, dataSpan) {
  const {
    name, dataMin: sMin, dataMax: sMax,
    posData, negData, posMax, negMax,
    viewPosMax, viewNegMax,
  } = signal;

  const centerY = laneLocalTop + SIGNAL_LANE_HEIGHT / 2;
  const halfH = SIGNAL_LANE_HEIGHT / 2;

  // Background
  const bgLeft = Math.max(0, dataToCanvasX(state.dataMin));
  const bgRight = Math.min(viewWidth, dataToCanvasX(state.dataMax));
  if (bgRight > bgLeft) {
    signalCtx.fillStyle = 'rgba(255, 255, 255, 0.04)';
    signalCtx.fillRect(bgLeft, laneLocalTop, bgRight - bgLeft, SIGNAL_LANE_HEIGHT);
    signalCtx.strokeStyle = 'rgba(255, 255, 255, 0.08)';
    signalCtx.lineWidth = 1;
    signalCtx.strokeRect(bgLeft, laneLocalTop, bgRight - bgLeft, SIGNAL_LANE_HEIGHT);
  }

  // Center axis
  signalCtx.beginPath();
  signalCtx.moveTo(bgLeft, centerY);
  signalCtx.lineTo(bgRight, centerY);
  signalCtx.strokeStyle = 'rgba(255, 255, 255, 0.35)';
  signalCtx.lineWidth = 1;
  signalCtx.stroke();

  const effPosMax = viewPosMax != null ? viewPosMax : posMax;
  const effNegMax = viewNegMax != null ? viewNegMax : negMax;

  // Determine visible data range from scroll position.
  const visStartPx = scrollLeft;
  const visEndPx = scrollLeft + viewWidth;
  const visDataMin = state.dataMin + ((visStartPx - margin.left) / plotWidth) * dataSpan;
  const visDataMax = state.dataMin + ((visEndPx - margin.left) / plotWidth) * dataSpan;

  const rangeMin = Math.max(sMin, Math.floor(visDataMin));
  const rangeMax = Math.min(sMax, Math.ceil(visDataMax));

  if (rangeMin <= rangeMax) {
    const visSpanBp = rangeMax - rangeMin + 1;
    // One bin per CSS pixel for the visible range
    const targetBins = Math.min(SIGNAL_MAX_BINS_PER_LANE, Math.max(50, viewWidth));
    const binBp = Math.max(1, Math.ceil(visSpanBp / targetBins));

    if (effPosMax > 0) {
      paintStepHist(posData, sMin, rangeMin, rangeMax, binBp,
        dataToCanvasX, centerY, -halfH / effPosMax, effPosMax,
        hexToRgba(signal.color, 0.55), signal.color);
    }
    if (negData && effNegMax > 0) {
      paintStepHist(negData, sMin, rangeMin, rangeMax, binBp,
        dataToCanvasX, centerY, halfH / effNegMax, effNegMax,
        hexToRgba(signal.color, 0.55), signal.color);
    }
  }

  // Lane label (sticky to left edge)
  const rangeText = negData ? `+${posMax} / -${negMax}` : `max ${posMax}`;
  const labelText = `${name}  (${rangeText})`;
  signalCtx.font = '700 11px system-ui, sans-serif';
  signalCtx.lineWidth = 3;
  signalCtx.strokeStyle = 'rgba(9, 12, 32, 0.75)';
  signalCtx.fillStyle = 'rgba(230, 242, 255, 0.9)';
  const labelX = Math.max(bgLeft + 6, 6);
  signalCtx.strokeText(labelText, labelX, laneLocalTop + 12);
  signalCtx.fillText(labelText, labelX, laneLocalTop + 12);
}

function paintStepHist(data, sMin, rangeMin, rangeMax, binBp,
  dataToCanvasX, baselineY, yPerUnit, cap, fillColor, strokeColor) {
  const startIdx = Math.max(0, rangeMin - sMin);
  const endIdx = Math.min(data.length - 1, rangeMax - sMin);
  if (startIdx > endIdx) return;

  const startBin = Math.floor(startIdx / binBp);
  const endBin = Math.floor(endIdx / binBp);

  signalCtx.beginPath();
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

    const x1 = dataToCanvasX(sMin + lo);
    const x2 = dataToCanvasX(sMin + hi + 1);
    const y = baselineY + clipped * yPerUnit;

    if (!started) {
      signalCtx.moveTo(x1, baselineY);
      signalCtx.lineTo(x1, y);
      started = true;
    } else {
      signalCtx.lineTo(x1, y);
    }
    signalCtx.lineTo(x2, y);
    lastX2 = x2;
  }

  if (started) {
    signalCtx.lineTo(lastX2, baselineY);
    signalCtx.closePath();
    signalCtx.fillStyle = fillColor;
    signalCtx.fill();
    signalCtx.strokeStyle = strokeColor;
    signalCtx.lineWidth = 0.8;
    signalCtx.stroke();
  }
}

// Re-renders just the signal lanes using cached layout. Called from the
// scroll handler (rAF-throttled) and after the user tweaks Y-axis caps.
function renderSignalsOnly() {
  if (!layoutState.margin || state.signals.length === 0) return;
  paintSignalCanvas();
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

  const headingRow = document.createElement('div');
  headingRow.className = 'signal-controls-heading-row';

  const heading = document.createElement('div');
  heading.className = 'signal-controls-heading';
  heading.textContent = 'Signal Y-axis caps (blank = auto)';
  headingRow.appendChild(heading);

  const masterWrap = document.createElement('label');
  masterWrap.className = 'signal-cap signal-cap-master';
  const masterLbl = document.createElement('span');
  masterLbl.textContent = 'Change all max';
  const masterInput = document.createElement('input');
  masterInput.type = 'number';
  masterInput.min = '1';
  masterInput.step = '1';
  masterInput.placeholder = 'all lanes';
  const masterValues = state.signals.map((signal) => signal.viewPosMax);
  if (
    masterValues.length > 0 &&
    masterValues.every((value) => value != null && value === masterValues[0])
  ) {
    masterInput.value = masterValues[0];
  }
  masterInput.addEventListener('input', () => {
    const raw = masterInput.value.trim();
    const v = raw === '' ? null : Math.max(1, Number(raw));
    const cap = Number.isFinite(v) ? v : null;
    for (const signal of state.signals) {
      signal.viewPosMax = cap;
      if (signal.negData) signal.viewNegMax = cap;
    }
    // Update individual inputs to reflect master value
    host.querySelectorAll('.signal-control-row input[type="number"]').forEach(inp => {
      inp.value = cap != null ? cap : '';
    });
    renderSignalsOnly();
  });
  masterWrap.appendChild(masterLbl);
  masterWrap.appendChild(masterInput);
  headingRow.appendChild(masterWrap);

  host.appendChild(headingRow);

  const orderedSignals = state.laneOrder
    .filter((lane) => lane.kind === 'signal')
    .map((lane) => state.signals[lane.index])
    .filter(Boolean);
  for (const signal of orderedSignals) {
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

    row.appendChild(makeSignalColorInput(signal));

    host.appendChild(row);
  }
}

function makeSignalColorInput(signal) {
  const wrap = document.createElement('label');
  wrap.className = 'color-control signal-color-control';
  const input = document.createElement('input');
  input.type = 'color';
  input.value = signal.color;
  input.title = `Change color for ${signal.name}`;
  input.addEventListener('input', (event) => {
    signal.color = event.target.value;
    renderSignalsOnly();
  });
  wrap.appendChild(input);
  return wrap;
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
  if (signal[field] != null) input.value = signal[field];
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

const colorInputs = [
  ['gene-color', 'gene'],
  ['peak-color', 'peak'],
];
for (const [id, colorName] of colorInputs) {
  document.getElementById(id).addEventListener('input', (event) => {
    state.colors[colorName] = event.target.value;
    applyColors();
  });
}

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

// Signal canvas repaints on scroll — rAF-throttled.
let signalScrollRafId = 0;
container.addEventListener('scroll', () => {
  if (state.signals.length === 0) return;
  if (signalScrollRafId) return;
  signalScrollRafId = requestAnimationFrame(() => {
    signalScrollRafId = 0;
    paintSignalCanvas();
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
let rulerState = { active: false };
let rulerDragState = null;

// Suppress browser context menu on the chart so right-click drag works.
container.addEventListener('contextmenu', (event) => {
  event.preventDefault();
});

container.addEventListener('mousedown', (event) => {
  // Right-click: start a ruler drag, clearing any existing ruler first.
  if (event.button === 2) {
    clearRuler();
    if (event.offsetY <= container.clientHeight - 12) {
      rulerDragState = {
        startX: event.clientX,
        scrollLeftAtStart: container.scrollLeft,
      };
    }
    event.preventDefault();
    return;
  }

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
  if (!dragState && !rulerDragState) return;

  const currentX = event.clientX;
  const rect = container.getBoundingClientRect();

  if (dragState) {
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
  }

  if (rulerDragState) {
    const startRel = rulerDragState.startX - rect.left;
    const currentRel = currentX - rect.left;

    const left = Math.max(0, Math.min(startRel, currentRel));
    const right = Math.min(container.clientWidth, Math.max(startRel, currentRel));
    const width = right - left;

    if (width > 4) {
      rulerOverlay.style.left = `${left + container.scrollLeft}px`;
      rulerOverlay.style.width = `${width}px`;
      rulerOverlay.hidden = false;

      // Live bp distance label, centered above the ruler box
      const svgWidth = totalWidthPx();
      const margin = layoutState.margin || { left: 24, right: 24 };
      const plotWidth = svgWidth - margin.left - margin.right;
      const dataSpan = state.dataMax - state.dataMin;
      const coordLeft = state.dataMin + ((container.scrollLeft + left - margin.left) / plotWidth) * dataSpan;
      const coordRight = state.dataMin + ((container.scrollLeft + right - margin.left) / plotWidth) * dataSpan;
      rulerLabel.textContent = formatBp(Math.abs(coordRight - coordLeft));
      rulerLabel.style.left = `${left + width / 2 + container.scrollLeft}px`;
      rulerLabel.hidden = false;
    } else {
      rulerOverlay.hidden = true;
      rulerLabel.hidden = true;
    }
  }
});

window.addEventListener('mouseup', (event) => {
  // Right-click release: finalize ruler if a meaningful drag occurred.
  if (event.button === 2) {
    if (!rulerDragState) return;

    const rect = container.getBoundingClientRect();
    const startRel = rulerDragState.startX - rect.left;
    const endRel = event.clientX - rect.left;

    const leftPx = Math.max(0, Math.min(startRel, endRel));
    const rightPx = Math.min(container.clientWidth, Math.max(startRel, endRel));
    const widthPx = rightPx - leftPx;

    rulerDragState = null;

    if (widthPx < 4) {
      // Too small to be intentional — ruler already cleared in mousedown.
      return;
    }

    // Convert pixel extents to bp coordinates.
    const svgWidth = totalWidthPx();
    const margin = layoutState.margin || { left: 24, right: 24 };
    const plotWidth = svgWidth - margin.left - margin.right;
    const dataSpan = state.dataMax - state.dataMin;
    const coordLeft = state.dataMin + ((container.scrollLeft + leftPx - margin.left) / plotWidth) * dataSpan;
    const coordRight = state.dataMin + ((container.scrollLeft + rightPx - margin.left) / plotWidth) * dataSpan;

    rulerLabel.textContent = formatBp(Math.abs(coordRight - coordLeft));
    rulerLabel.style.left = `${leftPx + widthPx / 2 + container.scrollLeft}px`;
    rulerLabel.hidden = false;
    rulerOverlay.hidden = false;
    rulerState.active = true;
    return;
  }

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

  // Convert pixel positions to data coordinates, accounting for plot margins
  const margin = layoutState.margin || {left: 24, right: 24};
  const plotWidth = layoutState.plotWidth || (totalWidthPx() - margin.left - margin.right);
  const dataSpan = state.dataMax - state.dataMin;

  const coordLeft =
    state.dataMin +
    ((container.scrollLeft + leftPx - margin.left) / plotWidth) * dataSpan;
  const coordRight =
    state.dataMin +
    ((container.scrollLeft + rightPx - margin.left) / plotWidth) * dataSpan;

  // Clamp to data range
  const clampedLeft = Math.max(state.dataMin, Math.min(coordLeft, coordLeft));
  const clampedRight = Math.max(state.dataMin, Math.max(coordLeft, coordRight));

  if (clampedRight <= clampedLeft) return;

  // Compute new zoom so that [coordLeft, coordRight] fills the viewport
  const selectedFraction = (clampedRight - clampedLeft) / dataSpan;
  const newZoom = clampZoom(1 / selectedFraction);

  state.zoom = newZoom;
  render();

  // Scroll so clampedLeft aligns with the left edge of the viewport
  const newSvgWidth = totalWidthPx();
  const newMargin = layoutState.margin || {left: 24, right: 24};
  const newPlotWidth = newSvgWidth - newMargin.left - newMargin.right;
  const newScrollLeft =
    newMargin.left + ((clampedLeft - state.dataMin) / dataSpan) * newPlotWidth;

  container.scrollLeft = Math.max(
    0,
    Math.min(newSvgWidth - container.clientWidth, newScrollLeft),
  );
});

// ---------------------------------------------------------------------------
// Lane reordering (drag handles)
//
// Every lane in the chart can be reordered by grabbing a small handle that
// appears on the left or right edge when the pointer is over that lane. The
// two gene-strand rows use one grouped handle so + and - move together.
// ---------------------------------------------------------------------------

const chartFrame = document.querySelector('.chart-frame');
const handlesContainer = document.getElementById('lane-handles');
const dropIndicator = document.getElementById('lane-drop-indicator');
let laneDragState = null;

function laneLabelFor(rect) {
  if (rect.kind === 'gene-group') return 'genes';
  if (rect.kind === 'gene+') return 'gene + strand';
  if (rect.kind === 'gene-') return 'gene − strand';
  if (rect.kind === 'signal') return state.signals[rect.index]?.name ?? 'signal';
  if (rect.kind === 'peak')   return state.peaks[rect.index]?.name ?? 'peak';
  return 'lane';
}

// The lane rects use SVG coordinates (y=0 at top of SVG). To translate
// them to chart-frame coordinates we add the SVG's top offset within
// chart-frame (i.e. any chart-scroll border / padding).
function svgTopInFrame() {
  const svgRect = svg.getBoundingClientRect();
  const frameRect = chartFrame.getBoundingClientRect();
  return svgRect.top - frameRect.top;
}

function buildLaneHandles(laneRects, margin) {
  handlesContainer.replaceChildren();
  const handleEls = new Map();
  if (!laneRects || laneRects.size === 0) {
    layoutState.handleEls = handleEls;
    return;
  }
  const yOffset = svgTopInFrame();
  const genePositive = laneRects.get('gene+');
  const geneNegative = laneRects.get('gene-');
  for (const [key, rect] of laneRects) {
    if (rect.kind === 'gene-') continue;
    const isGeneGroup = rect.kind === 'gene+' && geneNegative;
    const handleKey = isGeneGroup ? GENE_GROUP_KEY : key;
    const handleTop = isGeneGroup ? genePositive.top : rect.top;
    const handleHeight = isGeneGroup
      ? geneNegative.top + geneNegative.height - genePositive.top
      : rect.height;
    const els = [];
    const sides = ['left', 'right'];
    for (const side of sides) {
      const el = document.createElement('div');
      el.className = `lane-handle ${side}`;
      el.dataset.laneKey = handleKey;
      el.style.top = `${yOffset + handleTop}px`;
      el.style.height = `${handleHeight}px`;
      el.title = `Drag to reorder ${laneLabelFor(isGeneGroup ? { kind: 'gene-group' } : rect)}`;
      el.addEventListener('mousedown', (ev) => startLaneDrag(ev, handleKey));
      handlesContainer.appendChild(el);
      els.push(el);
    }
    handleEls.set(handleKey, els);
  }
  layoutState.handleEls = handleEls;
}

function setHandleVisibility(hoveredKey) {
  if (!layoutState.handleEls) return;
  for (const [key, els] of layoutState.handleEls) {
    const on = key === hoveredKey;
    for (const el of els) el.classList.toggle('visible', on);
  }
}

function hideAllHandles() {
  setHandleVisibility(null);
}

// Determine which lane the pointer is over (in chart-frame Y coords).
function laneKeyAtY(y) {
  if (!layoutState.laneRects) return null;
  const yOffset = svgTopInFrame();
  const genePositive = layoutState.laneRects.get('gene+');
  const geneNegative = layoutState.laneRects.get('gene-');
  if (genePositive && geneNegative) {
    const geneTop = yOffset + genePositive.top;
    const geneBottom = yOffset + geneNegative.top + geneNegative.height;
    if (y >= geneTop && y <= geneBottom) return GENE_GROUP_KEY;
  }
  for (const [key, rect] of layoutState.laneRects) {
    if (rect.kind === 'gene+' || rect.kind === 'gene-') continue;
    const top = yOffset + rect.top;
    if (y >= top && y <= top + rect.height) return key;
  }
  return null;
}

chartFrame.addEventListener('mousemove', (ev) => {
  if (laneDragState) return; // hover suppressed during drag
  const rect = chartFrame.getBoundingClientRect();
  const y = ev.clientY - rect.top;
  setHandleVisibility(laneKeyAtY(y));
});

chartFrame.addEventListener('mouseleave', () => {
  if (laneDragState) return;
  hideAllHandles();
});

function startLaneDrag(ev, laneKey) {
  if (ev.button !== 0) return;
  ev.preventDefault();
  ev.stopPropagation();
  const fromIdx = laneKey === GENE_GROUP_KEY
    ? state.laneOrder.findIndex((l) => l.kind === 'gene+')
    : state.laneOrder.findIndex((l) => l.key === laneKey);
  if (fromIdx === -1) return;

  const handleEls = layoutState.handleEls?.get(laneKey) ?? [];
  for (const el of handleEls) el.classList.add('dragging', 'visible');

  laneDragState = {
    laneKey,
    fromIdx,
    groupSize: laneKey === GENE_GROUP_KEY ? 2 : 1,
    dropIdx: fromIdx,
  };
  document.body.style.cursor = 'grabbing';

  const rect = chartFrame.getBoundingClientRect();
  const y = ev.clientY - rect.top;
  updateDropTarget(y);

  window.addEventListener('mousemove', onLaneDragMove);
  window.addEventListener('mouseup', endLaneDrag);
}

function onLaneDragMove(ev) {
  if (!laneDragState) return;
  const rect = chartFrame.getBoundingClientRect();
  const y = ev.clientY - rect.top;
  updateDropTarget(y);
}

function updateDropTarget(y) {
  const dropIdx = computeDropIndex(y);
  laneDragState.dropIdx = dropIdx;
  const insertY = computeInsertY(dropIdx);
  dropIndicator.style.top = `${insertY}px`;
  dropIndicator.hidden = false;
}

// Given a cursor Y (chart-frame coords), find the insertion index in
// state.laneOrder where the dragged lane should land. An index of N
// means "after the last lane". The gene rows are treated as one slot.
function computeDropIndex(y) {
  const yOffset = svgTopInFrame();
  const order = state.laneOrder;
  let idx = order.length;
  for (let i = 0; i < order.length; i++) {
    if (order[i].kind === 'gene-') continue;
    const r = layoutState.laneRects.get(order[i].key);
    if (!r) continue;
    const geneNegative = order[i].kind === 'gene+' ? layoutState.laneRects.get('gene-') : null;
    const bottom = geneNegative ? geneNegative.top + geneNegative.height : r.top + r.height;
    const midY = yOffset + (r.top + bottom) / 2;
    if (y < midY) { idx = i; break; }
  }
  return idx === 1 ? 2 : idx;
}

// Y (in chart-frame coords) where the drop-indicator line should render
// for a given insertion index — halfway between the neighbouring lanes.
function computeInsertY(idx) {
  const yOffset = svgTopInFrame();
  const order = state.laneOrder;
  if (order.length === 0) return yOffset;
  if (idx <= 0) {
    const first = layoutState.laneRects.get(order[0].key);
    return yOffset + first.top - LANE_GAP / 2;
  }
  if (idx >= order.length) {
    const last = layoutState.laneRects.get(order[order.length - 1].key);
    return yOffset + last.top + last.height + LANE_GAP / 2;
  }
  const prev = layoutState.laneRects.get(order[idx - 1].key);
  const next = layoutState.laneRects.get(order[idx].key);
  return yOffset + (prev.top + prev.height + next.top) / 2;
}

function endLaneDrag() {
  window.removeEventListener('mousemove', onLaneDragMove);
  window.removeEventListener('mouseup', endLaneDrag);
  document.body.style.cursor = '';
  dropIndicator.hidden = true;

  const st = laneDragState;
  laneDragState = null;
  if (!st) return;

  const handleEls = layoutState.handleEls?.get(st.laneKey) ?? [];
  for (const el of handleEls) el.classList.remove('dragging');

  const { fromIdx, dropIdx, groupSize } = st;
  // A grouped gene drag occupies two adjacent insertion slots.
  if (dropIdx === fromIdx || dropIdx === fromIdx + groupSize) {
    hideAllHandles();
    return;
  }
  const items = state.laneOrder.splice(fromIdx, groupSize);
  const adjusted = dropIdx > fromIdx ? dropIdx - groupSize : dropIdx;
  state.laneOrder.splice(adjusted, 0, ...items);
  render();
  buildSignalControls();
  hideAllHandles();
}
