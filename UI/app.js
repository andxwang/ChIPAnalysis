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

const state = {
  genes: [],
  peaks: [],
  dataMin: 0,
  dataMax: 1,
  zoom: 1,
};

// ---------------------------------------------------------------------------
// Data loading
// ---------------------------------------------------------------------------

async function loadDataFromFiles() {
  const genesFile = document.getElementById('genes-file').files[0];
  const peaksFile = document.getElementById('peaks-file').files[0];

  if (!genesFile || !peaksFile) {
    status.textContent = 'Please select both files first.';
    return;
  }

  try {
    const [genesText, peaksText] = await Promise.all([
      readFileSmart(genesFile),
      readFileSmart(peaksFile),
    ]);

    state.genes = parseGenes(genesText);
    state.peaks = parsePeaks(peaksText);

    if (state.genes.length === 0 && state.peaks.length === 0) {
      status.textContent = 'No features found in the provided files.';
      return;
    }

    const { min, max } = computeDataExtent(state.genes, state.peaks);
    // Pad by 1% each side so features never sit flush against the edges.
    const pad = Math.max(1, Math.round((max - min) * 0.01));
    state.dataMin = Math.max(0, min - pad);
    state.dataMax = max + pad;
    state.zoom = clampZoom(ZOOM_STEP ** 8);

    render();
    container.scrollLeft = 0;
    status.textContent = `Loaded ${state.genes.length} genes and ${state.peaks.length} peaks.`;
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

function extractName(attr) {
  const match = attr.match(/name=([^;"]+)/i);
  return match ? match[1].trim() : 'unknown';
}

function computeDataExtent(genes, peaks) {
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
  const height = 340;
  const margin = { top: 24, right: 24, bottom: 46, left: 24 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.setAttribute('width', width);
  svg.setAttribute('height', height);
  svg.style.width = `${width}px`;
  svg.style.height = `${height}px`;
  svg.replaceChildren();

  const xOf = (coord) =>
    margin.left + ((coord - state.dataMin) / dataSpan) * plotWidth;

  // Lane geometry: forward strand above axis, reverse strand below.
  const axisY = margin.top + plotHeight / 2;
  const laneHeight = 22;
  const laneGap = 10;
  const forwardY = axisY - laneGap - laneHeight;
  const reverseY = axisY + laneGap;
  const peakY = axisY + laneGap + laneHeight + 14;
  const peakHeight = 18;

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
  if (state.genes.length === 0 && state.peaks.length === 0) return;
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