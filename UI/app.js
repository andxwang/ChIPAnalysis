const genesPath = '../MabATCC19977_gff.gff';
const peaksPath = '../SigHP1/SigHP1_FDR0.01combo.gff';

const svg = document.getElementById('chart');
const status = document.getElementById('status');
const container = document.getElementById('chart-container');

const state = {
  genes: [],
  peaks: [],
  viewStart: 0,
  viewEnd: 0,
  zoom: 1,
  initialWindow: 250000,
  isSyncingScroll: false,
};

function getMaxCoord() {
  return Math.max(
    ...state.genes.flatMap((g) => [g.start, g.end]),
    ...state.peaks.flatMap((p) => [p.start, p.end]),
  );
}

async function loadData() {
  const [genesText, peaksText] = await Promise.all([
    fetch(genesPath).then((res) => res.text()),
    fetch(peaksPath).then((res) => res.text()),
  ]);

  state.genes = parseGenes(genesText);
  state.peaks = parsePeaks(peaksText);

  const maxCoord = getMaxCoord();
  const initialWindow = Math.max(250000, Math.round(maxCoord * 0.1));

  state.initialWindow = initialWindow;
  state.viewStart = 0;
  state.viewEnd = initialWindow;
  state.zoom = 1;

  render();
  status.textContent = `Loaded ${state.genes.length} genes and ${state.peaks.length} peaks`;
}

function parseGenes(text) {
  return text
    .trim()
    .split(/\r?\n/)
    .filter((line) => line && !line.startsWith('#'))
    .map((line) => {
      const cols = line.split('\t');
      const start = Number(cols[3]);
      const end = Number(cols[4]);
      const direction = cols[6];
      const attr = cols[8] || '';
      const name = extractName(attr);
      return { start, end, direction, name };
    })
    .filter((row) => Number.isFinite(row.start) && Number.isFinite(row.end));
}

function parsePeaks(text) {
  return text
    .trim()
    .split(/\r?\n/)
    .filter((line) => line && !line.startsWith('#'))
    .map((line) => {
      const cols = line.split('\t');
      const start = Number(cols[3]);
      const end = Number(cols[4]);
      const score = Number(cols[5]);
      return { start, end, score };
    })
    .filter((row) => Number.isFinite(row.start) && Number.isFinite(row.end));
}

function extractName(attr) {
  const match = attr.match(/name=([^;]+)/);
  return match ? match[1] : 'unknown';
}

function getNiceTickStep(span, targetTicks = 10) {
  if (span <= 0) return 1;

  const roughStep = span / targetTicks;
  const magnitude = 10 ** Math.floor(Math.log10(roughStep));
  const candidates = [1, 2, 2.5, 5, 10].map((value) => value * magnitude);

  let bestStep = candidates[0];
  let bestDiff = Number.POSITIVE_INFINITY;

  candidates.forEach((step) => {
    const diff = Math.abs(step - roughStep);
    if (diff < bestDiff || (diff === bestDiff && step < bestStep)) {
      bestStep = step;
      bestDiff = diff;
    }
  });

  return bestStep;
}

function getTickValues(start, end, targetTicks = 10) {
  const span = Math.max(1, end - start);
  const step = getNiceTickStep(span, targetTicks);
  const first = Math.floor(start / step) * step;
  const last = Math.ceil(end / step) * step;
  const ticks = [];
  for (let value = first; value <= last; value += step) {
    ticks.push(value);
  }
  return ticks;
}

function getChartWidth() {
  const viewportWidth = Math.max(900, container.clientWidth || 1200);
  const span = Math.max(1, state.viewEnd - state.viewStart);
  const minViewportPadding = 1400;
  const spanBasedWidth = Math.round(viewportWidth + Math.max(800, span / 40));
  return Math.max(viewportWidth + minViewportPadding, Math.round(viewportWidth * 2.6), spanBasedWidth);
}

function syncScrollPosition() {
  const maxScroll = Math.max(0, getChartWidth() - container.clientWidth);
  if (maxScroll <= 0) return;

  const span = Math.max(1, state.viewEnd - state.viewStart);
  const maxCoord = Math.max(1, getMaxCoord());
  const maxStart = Math.max(0, maxCoord - span);
  const ratio = maxStart > 0 ? state.viewStart / maxStart : 0;
  const nextScrollLeft = Math.max(0, Math.min(maxScroll, ratio * maxScroll));

  state.isSyncingScroll = true;
  container.scrollLeft = nextScrollLeft;
  requestAnimationFrame(() => {
    state.isSyncingScroll = false;
  });
}

function render({ preserveScroll = false } = {}) {
  const span = Math.max(1, state.viewEnd - state.viewStart);
  const width = getChartWidth();
  const height = 520;
  const margin = { top: 40, right: 40, bottom: 70, left: 70 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;

  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.setAttribute('width', width);
  svg.setAttribute('height', height);
  svg.style.width = `${width}px`;
  svg.style.maxWidth = 'none';
  svg.style.minWidth = `${width}px`;
  svg.innerHTML = '';

  const x = (coord) => margin.left + ((coord - state.viewStart) / span) * plotWidth;

  const axisHeight = margin.top + plotHeight / 2;
  const targetTicks = Math.max(6, Math.min(12, Math.floor(plotWidth / 110)));
  const ticks = getTickValues(state.viewStart, state.viewEnd, targetTicks);
  ticks.forEach((coord) => {
    const xPos = x(coord);
    const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    line.setAttribute('x1', xPos);
    line.setAttribute('x2', xPos);
    line.setAttribute('y1', margin.top);
    line.setAttribute('y2', margin.top + plotHeight);
    line.setAttribute('class', 'grid-line');
    svg.appendChild(line);

    const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    label.setAttribute('x', xPos);
    label.setAttribute('y', margin.top + plotHeight + 24);
    label.setAttribute('text-anchor', 'middle');
    label.setAttribute('class', 'tick-label');
    label.textContent = coord.toLocaleString();
    svg.appendChild(label);
  });

  const axis = document.createElementNS('http://www.w3.org/2000/svg', 'line');
  axis.setAttribute('x1', margin.left);
  axis.setAttribute('x2', margin.left + plotWidth);
  axis.setAttribute('y1', axisHeight);
  axis.setAttribute('y2', axisHeight);
  axis.setAttribute('class', 'axis-line');
  svg.appendChild(axis);

  const axisLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
  axisLabel.setAttribute('x', margin.left + plotWidth / 2);
  axisLabel.setAttribute('y', height - 16);
  axisLabel.setAttribute('text-anchor', 'middle');
  axisLabel.setAttribute('class', 'axis-label');
  axisLabel.textContent = 'Coordinate (bp)';
  svg.appendChild(axisLabel);

  const geneTopY = axisHeight - 34;
  const geneBottomY = axisHeight + 8;
  const peakY = axisHeight + 58;

  state.genes.forEach((gene) => {
    const geneStart = x(gene.start);
    const geneEnd = x(gene.end);
    const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    const rectWidth = Math.max(3, Math.abs(geneEnd - geneStart));
    const rowY = gene.direction === '+' ? geneTopY - 10 : geneBottomY - 10;
    rect.setAttribute('x', Math.min(geneStart, geneEnd));
    rect.setAttribute('y', rowY);
    rect.setAttribute('width', rectWidth);
    rect.setAttribute('height', 20);
    rect.setAttribute('rx', 5);
    rect.setAttribute('class', 'gene-bar');
    svg.appendChild(rect);

    const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    label.setAttribute('x', Math.min(geneStart, geneEnd) + 4);
    label.setAttribute('y', rowY - 6);
    label.setAttribute('class', 'gene-label');
    label.textContent = gene.name;
    svg.appendChild(label);
  });

  state.peaks.forEach((peak) => {
    const peakStart = x(peak.start);
    const peakEnd = x(peak.end);
    const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    const rectWidth = Math.max(3, Math.abs(peakEnd - peakStart));
    rect.setAttribute('x', Math.min(peakStart, peakEnd));
    rect.setAttribute('y', peakY - 10);
    rect.setAttribute('width', rectWidth);
    rect.setAttribute('height', 20);
    rect.setAttribute('rx', 5);
    rect.setAttribute('class', 'peak-bar');
    svg.appendChild(rect);

    const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    label.setAttribute('x', Math.min(peakStart, peakEnd) + 4);
    label.setAttribute('y', peakY + 24);
    label.setAttribute('class', 'peak-label');
    label.textContent = `Peak ${peak.score.toFixed(3)}`;
    svg.appendChild(label);
  });

  if (!preserveScroll) {
    syncScrollPosition();
  }
}

function panBy(direction) {
  const span = state.viewEnd - state.viewStart;
  const maxCoord = getMaxCoord();
  const step = Math.max(1000, span / 3);
  if (direction === 'left') {
    state.viewStart = Math.max(0, state.viewStart - step);
  } else {
    state.viewStart = Math.min(maxCoord - span, state.viewStart + step);
  }
  state.viewEnd = state.viewStart + span;
  render();
}

function zoomBy(direction) {
  const maxCoord = getMaxCoord();
  const span = state.viewEnd - state.viewStart;
  const center = (state.viewStart + state.viewEnd) / 2;

  if (direction === 'in') {
    const nextSpan = Math.max(1000, span * 0.6);
    state.viewStart = Math.max(0, center - nextSpan / 2);
    state.viewEnd = Math.min(maxCoord, center + nextSpan / 2);
  } else {
    const nextSpan = Math.min(maxCoord, Math.max(1000, span * 1.4));
    state.viewStart = Math.max(0, center - nextSpan / 2);
    state.viewEnd = Math.min(maxCoord, center + nextSpan / 2);
  }

  render();
}

function resetView() {
  const maxCoord = getMaxCoord();
  const initialWindow = Math.max(250000, Math.round(maxCoord * 0.1));
  state.initialWindow = initialWindow;
  state.viewStart = 0;
  state.viewEnd = initialWindow;
  render();
}

document.getElementById('zoom-in').addEventListener('click', () => zoomBy('in'));
document.getElementById('zoom-out').addEventListener('click', () => zoomBy('out'));
document.getElementById('reset').addEventListener('click', resetView);
document.getElementById('pan-left').addEventListener('click', () => panBy('left'));
document.getElementById('pan-right').addEventListener('click', () => panBy('right'));

window.addEventListener('resize', render);

container.addEventListener('wheel', (event) => {
  event.preventDefault();
  if (event.deltaY < 0) {
    zoomBy('in');
  } else {
    zoomBy('out');
  }
}, { passive: false });

container.addEventListener('scroll', () => {
  if (state.isSyncingScroll) return;

  const maxScroll = Math.max(0, getChartWidth() - container.clientWidth);
  if (maxScroll <= 0) return;

  const span = Math.max(1, state.viewEnd - state.viewStart);
  const maxCoord = Math.max(1, getMaxCoord());
  const maxStart = Math.max(0, maxCoord - span);
  const ratio = maxScroll > 0 ? container.scrollLeft / maxScroll : 0;
  const nextStart = Math.max(0, Math.min(maxStart, ratio * maxStart));

  if (Math.abs(nextStart - state.viewStart) > 1) {
    state.viewStart = nextStart;
    state.viewEnd = state.viewStart + span;
    render({ preserveScroll: true });
  }
});

loadData();
