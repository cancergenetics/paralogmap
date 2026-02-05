const state = {
  aliasMaps: null,
  pairBuckets: null,
  diseaseBuckets: null,
  cellLineDisplay: new Map(),
  cellLineLineage: new Map(),
  cellLineDisease: new Map(),
  geneMeta: new Map(),
  cellLineIndex: [],
  geneIndex: [],
  pairIndex: [],
  diseaseIndex: [],
  topPairsByCellLine: null,
  topZPairsByCellLine: null,
  bucketCache: new Map(),
  bucketZCache: new Map(),
  geneBucketCache: new Map(),
  diseaseBucketCache: new Map(),
  suggestions: {
    cell: [],
    gene: [],
    pair: [],
    disease: [],
  },
};

const statusEl = document.getElementById("status");
const resultsEl = document.getElementById("results");
const inputEl = document.getElementById("searchInput");
const buttonEl = document.getElementById("searchBtn");
const suggestionsEl = document.getElementById("suggestions");
const searchRowEl = document.getElementById("searchRow");
const diseaseSelectWrapEl = document.getElementById("diseaseSelectWrap");
const diseaseSelectEl = document.getElementById("diseaseSelect");
const aboutBtn = document.getElementById("aboutBtn");
const aboutPanel = document.getElementById("aboutPanel");
const tabEls = document.querySelectorAll(".tab");

let currentMode = "cell";

function normalizeAlias(value) {
  return value.toUpperCase().replace(/[^A-Z0-9]/g, "");
}

function geneBucket(symbol) {
  const norm = normalizeAlias(symbol);
  if (!norm) return "00";
  const bucket = norm.charCodeAt(0) % 64;
  return String(bucket).padStart(2, "0");
}

function formatScore(value) {
  if (value === null || value === undefined) return "";
  return Number(value).toFixed(4);
}

function formatZ(value) {
  if (value === null || value === undefined) return "";
  return Number(value).toFixed(2);
}

function formatP(value) {
  if (value === null || value === undefined) return "";
  const num = Number(value);
  if (!Number.isFinite(num)) return "";
  if (num === 0) return "0";
  if (num < 0.001) return num.toExponential(2);
  return num.toFixed(3);
}

function cellToText(cell) {
  if (cell instanceof Node) {
    return cell.getAttribute("data-value") || cell.textContent || "";
  }
  if (cell === null || cell === undefined) return "";
  return String(cell);
}

function downloadCsv(filename, headers, rows) {
  const lines = [];
  const esc = (value) => {
    const text = cellToText(value);
    if (/[",\n]/.test(text)) {
      return `"${text.replace(/"/g, '""')}"`;
    }
    return text;
  };
  lines.push(headers.map(esc).join(","));
  rows.forEach((row) => {
    lines.push(row.map(esc).join(","));
  });
  const blob = new Blob([lines.join("\n")], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}

function setStatus(message) {
  statusEl.textContent = message;
  statusEl.style.display = "block";
}

function clearStatus() {
  statusEl.textContent = "";
  statusEl.style.display = "none";
}

function clearResults() {
  resultsEl.innerHTML = "";
}

function renderSection(title, meta, headers, rows, downloadConfig = null) {
  const section = document.createElement("div");
  section.className = "section";

  const headerRow = document.createElement("div");
  headerRow.className = "section-header";
  const h2 = document.createElement("h2");
  h2.textContent = title;
  headerRow.appendChild(h2);

  if (downloadConfig) {
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "download-btn";
    btn.textContent = "Download CSV";
    btn.addEventListener("click", () => {
      downloadCsv(downloadConfig.filename, downloadConfig.headers, downloadConfig.rows);
    });
    headerRow.appendChild(btn);
  }

  section.appendChild(headerRow);

  if (meta) {
    const metaEl = document.createElement("div");
    metaEl.className = "meta";
    if (typeof meta === "string") {
      metaEl.textContent = meta;
    } else {
      metaEl.appendChild(meta);
    }
    section.appendChild(metaEl);
  }

  const list = document.createElement("div");
  list.className = "list";

  const listHeaderRow = document.createElement("div");
  listHeaderRow.className = "row header";
  headers.forEach((label) => {
    const cell = document.createElement("div");
    cell.textContent = label;
    listHeaderRow.appendChild(cell);
  });
  list.appendChild(listHeaderRow);

  rows.forEach((row) => {
    const rowEl = document.createElement("div");
    rowEl.className = "row";
    row.forEach((cellValue, index) => {
      const cell = document.createElement("div");
      if (cellValue instanceof Node) {
        cell.appendChild(cellValue);
      } else {
        cell.textContent = cellValue;
      }
      if (index === 0) {
        cell.classList.add("mono");
      }
      rowEl.appendChild(cell);
    });
    list.appendChild(rowEl);
  });

  section.appendChild(list);
  return section;
}

function formatCellLine(depmapId) {
  const display = state.cellLineDisplay.get(depmapId);
  if (display && display !== depmapId) {
    return `${display} (${depmapId})`;
  }
  return depmapId;
}

function formatLineage(depmapId) {
  return state.cellLineLineage.get(depmapId) || "Unknown";
}

function formatPrimaryDisease(depmapId) {
  return state.cellLineDisease.get(depmapId) || "Unknown";
}

function geneDisplay(symbol) {
  const meta = state.geneMeta.get(symbol);
  if (meta && meta.name) {
    return `${symbol} — ${meta.name}`;
  }
  return symbol;
}

function geneFullName(symbol) {
  const meta = state.geneMeta.get(symbol);
  return (meta && meta.name) || "";
}

function geneFullNameLine(symbol) {
  const name = geneFullName(symbol);
  return name ? `${symbol}: ${name}` : symbol;
}

function createHgncLink(symbol) {
  const meta = state.geneMeta.get(symbol);
  if (!meta || !meta.hgnc_id) {
    return null;
  }
  const link = document.createElement("a");
  link.className = "link";
  link.href = `https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/${meta.hgnc_id}`;
  link.target = "_blank";
  link.rel = "noopener";
  link.textContent = `${symbol} (HGNC)`;
  return link;
}

function createPairAnchor(pairId) {
  const link = document.createElement("a");
  link.className = "link";
  link.href = "#";
  link.textContent = pairId;
  link.addEventListener("click", (event) => {
    event.preventDefault();
    setMode("pair", { skipSuggestions: true });
    inputEl.value = pairId;
    renderPair(pairId);
  });
  return link;
}

function createCellLineAnchor(depmapId) {
  const link = document.createElement("a");
  link.className = "link";
  link.href = "#";
  link.textContent = formatCellLine(depmapId);
  link.addEventListener("click", (event) => {
    event.preventDefault();
    setMode("cell", { skipSuggestions: true });
    inputEl.value = depmapId;
    renderCellLine(depmapId);
  });
  return link;
}

function showSuggestions(items) {
  suggestionsEl.innerHTML = "";
  if (!items.length) {
    suggestionsEl.classList.remove("active");
    return;
  }
  items.forEach((item) => {
    const div = document.createElement("div");
    div.className = "suggestion-item";
    div.textContent = item.label;
    div.addEventListener("click", () => {
      inputEl.value = item.value;
      suggestionsEl.classList.remove("active");
      handleSearch();
    });
    suggestionsEl.appendChild(div);
  });
  suggestionsEl.classList.add("active");
}

function hideSuggestions() {
  suggestionsEl.classList.remove("active");
  suggestionsEl.innerHTML = "";
}

function buildSuggestionEntry({ value, label, aliases }) {
  return {
    value,
    label,
    aliases: aliases.map(normalizeAlias),
  };
}

function updateSuggestions() {
  const query = inputEl.value.trim();
  if (currentMode === "disease") {
    hideSuggestions();
    return;
  }
  if (!query) {
    hideSuggestions();
    return;
  }
  const norm = normalizeAlias(query);
  if (!norm) {
    hideSuggestions();
    return;
  }
  const list = state.suggestions[currentMode] || [];
  const matches = [];
  for (const item of list) {
    if (item.aliases.some((alias) => alias.startsWith(norm))) {
      matches.push(item);
    }
    if (matches.length >= 10) break;
  }
  showSuggestions(matches);
}

function updateSearchVisibility(mode) {
  document.body.classList.toggle("mode-disease", mode === "disease");
  if (!searchRowEl || !diseaseSelectWrapEl) return;
  if (mode === "disease") {
    diseaseSelectWrapEl.removeAttribute("hidden");
    hideSuggestions();
  } else {
    diseaseSelectWrapEl.setAttribute("hidden", "");
  }
}

function buildDiseaseSelect() {
  if (!diseaseSelectEl) return;
  diseaseSelectEl.innerHTML = "";
  const placeholder = document.createElement("option");
  placeholder.value = "";
  placeholder.textContent = "Select a disease…";
  diseaseSelectEl.appendChild(placeholder);
  state.diseaseIndex.forEach((entry) => {
    const opt = document.createElement("option");
    opt.value = entry.disease;
    opt.textContent = entry.disease;
    diseaseSelectEl.appendChild(opt);
  });
}

async function fetchJson(path) {
  const res = await fetch(path);
  if (!res.ok) {
    throw new Error(`Failed to load ${path}`);
  }
  return res.json();
}

async function ensureCellLineData() {
  if (state.topPairsByCellLine && state.topZPairsByCellLine) return;
  const [topPairs, topZ] = await Promise.all([
    fetchJson("dist/top/top_pairs_by_cell_line.json"),
    fetchJson("dist/top/top_zpairs_by_cell_line.json"),
  ]);
  state.topPairsByCellLine = topPairs;
  state.topZPairsByCellLine = topZ;
}

async function loadBucket(bucket) {
  if (!state.bucketCache.has(bucket)) {
    const data = await fetchJson(`dist/top/pair_topcells_bucket_${bucket}.json`);
    state.bucketCache.set(bucket, data);
  }
  return state.bucketCache.get(bucket);
}

async function loadBucketZ(bucket) {
  if (!state.bucketZCache.has(bucket)) {
    const data = await fetchJson(`dist/top/pair_topcells_z_bucket_${bucket}.json`);
    state.bucketZCache.set(bucket, data);
  }
  return state.bucketZCache.get(bucket);
}

async function loadGeneBucket(bucket) {
  if (!state.geneBucketCache.has(bucket)) {
    const data = await fetchJson(`dist/agg/gene_pairs_by_mean_bucket_${bucket}.json`);
    state.geneBucketCache.set(bucket, data);
  }
  return state.geneBucketCache.get(bucket);
}

async function loadDiseaseBucket(bucket) {
  if (!state.diseaseBucketCache.has(bucket)) {
    const data = await fetchJson(`dist/disease/disease_pairs_bucket_${bucket}.json`);
    state.diseaseBucketCache.set(bucket, data);
  }
  return state.diseaseBucketCache.get(bucket);
}

async function renderCellLine(depmapId) {
  await ensureCellLineData();
  clearResults();
  clearStatus();
  hideSuggestions();

  const display = formatCellLine(depmapId);
  const topPairs = state.topPairsByCellLine[depmapId] || [];
  const topZPairs = state.topZPairsByCellLine[depmapId] || [];
  const info = document.createElement("div");
  info.className = "info";
  const label = document.createElement("strong");
  label.textContent = `Cell line: ${display}`;
  info.appendChild(label);

  const lineage = document.createElement("span");
  lineage.textContent = `Lineage: ${formatLineage(depmapId)}`;
  info.appendChild(lineage);

  const disease = document.createElement("span");
  disease.textContent = `Primary disease: ${formatPrimaryDisease(depmapId)}`;
  info.appendChild(disease);

  const depmapLink = document.createElement("a");
  depmapLink.className = "link";
  depmapLink.href = `https://depmap.org/portal/cell_line/${depmapId}?tab=overview`;
  depmapLink.target = "_blank";
  depmapLink.rel = "noopener";
  depmapLink.textContent = "Open in DepMap";
  info.appendChild(depmapLink);

  resultsEl.appendChild(info);

  const grid = document.createElement("div");
  grid.className = "grid-two";
  grid.appendChild(
    renderSection(
      "Top 100 gene pairs by raw score",
      null,
      ["Pair", "Score", ""],
      topPairs.map((item) => [createPairAnchor(item.pair_id), formatScore(item.score), ""]),
      {
        filename: `${depmapId}_top_pairs_by_score.csv`,
        headers: ["pair_id", "score"],
        rows: topPairs.map((item) => [item.pair_id, item.score]),
      }
    )
  );

  grid.appendChild(
    renderSection(
      "Top 100 gene pairs by differential z-score",
      null,
      ["Pair", "Score", "Z"],
      topZPairs.map((item) => [createPairAnchor(item.pair_id), formatScore(item.score), formatZ(item.z)]),
      {
        filename: `${depmapId}_top_pairs_by_z.csv`,
        headers: ["pair_id", "score", "z"],
        rows: topZPairs.map((item) => [item.pair_id, item.score, item.z]),
      }
    )
  );

  resultsEl.appendChild(grid);
}

async function renderPair(pairId) {
  const bucket = state.pairBuckets[pairId];
  if (!bucket) {
    setStatus("No bucket found for this pair.");
    return;
  }

  clearResults();
  hideSuggestions();
  setStatus("Loading pair buckets…");

  const byScore = await loadBucket(bucket);

  clearStatus();

  const listScore = byScore[pairId] || [];

  const info = document.createElement("div");
  info.className = "info";
  const label = document.createElement("strong");
  label.textContent = `Pair: ${pairId}`;
  info.appendChild(label);

  const [gene1, gene2] = pairId.split("__");
  const names = document.createElement("span");
  names.textContent = `${geneFullNameLine(gene1)} • ${geneFullNameLine(gene2)}`;
  info.appendChild(names);
  const link1 = createHgncLink(gene1);
  const link2 = createHgncLink(gene2);
  if (link1) info.appendChild(link1);
  if (link2) info.appendChild(link2);

  resultsEl.appendChild(info);

  resultsEl.appendChild(
    renderSection(
      "Top 100 cell lines by raw score",
      null,
      ["Cell line", "Lineage", "Score"],
      listScore.map((item) => [
        createCellLineAnchor(item.depmap_id),
        formatLineage(item.depmap_id),
        formatScore(item.score),
      ]),
      {
        filename: `${pairId}_top_cell_lines_by_score.csv`,
        headers: ["depmap_id", "lineage", "score"],
        rows: listScore.map((item) => [
          item.depmap_id,
          formatLineage(item.depmap_id),
          item.score,
        ]),
      }
    )
  );
}

async function renderGene(symbol) {
  const bucket = geneBucket(symbol);
  clearResults();
  setStatus("Loading gene bucket…");
  hideSuggestions();

  const data = await loadGeneBucket(bucket);
  clearStatus();

  const list = (data && data[symbol]) || [];
  const info = document.createElement("div");
  info.className = "info";
  const label = document.createElement("strong");
  label.textContent = `Gene: ${symbol}`;
  info.appendChild(label);
  const meta = state.geneMeta.get(symbol);
  if (meta && meta.name) {
    const name = document.createElement("span");
    name.textContent = `Full name: ${meta.name}`;
    info.appendChild(name);
  }
  const link = createHgncLink(symbol);
  if (link) info.appendChild(link);
  resultsEl.appendChild(info);

  resultsEl.appendChild(
    renderSection(
      "Gene pairs sorted by mean score",
      null,
      ["Pair", "Partner", "Mean"],
      list.map((item) => [createPairAnchor(item.pair_id), item.partner, formatScore(item.mean)]),
      {
        filename: `${symbol}_pairs_by_mean.csv`,
        headers: ["pair_id", "partner", "mean"],
        rows: list.map((item) => [item.pair_id, item.partner, item.mean]),
      }
    )
  );
}

async function renderDisease(disease) {
  const bucket = state.diseaseBuckets[disease];
  if (!bucket) {
    setStatus("No bucket found for this disease.");
    return;
  }
  clearResults();
  hideSuggestions();
  setStatus("Loading disease results…");

  const data = await loadDiseaseBucket(bucket);
  clearStatus();
  if (diseaseSelectEl) {
    diseaseSelectEl.value = disease;
  }

  const list = (data && data[disease]) || [];
  const info = document.createElement("div");
  info.className = "info";
  const label = document.createElement("strong");
  label.textContent = `Disease: ${disease}`;
  info.appendChild(label);
  resultsEl.appendChild(info);

  resultsEl.appendChild(
    renderSection(
      "Disease-specific paralog pair dependencies",
      null,
      ["Pair", "Median (disease)", "Median (other)", "Diff", "P-value"],
      list.map((item) => [
        createPairAnchor(item.pair_id),
        formatScore(item.median_disease),
        formatScore(item.median_other),
        formatScore(item.diff),
        formatP(item.p_value),
      ]),
      {
        filename: `${disease.replace(/\\s+/g, "_")}_disease_pairs.csv`,
        headers: ["pair_id", "median_disease", "median_other", "diff", "p_value"],
        rows: list.map((item) => [
          item.pair_id,
          item.median_disease,
          item.median_other,
          item.diff,
          item.p_value,
        ]),
      }
    )
  );
}

async function handleSearch() {
  hideSuggestions();
  const query = inputEl.value.trim();
  if (!query) {
    setStatus("Enter a search term to begin.");
    return;
  }
  if (!state.aliasMaps) {
    setStatus("Index not loaded yet.");
    return;
  }

  const norm = normalizeAlias(query);
  let resolved = null;

  if (currentMode === "cell") {
    resolved = state.aliasMaps.cell_lines[norm];
    if (!resolved) {
      setStatus("No matching cell line alias found.");
      return;
    }
    setStatus("Loading cell line results…");
    await renderCellLine(resolved);
    return;
  }

  if (currentMode === "pair") {
    resolved = state.aliasMaps.pairs[norm];
    if (!resolved) {
      setStatus("No matching gene pair alias found.");
      return;
    }
    await renderPair(resolved);
    return;
  }

  if (currentMode === "gene") {
    resolved = state.aliasMaps.genes[norm];
    if (!resolved) {
      setStatus("No matching gene alias found.");
      return;
    }
    await renderGene(resolved);
    return;
  }

  if (currentMode === "disease") {
    setStatus("Select a disease from the dropdown above.");
    return;
  }
}

function setMode(mode, options = {}) {
  currentMode = mode;
  tabEls.forEach((tab) => {
    tab.classList.toggle("active", tab.dataset.mode === mode);
  });
  updateSearchVisibility(mode);
  if (!options.skipStatus) {
    clearResults();
    if (mode === "disease") {
      setStatus("Select a disease from the dropdown above.");
    } else {
      setStatus("Enter a search term to begin.");
    }
  }
  hideSuggestions();
  if (!options.skipSuggestions) {
    updateSuggestions();
  }
  if (mode === "disease" && !options.skipAutoRender && diseaseSelectEl && diseaseSelectEl.value) {
    renderDisease(diseaseSelectEl.value);
  }
}

async function init() {
  try {
    const [aliasMaps, cellLines, pairBuckets, genes, pairs, diseases, diseaseBuckets] = await Promise.all([
      fetchJson("dist/index/alias_maps.json"),
      fetchJson("dist/index/cell_lines.json"),
      fetchJson("dist/index/pair_buckets.json"),
      fetchJson("dist/index/genes.json"),
      fetchJson("dist/index/pairs.json"),
      fetchJson("dist/index/diseases.json"),
      fetchJson("dist/index/disease_buckets.json"),
    ]);

    state.aliasMaps = aliasMaps;
    state.pairBuckets = pairBuckets;
    state.diseaseBuckets = diseaseBuckets;
    state.cellLineIndex = cellLines;
    state.geneIndex = genes;
    state.pairIndex = pairs;
    state.diseaseIndex = diseases;
    state.cellLineDisplay = new Map(
      cellLines.map((entry) => [entry.depmap_id, entry.display])
    );
    state.cellLineLineage = new Map(
      cellLines.map((entry) => [entry.depmap_id, entry.lineage || ""])
    );
    state.cellLineDisease = new Map(
      cellLines.map((entry) => [entry.depmap_id, entry.primary_disease || ""])
    );
    state.geneMeta = new Map(
      genes.map((entry) => [entry.symbol, { name: entry.name, hgnc_id: entry.hgnc_id }])
    );

    state.suggestions.cell = cellLines.map((entry) =>
      buildSuggestionEntry({
        value: entry.depmap_id,
        label: entry.display && entry.display !== entry.depmap_id
          ? `${entry.display} (${entry.depmap_id})`
          : entry.depmap_id,
        aliases: entry.aliases || [entry.depmap_id],
      })
    );

    state.suggestions.gene = genes.map((entry) =>
      buildSuggestionEntry({
        value: entry.symbol,
        label: entry.name ? `${entry.symbol} — ${entry.name}` : entry.symbol,
        aliases: entry.aliases || [entry.symbol],
      })
    );

    state.suggestions.pair = pairs.map((entry) =>
      buildSuggestionEntry({
        value: entry.pair_id,
        label: entry.pair_id,
        aliases: [
          entry.pair_id,
          `${entry.gene1}_${entry.gene2}`,
          `${entry.gene2}_${entry.gene1}`,
          `${entry.gene2}__${entry.gene1}`,
        ],
      })
    );

    state.suggestions.disease = diseases.map((entry) =>
      buildSuggestionEntry({
        value: entry.disease,
        label: entry.disease,
        aliases: entry.aliases || [entry.disease],
      })
    );

    buildDiseaseSelect();
    updateSearchVisibility(currentMode);

    setStatus("Enter a search term to begin.");
  } catch (err) {
    console.error(err);
    setStatus("Failed to load index files. Check dist/ output.");
  }
}

buttonEl.addEventListener("click", handleSearch);
inputEl.addEventListener("keydown", (event) => {
  if (event.key === "Enter") {
    handleSearch();
  }
});
inputEl.addEventListener("input", updateSuggestions);

tabEls.forEach((tab) => {
  tab.addEventListener("click", () => setMode(tab.dataset.mode));
});

if (diseaseSelectEl) {
  diseaseSelectEl.addEventListener("change", () => {
    const disease = diseaseSelectEl.value;
    if (!disease) return;
    setMode("disease", { skipSuggestions: true, skipStatus: true, skipAutoRender: true });
    renderDisease(disease);
  });
}

if (aboutBtn && aboutPanel) {
  aboutBtn.addEventListener("click", () => {
    const isHidden = aboutPanel.hasAttribute("hidden");
    if (isHidden) {
      aboutPanel.removeAttribute("hidden");
    } else {
      aboutPanel.setAttribute("hidden", "");
    }
  });
}

init();
