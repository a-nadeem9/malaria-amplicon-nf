const STAGES = [
  "prepare_inputs",
  "dada2",
  "prepare_stage2",
  "asv_mapping",
  "prepare_stage3",
  "cigar_check",
  "asv_to_cigar",
  "cdhit",
  "report"
];

const STAGE_LABELS = {
  prepare_inputs: "Checking sample sheet",
  dada2: "Running DADA2",
  prepare_stage2: "Cleaning ASV table",
  asv_mapping: "Mapping ASVs",
  prepare_stage3: "Preparing CIGAR inputs",
  cigar_check: "Checking CIGAR inputs",
  asv_to_cigar: "Converting ASVs to CIGAR",
  cdhit: "Clustering ASV sequences",
  report: "Writing report"
};

const SCAN_STAGE_LABELS = [
  "Discovering FASTQ files",
  "Pairing forward and reverse reads",
  "Parsing sample identifiers",
  "Matching optional metadata",
  "Validating sample records"
];
const SCAN_PAGE_SIZE = 25;
const SAMPLE_PREVIEW_LIMIT = 50;

const STORE_KEY = "malaria-amplicon-nf.flask.settings";
const EMPTY_LOG_TEXT = "> Waiting for run.\n> Live progress appears here while malaria-amplicon-nf is active.";
const ASV_FILTERING_SUMMARY_LABEL = "ASV filtering summary";
let pollTimer = null;
let scanInFlight = false;
let scanReady = false;
let logInFlight = false;
let latestStatusPayload = null;
let latestResultsPayload = null;
let latestAnalysisInputStatus = null;
let lastRunStatus = "";
let completedRedirectKey = "";
let displayedProgressPercent = 0;
let activeOutdir = "";
let pathStyle = "";
let workspaceRoot = "";
const NATIVE_PICKER_TIMEOUT_MS = 120000;
let followLog = true;
let latestDinemitesReadiness = null;
let preferencesSaveTimer = null;
let latestDciferReadiness = null;
let scanReviewRows = [];
let scanReviewFilter = "all";
let scanReviewLibrary = "all";
let scanReviewPage = 1;
let detectedLibraryCounts = {};
let restoredRunLibraries = null;
let metadataContract = {
  schema_version: 1,
  columns: {},
  detection_value_map: {},
  excluded_status_values: []
};
let latestMetadataInspection = null;
let metadataContractValid = true;

function $(selector) {
  return document.querySelector(selector);
}

function text(node, value) {
  if (!node) return;
  node.textContent = value == null ? "" : String(value);
}

function escapeHtml(value) {
  return String(value == null ? "" : value).replace(/[&<>"']/g, (char) => ({
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    "\"": "&quot;",
    "'": "&#39;"
  }[char]));
}

function userMessage(value) {
  let message = String(value || "").trim();
  if (!message) return "";
  message = message.replace(
    /seqtab_cigar\.tsv is required:\s+\S*?run_dada2[\\/]+seqtab_cigar\.tsv/gi,
    "Run the main pipeline first. run_dada2/seqtab_cigar.tsv is required"
  );
  message = message.replace(
    /Sample sheet not found:\s+\S*?([^\\/]+\.csv)\b/gi,
    "Sample sheet not found: $1"
  );
  return message;
}

function terminalLineClass(line) {
  const trimmed = String(line || "").trim();
  if (trimmed.startsWith("[OK]")) return "log-ok";
  if (trimmed.startsWith("[WARN]") || trimmed.startsWith("WARN:")) return "log-warn";
  if (trimmed.startsWith("[ERROR]") || trimmed.startsWith("ERROR") || trimmed.includes("failed")) return "log-error";
  if (trimmed.startsWith(">")) return "log-muted";
  return "";
}

function renderTerminalLog(node, value) {
  if (!node) return;
  node.replaceChildren();
  const lines = String(value == null ? "" : value).split("\n");
  lines.forEach((line, index) => {
    const span = document.createElement("span");
    const lineClass = terminalLineClass(line);
    if (lineClass) span.className = lineClass;
    span.textContent = line;
    node.appendChild(span);
    if (index < lines.length - 1) node.appendChild(document.createTextNode("\n"));
  });
}

function scrollLogToBottom(node) {
  if (!node) return;
  requestAnimationFrame(() => {
    node.scrollTop = node.scrollHeight;
  });
}

function getSettings() {
  try {
    return JSON.parse(localStorage.getItem(STORE_KEY) || "{}");
  } catch (_error) {
    return {};
  }
}

async function loadPersistedSettings() {
  const local = getSettings();
  try {
    const payload = await fetchJson("/api/preferences");
    return {...local, ...(payload.settings || {})};
  } catch (_error) {
    return local;
  }
}

function persistSettings(settings) {
  if (preferencesSaveTimer) window.clearTimeout(preferencesSaveTimer);
  preferencesSaveTimer = window.setTimeout(() => {
    preferencesSaveTimer = null;
    postJson("/api/preferences", {settings}).catch(() => {});
  }, 250);
}

function normalizeDrivePath(value) {
  const raw = String(value || "").trim();
  if (pathStyle !== "wsl") return raw;
  const match = raw.match(/^([A-Za-z]):[\\/]?(.*)$/);
  if (!match) return raw;
  const drive = match[1].toLowerCase();
  const rest = match[2].replaceAll("\\", "/");
  return `/mnt/${drive}/${rest}`;
}

function normalizePathInput(selector) {
  const node = $(selector);
  if (node && node.value) node.value = normalizeDrivePath(node.value);
}

function normalizePathInputs() {
  ["#fastq-dir", "#metadata-path", "#kelt-barcode-map", "#run-samples", "#outdir", "#results-outdir"].forEach(normalizePathInput);
  syncPathTitles();
  updateSamplePathHelp();
}

function syncPathTitles() {
  ["#fastq-dir", "#metadata-path", "#kelt-barcode-map", "#run-samples", "#outdir", "#results-outdir"].forEach((selector) => {
    const node = $(selector);
    if (node) node.title = node.value || "";
  });
}

function isAbsolutePath(value) {
  const raw = String(value || "").trim();
  return raw.startsWith("/") || /^[A-Za-z]:[\\/]/.test(raw);
}

function joinPath(parent, child) {
  const base = String(parent || "").trim();
  const name = String(child || "").trim();
  if (!base) return name;
  const separator = base.includes("\\") && !base.includes("/") ? "\\" : "/";
  if (base === ".") return name.replace(/^[\\/]+/, "");
  return `${base.replace(/[\\/]+$/, "")}${separator}${name.replace(/^[\\/]+/, "")}`;
}

function pathBasename(value, fallback = "") {
  const raw = String(value || "").trim().replace(/[\\/]+$/, "");
  if (!raw) return fallback;
  const parts = raw.split(/[\\/]+/).filter(Boolean);
  return parts.length ? parts[parts.length - 1] : fallback;
}

function pathDirname(value, fallback = ".") {
  const raw = String(value || "").trim();
  if (!raw) return fallback;
  const stripped = raw.replace(/[\\/]+$/, "");
  const index = Math.max(stripped.lastIndexOf("/"), stripped.lastIndexOf("\\"));
  if (index < 0) return fallback;
  if (/^[A-Za-z]:[\\/]?/.test(stripped) && index <= 2) return stripped.slice(0, 3);
  if (index === 0) return stripped.slice(0, 1);
  return stripped.slice(0, index) || fallback;
}

function syncGeneratedSampleSheetPath() {
  const path = activeOutdir ? joinPath(activeOutdir, "samples.csv") : "samples.csv";
  $("#run-samples").value = path;
  return path;
}

function describeSampleSheetPath(value, contextPath, prefix) {
  const path = String(value || "samples.csv").trim() || "samples.csv";
  const context = String(contextPath || "").trim();
  if (isAbsolutePath(path)) return `${prefix}: ${path}`;
  if (context) {
    const normalizedPath = path.replace(/^\.[\\/]+/, "");
    return `${prefix}: ${joinPath(context, normalizedPath)}`;
  }
  return `${prefix}: app workspace/${path}`;
}

function updateSamplePathHelp() {
  const samplePath = activeOutdir ? joinPath(activeOutdir, "samples.csv") : "samples.csv";
  const help = activeOutdir
    ? describeSampleSheetPath(samplePath, workspaceRoot, "Using sample sheet")
    : "samples.csv is created automatically inside the run folder when you start.";
  text(
    $("#run-samples-help"),
    help
  );
}

function analysisMinAbundancePct() {
  return Number($("#analysis-min-abundance-pct")?.value || 1);
}

function analysisAbundanceDenominator() {
  return $("#analysis-abundance-denominator")?.value || "locus";
}

function updateAnalysisFilterSummary() {
  const summary = $("#analysis-filter-summary");
  if (!summary) return;
  const abundance = analysisMinAbundancePct().toLocaleString(undefined, {
    maximumFractionDigits: 2
  });
  const denominator = analysisAbundanceDenominator() === "sample"
    ? "all reads in the sample"
    : "reads at the same locus";
  text(summary, `Keep alleles that pass ${abundance}% of ${denominator} in every technical replicate, merge their counts, then require at least 100 reads at the merged locus.`);
}

function analysisCdhitMode() {
  return $("#dinemites-analysis-mode")?.value || $("#dcifer-analysis-mode")?.value || "off";
}

function normalizeAnalysisMode(mode) {
  if (mode === "summed" || mode === "cdhit98") return "summed";
  return "off";
}

function analysisApiMode() {
  const selected = normalizeAnalysisMode(analysisCdhitMode());
  if (selected === "summed") return "cdhit98";
  return "primary";
}

function syncAnalysisInputChoice(mode) {
  const normalized = normalizeAnalysisMode(mode);
  document.querySelectorAll('input[name="analysis-input-choice"]').forEach((radio) => {
    radio.checked = radio.value === normalized;
  });
  $("#analysis-build-primary")?.classList.toggle("primary", normalized === "off");
  $("#analysis-build-table")?.classList.toggle("primary", normalized === "summed");
}

function setAnalysisMode(mode, persist = true) {
  const normalized = normalizeAnalysisMode(mode);
  syncAnalysisInputChoice(normalized);
  ["#dinemites-analysis-mode", "#dcifer-analysis-mode"].forEach((selector) => {
    const select = $(selector);
    if (select) select.value = normalized;
  });
  const cdhit = normalized === "summed";
  text(
    $("#dinemites-analysis-mode-note"),
    cdhit
      ? "CD-HIT comparison results are saved separately in dinemites_cdhit98."
      : "Primary results are saved in the dinemites folder."
  );
  text(
    $("#dcifer-analysis-mode-note"),
    cdhit
      ? "CD-HIT comparison results are saved separately in dcifer_cdhit98."
      : "Primary results are saved in the dcifer folder."
  );
  text($("#dinemites-run"), cdhit ? "Run CD-HIT comparison" : "Run DINEMITES");
  text($("#dcifer-run"), cdhit ? "Run CD-HIT comparison" : "Run Dcifer");
  if (persist) saveSettings();
}

function setAnalysisModeControlsDisabled(disabled) {
  ["#dinemites-analysis-mode", "#dcifer-analysis-mode"].forEach((selector) => {
    const select = $(selector);
    if (select) select.disabled = disabled;
  });
}

function saveSettings() {
  normalizePathInputs();
  const sharedMinAbundancePct = $("#analysis-min-abundance-pct") ? $("#analysis-min-abundance-pct").value : "1";
  const sharedAbundanceDenominator = $("#analysis-abundance-denominator") ? $("#analysis-abundance-denominator").value : "locus";
  const sharedCdhitMode = analysisCdhitMode();
  const settings = {
    metadataDateOrder: $("#metadata-date-order") ? $("#metadata-date-order").value : "auto",
    fallbackCollectionYear: $("#fallback-collection-year") ? $("#fallback-collection-year").value : "",
    fallbackCollectionDay: $("#fallback-collection-day") ? $("#fallback-collection-day").value : "27",
    resumeRun: false,
    dryRun: false,
    cpus: $("#cpus").value,
    memory: $("#memory").value,
    analysisMinAbundancePct: sharedMinAbundancePct,
    analysisAbundanceDenominator: sharedAbundanceDenominator,
    analysisCdhitMode: sharedCdhitMode,
    dinemitesEnabled: true,
    dinemitesModel: $("#dinemites-model") ? $("#dinemites-model").value : "simple",
    dinemitesNLags: $("#dinemites-n-lags") ? $("#dinemites-n-lags").value : "3",
    dinemitesTLag: $("#dinemites-t-lag") ? $("#dinemites-t-lag").value : "90",
    dinemitesMinAbundancePct: sharedMinAbundancePct,
    dinemitesAbundanceDenominator: sharedAbundanceDenominator,
    dinemitesNoDayCutoff: $("#dinemites-no-day-cutoff") ? $("#dinemites-no-day-cutoff").checked : true,
    dinemitesSeed: $("#dinemites-seed") ? $("#dinemites-seed").value : "1",
    dinemitesRefresh: $("#dinemites-refresh-interval") ? $("#dinemites-refresh-interval").value : "100",
    dinemitesBayesianLagDays: $("#dinemites-bayesian-lag-days") ? $("#dinemites-bayesian-lag-days").value : "30",
    dinemitesBayesianChains: $("#dinemites-bayesian-chains") ? $("#dinemites-bayesian-chains").value : "4",
    dinemitesBayesianParallelChains: $("#dinemites-bayesian-parallel-chains") ? $("#dinemites-bayesian-parallel-chains").value : "2",
    dinemitesBayesianWarmup: $("#dinemites-bayesian-warmup") ? $("#dinemites-bayesian-warmup").value : "500",
    dinemitesBayesianSampling: $("#dinemites-bayesian-sampling") ? $("#dinemites-bayesian-sampling").value : "500",
    dinemitesBayesianAdaptDelta: $("#dinemites-bayesian-adapt-delta") ? $("#dinemites-bayesian-adapt-delta").value : "0.99",
    dinemitesUseSeasonCovariate: $("#dinemites-covariate-season") ? $("#dinemites-covariate-season").checked : true,
    dinemitesUseAgeCovariate: $("#dinemites-covariate-age") ? $("#dinemites-covariate-age").checked : false,
    dinemitesUseGenderCovariate: $("#dinemites-covariate-gender") ? $("#dinemites-covariate-gender").checked : false,
    dinemitesCustomCovariates: $("#dinemites-infection-covariates") ? $("#dinemites-infection-covariates").value : "",
    dinemitesBayesianDropOut: $("#dinemites-bayesian-drop-out") ? $("#dinemites-bayesian-drop-out").checked : false,
    dciferEnabled: true,
    dciferMinAbundancePct: sharedMinAbundancePct,
    dciferAbundanceDenominator: sharedAbundanceDenominator,
    dciferCoiLrank: $("#dcifer-coi-lrank") ? $("#dcifer-coi-lrank").value : "2",
    dciferIbdGridNr: $("#dcifer-ibd-grid-nr") ? $("#dcifer-ibd-grid-nr").value : "1000",
    dciferAlpha: $("#dcifer-alpha") ? $("#dcifer-alpha").value : "0.05"
  };
  localStorage.setItem(STORE_KEY, JSON.stringify(settings));
  persistSettings(settings);
}

function restoreSettings(settings = getSettings()) {
  if (settings.metadataDateOrder && $("#metadata-date-order")) $("#metadata-date-order").value = settings.metadataDateOrder;
  if (settings.fallbackCollectionYear != null && $("#fallback-collection-year")) {
    $("#fallback-collection-year").value = settings.fallbackCollectionYear;
  }
  if (settings.fallbackCollectionDay != null && $("#fallback-collection-day")) {
    $("#fallback-collection-day").value = settings.fallbackCollectionDay;
  }
  $("#resume-run").checked = false;
  $("#dry-run").checked = false;
  if (settings.cpus != null) $("#cpus").value = settings.cpus;
  if (settings.memory != null) $("#memory").value = settings.memory;
  const sharedMinAbundancePct = settings.analysisMinAbundancePct
    ?? settings.dinemitesMinAbundancePct
    ?? settings.dciferMinAbundancePct;
  const sharedAbundanceDenominator = settings.analysisAbundanceDenominator
    || settings.dinemitesAbundanceDenominator
    || settings.dciferAbundanceDenominator;
  if (sharedMinAbundancePct != null && $("#analysis-min-abundance-pct")) {
    $("#analysis-min-abundance-pct").value = sharedMinAbundancePct;
  }
  if (sharedAbundanceDenominator && $("#analysis-abundance-denominator")) {
    $("#analysis-abundance-denominator").value = sharedAbundanceDenominator;
  }
  const restoredAnalysisMode = settings.analysisCdhitMode === "summed"
    ? "summed"
    : "off";
  setAnalysisMode(restoredAnalysisMode, false);
  if ($("#dinemites-enable")) $("#dinemites-enable").checked = true;
  if (settings.dinemitesModel && $("#dinemites-model")) $("#dinemites-model").value = settings.dinemitesModel;
  if (settings.dinemitesNLags != null && $("#dinemites-n-lags")) $("#dinemites-n-lags").value = settings.dinemitesNLags;
  if (settings.dinemitesTLag != null && $("#dinemites-t-lag")) $("#dinemites-t-lag").value = settings.dinemitesTLag;
  if (typeof settings.dinemitesNoDayCutoff === "boolean" && $("#dinemites-no-day-cutoff")) {
    $("#dinemites-no-day-cutoff").checked = settings.dinemitesNoDayCutoff;
  }
  if (settings.dinemitesSeed != null && $("#dinemites-seed")) $("#dinemites-seed").value = settings.dinemitesSeed;
  if (settings.dinemitesRefresh != null && $("#dinemites-refresh-interval")) $("#dinemites-refresh-interval").value = settings.dinemitesRefresh;
  if (settings.dinemitesBayesianLagDays != null && $("#dinemites-bayesian-lag-days")) $("#dinemites-bayesian-lag-days").value = settings.dinemitesBayesianLagDays;
  if (settings.dinemitesBayesianChains != null && $("#dinemites-bayesian-chains")) $("#dinemites-bayesian-chains").value = settings.dinemitesBayesianChains;
  if (settings.dinemitesBayesianParallelChains != null && $("#dinemites-bayesian-parallel-chains")) {
    $("#dinemites-bayesian-parallel-chains").value = settings.dinemitesBayesianParallelChains;
  }
  if (settings.dinemitesBayesianWarmup != null && $("#dinemites-bayesian-warmup")) $("#dinemites-bayesian-warmup").value = settings.dinemitesBayesianWarmup;
  if (settings.dinemitesBayesianSampling != null && $("#dinemites-bayesian-sampling")) $("#dinemites-bayesian-sampling").value = settings.dinemitesBayesianSampling;
  if (settings.dinemitesBayesianAdaptDelta != null && $("#dinemites-bayesian-adapt-delta")) {
    $("#dinemites-bayesian-adapt-delta").value = settings.dinemitesBayesianAdaptDelta;
  }
  const legacyCovariates = String(settings.dinemitesInfectionCovariates || "");
  if ($("#dinemites-covariate-season")) {
    $("#dinemites-covariate-season").checked = typeof settings.dinemitesUseSeasonCovariate === "boolean"
      ? settings.dinemitesUseSeasonCovariate
      : !legacyCovariates || legacyCovariates === "auto" || legacyCovariates.includes("covariate_season");
  }
  if ($("#dinemites-covariate-age")) {
    $("#dinemites-covariate-age").checked = typeof settings.dinemitesUseAgeCovariate === "boolean"
      ? settings.dinemitesUseAgeCovariate
      : legacyCovariates.includes("covariate_age");
  }
  if ($("#dinemites-covariate-gender")) {
    $("#dinemites-covariate-gender").checked = typeof settings.dinemitesUseGenderCovariate === "boolean"
      ? settings.dinemitesUseGenderCovariate
      : legacyCovariates.includes("covariate_gender");
  }
  if ($("#dinemites-infection-covariates")) {
    $("#dinemites-infection-covariates").value = settings.dinemitesCustomCovariates != null
      ? settings.dinemitesCustomCovariates
      : "";
  }
  if (typeof settings.dinemitesBayesianDropOut === "boolean" && $("#dinemites-bayesian-drop-out")) {
    $("#dinemites-bayesian-drop-out").checked = settings.dinemitesBayesianDropOut;
  }
  if ($("#dcifer-enable")) $("#dcifer-enable").checked = true;
  if (settings.dciferCoiLrank != null && $("#dcifer-coi-lrank")) $("#dcifer-coi-lrank").value = settings.dciferCoiLrank;
  if (settings.dciferIbdGridNr != null && $("#dcifer-ibd-grid-nr")) $("#dcifer-ibd-grid-nr").value = settings.dciferIbdGridNr;
  if (settings.dciferAlpha != null && $("#dcifer-alpha")) $("#dcifer-alpha").value = settings.dciferAlpha;
  syncGeneratedSampleSheetPath();
  normalizePathInputs();
}

async function fetchJson(url, options = {}, settings = {}) {
  const {allowAppError = false} = settings;
  const response = await fetch(url, {
    headers: {"Content-Type": "application/json"},
    ...options
  });
  const payload = await response.json().catch(() => ({}));
  if (!response.ok || (!allowAppError && payload.ok === false)) {
    throw new Error(payload.error || `Request failed: ${response.status}`);
  }
  return payload;
}

function postJson(url, body, settings = {}) {
  return fetchJson(url, {
    method: "POST",
    body: JSON.stringify(body)
  }, settings);
}

async function chooseNativeFolder({initial, prompt, allowNewFolder = false}) {
  return postJson("/api/select-folder", {
    initial,
    prompt,
    allow_new_folder: allowNewFolder
  });
}

async function chooseNativeFile({initial, prompt, kind = "metadata"}) {
  return postJson("/api/select-file", {initial, prompt, kind});
}

function withTimeout(promise, ms, message) {
  let timeoutId;
  const timeout = new Promise((_, reject) => {
    timeoutId = window.setTimeout(() => reject(new Error(message)), ms);
  });
  return Promise.race([promise, timeout]).finally(() => window.clearTimeout(timeoutId));
}

function delay(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

async function withButtonFeedback(button, busyText, action) {
  if (!button) return action();
  const oldText = button.textContent;
  const oldLabel = button.getAttribute("aria-label");
  const oldTitle = button.getAttribute("title");
  const iconOnly = button.classList.contains("icon-button");
  button.disabled = true;
  if (iconOnly) {
    button.classList.add("is-refreshing");
    button.setAttribute("aria-label", busyText);
    button.setAttribute("title", busyText);
  } else {
    text(button, busyText);
  }
  try {
    await Promise.all([action(), delay(900)]);
  } finally {
    if (iconOnly) {
      button.classList.remove("is-refreshing");
      if (oldLabel) button.setAttribute("aria-label", oldLabel);
      if (oldTitle) button.setAttribute("title", oldTitle);
    } else {
      text(button, oldText);
    }
    button.disabled = false;
  }
}

function bindPlotWheelScrolling() {
  document.addEventListener("wheel", (event) => {
    const target = event.target instanceof Element ? event.target : null;
    const plotCard = target?.closest(".dinemites-plot-card, .dcifer-heatmap-card");
    if (!plotCard || event.defaultPrevented || event.shiftKey) return;
    if (Math.abs(event.deltaY) <= Math.abs(event.deltaX)) return;

    const scrollRoot = document.scrollingElement || document.documentElement;
    const maxScrollTop = scrollRoot.scrollHeight - window.innerHeight;
    const scrollingDown = event.deltaY > 0;
    const canScrollPage = scrollingDown
      ? window.scrollY < maxScrollTop - 1
      : window.scrollY > 1;
    if (!canScrollPage) return;

    const deltaScale = event.deltaMode === 1
      ? 24
      : event.deltaMode === 2
        ? window.innerHeight
        : 1;
    event.preventDefault();
    window.scrollBy({top: event.deltaY * deltaScale, left: 0, behavior: "auto"});
  }, {passive: false});
}

function setPill(node, label, status) {
  if (!node) return;
  const value = String(label == null ? "" : label);
  const quietDefault = !status && ["pending", "Pending", "No scan yet", "No preview"].includes(value);
  node.classList.remove("status-pill", "status-meta", "ok", "warn", "bad");
  text(node, value);
  if (quietDefault) {
    node.hidden = true;
    node.classList.add("status-meta");
    return;
  }
  node.hidden = false;
  node.classList.add(status ? "status-pill" : "status-meta");
  if (status) node.classList.add(status);
}

function payloadStatus(payload) {
  const state = payload?.state || {};
  const summary = payload?.summary || {};
  const status = state.status || summary.status || "pending";
  if (payload?.active && status === "pending") return "starting";
  return status;
}

function isActiveStatus(status) {
  return status === "starting" || status === "running";
}

function checkStatusReady() {
  return ($("#check-status")?.textContent || "").trim().toLowerCase() === "ready";
}

function setScanReady(ready) {
  scanReady = Boolean(ready);
  const button = $("#proceed-run-button");
  if (button) {
    button.disabled = !scanReady;
    button.title = scanReady ? "Continue to run setup." : "Scan a FASTQ folder successfully first.";
  }
  updateRunButtonAvailability();
}

function runStartBlocker(active = isActiveStatus(lastRunStatus)) {
  if (active) return "A workflow is already running.";
  if (!scanReady) return "Return to Configuration and scan the FASTQ folder before starting.";
  const availableLibraries = libraryCheckboxes();
  if (availableLibraries.length && !checkedRunLibraries().length) {
    return "Select at least one sequencing library.";
  }
  const outputFolder = String($("#outdir")?.value || "").trim();
  if (!outputFolder) return "Choose an output folder before starting.";
  if (!isAbsolutePath(outputFolder)) return "Choose a full output folder path using Choose folder.";
  return "";
}

function updateRunButtonAvailability(active = isActiveStatus(lastRunStatus)) {
  const button = $("#run-button");
  if (!button) return;
  const blocker = runStartBlocker(active);
  button.disabled = Boolean(blocker);
  button.title = blocker || "Start the workflow.";

  if (!active) {
    const message = $("#run-readiness");
    if (message) {
      text(message, blocker || "Ready to start.");
      message.className = `field-help ${blocker ? "warn" : "ok"}`;
    }
  }
}

function updateScanButtonAvailability() {
  const button = $("#scan-button");
  if (!button) return;
  const fastqFolder = String($("#fastq-dir")?.value || "").trim();
  button.disabled = scanInFlight || !fastqFolder;
  button.title = fastqFolder ? "Scan the selected FASTQ folder." : "Choose a FASTQ folder first.";
}

function stageClass(status) {
  if (status === "complete") return "complete";
  if (status === "started" || status === "running") return "running";
  if (status === "failed" || status === "error") return "failed";
  return "pending";
}

function stageStatusLabel(status) {
  const klass = stageClass(status);
  if (klass === "running") return "running";
  if (klass === "failed") return "failed";
  return "";
}

function setProgressDisplay(value) {
  const percent = Math.max(0, Math.min(100, Math.round(Number(value) || 0)));
  displayedProgressPercent = percent;
  const progressFill = $("#progress-fill");
  const globalFill = $("#global-progress-fill");
  const progressBar = progressFill?.parentElement;
  if (progressFill) progressFill.style.width = `${percent}%`;
  if (globalFill) globalFill.style.width = `${percent}%`;
  if (progressBar) progressBar.setAttribute("aria-valuenow", String(percent));
  text($("#progress-percent"), `${percent}%`);
  text($("#pipeline-percent"), `${percent}%`);
}

function failureMessageFromLog(raw) {
  const value = String(raw || "");
  if (/Could not find or load main class/i.test(value)) {
    return "Nextflow could not start because a Windows path was parsed incorrectly.";
  }
  if (/unrecognized arguments:.*(?:Users|mnt\/c)/i.test(value)) {
    return "A Windows file path could not be parsed correctly by the workflow.";
  }
  if (/one or both of the fastq files not found/i.test(value)) {
    return "One or more FASTQ files could not be opened. Re-scan the source folder and try again.";
  }
  const processMatch = value.match(/Error executing process > '([^']+)'/i);
  if (processMatch) {
    return `${processMatch[1].replace(/\s*\([^)]*\)$/, "")} stopped with an error. Open the technical log for the command output.`;
  }
  return "The workflow stopped unexpectedly. Open the technical log for the command output.";
}

function updateRunFailureDetail(raw) {
  text($("#run-failure-detail"), failureMessageFromLog(raw));
}

function animateProgressTo(percent) {
  const target = Math.max(0, Math.min(100, Math.round(Number(percent) || 0)));
  if (target === displayedProgressPercent) return;
  setProgressDisplay(target);
}

function compactTerminalLog(raw, statusPayload) {
  const value = String(raw || "").replace(/\r\n/g, "\n").trimEnd();
  const status = payloadStatus(statusPayload);
  if (!value) {
    return isActiveStatus(status)
      ? "> Starting malaria-amplicon-nf.\n> Live Nextflow progress will appear here."
      : EMPTY_LOG_TEXT;
  }

  const lines = value.split("\n");
  const nextflowStart = lines.findIndex((line) => {
    const trimmed = line.trim();
    return (
      trimmed.includes("N E X T F L O W") ||
      trimmed.startsWith("WARN:") ||
      trimmed.startsWith("Launching `") ||
      trimmed.startsWith("executor >") ||
      /^\[[0-9a-f]{2}\//i.test(trimmed)
    );
  });
  if (nextflowStart >= 0) {
    return lines.slice(nextflowStart).slice(-220).join("\n").trimEnd();
  }

  const filtered = lines.filter((line) => {
    const trimmed = line.trim();
    if (!trimmed) return false;
    if (trimmed.startsWith("[SIMPLseq/App]")) return false;
    if (trimmed.startsWith("[OK]")) return false;
    if (trimmed.startsWith("[WARN] Large/high-depth datasets")) return false;
    if (trimmed.startsWith("[SIMPLseq] preparing")) return false;
    if (trimmed.startsWith("[SIMPLseq] wrote runtime versions")) return false;
    if (trimmed.startsWith("[SIMPLseq] wrote input FASTQ MD5 table")) return false;
    return true;
  });
  const cleaned = filtered.length ? filtered : lines.slice(-40);
  return cleaned.slice(-180).join("\n").trimEnd() || EMPTY_LOG_TEXT;
}

function setStep(name, active) {
  const node = document.querySelector(`.step[data-step="${name}"]`);
  if (node) node.classList.toggle("is-active", active);
}

function selectTab(requestedName) {
  const name = requestedName === "analysis" ? "qc" : requestedName;
  const downstreamTabs = new Set(["qc", "dinemites", "dcifer"]);
  const appShell = document.querySelector(".app-shell");
  if (appShell) appShell.classList.toggle("is-downstream-mode", downstreamTabs.has(name));
  const currentPanel = document.querySelector(".tab-panel.is-active");
  const nextPanelId = `tab-${name}`;
  const isPanelChange = !currentPanel || currentPanel.id !== nextPanelId;
  document.querySelectorAll(".tab").forEach((tab) => {
    const active = tab.dataset.tab === name || (tab.dataset.tab === "analysis" && downstreamTabs.has(name));
    tab.classList.toggle("is-active", active);
    if (active) tab.setAttribute("aria-current", "step");
    else tab.removeAttribute("aria-current");
  });
  document.querySelectorAll(".tab-panel").forEach((panel) => {
    const active = panel.id === nextPanelId;
    panel.classList.toggle("is-active", active);
    panel.classList.remove("is-entering");
    if (active && isPanelChange) {
      requestAnimationFrame(() => {
        panel.classList.add("is-entering");
        window.setTimeout(() => panel.classList.remove("is-entering"), 260);
      });
    }
  });
  const subtabs = $("#downstream-subtabs");
  if (subtabs) subtabs.hidden = !downstreamTabs.has(name);
  const downstreamOrder = ["qc", "dinemites", "dcifer"];
  const activeDownstreamIndex = downstreamOrder.indexOf(name);
  document.querySelectorAll(".downstream-subtab").forEach((tab) => {
    const tabIndex = downstreamOrder.indexOf(tab.dataset.downstreamTab);
    const isActive = tab.dataset.downstreamTab === name;
    tab.classList.toggle("is-active", isActive);
    tab.classList.toggle("is-past", activeDownstreamIndex > -1 && tabIndex < activeDownstreamIndex);
    if (isActive) tab.setAttribute("aria-current", "step");
    else tab.removeAttribute("aria-current");
  });
  if (name === "results" || name === "qc") {
    refreshResults().catch(() => {
      updateResultsRunSummary({}, activeOutdir || $("#results-outdir")?.value || "");
    });
  } else if (name === "dinemites") {
    loadDinemitesResults();
  } else if (name === "dcifer") {
    loadDciferResults();
  }
}

function setFolderMessage(message, className = "") {
  const node = $("#folder-message");
  if (!node) return;
  node.className = `field-help ${className}`.trim();
  text(node, message);
}

function setMetadataMessage(message, className = "") {
  const node = $("#metadata-path-help");
  if (!node) return;
  node.className = `field-help ${className}`.trim();
  text(node, message);
}

function resetMetadataContract() {
  metadataContract = {
    schema_version: 1,
    columns: {},
    detection_value_map: {},
    excluded_status_values: []
  };
  latestMetadataInspection = null;
  metadataContractValid = true;
}

function metadataContractPayload() {
  return {
    schema_version: 1,
    columns: {...(metadataContract.columns || {})},
    detection_value_map: {...(metadataContract.detection_value_map || {})},
    excluded_status_values: [...(metadataContract.excluded_status_values || [])]
  };
}

function populateMetadataSheets(payload) {
  const select = $("#metadata-sheet");
  if (!select) return;
  const sheets = Array.isArray(payload.sheets) ? payload.sheets : [];
  const selected = String(payload.sheet || "");
  select.replaceChildren();
  if (!sheets.length) {
    const option = document.createElement("option");
    option.value = "";
    text(option, "Not applicable");
    select.appendChild(option);
    select.disabled = true;
    return;
  }
  sheets.forEach((sheet) => {
    const option = document.createElement("option");
    option.value = sheet;
    text(option, sheet);
    option.selected = sheet === selected;
    select.appendChild(option);
  });
  select.disabled = sheets.length < 2;
}

function populateMetadataColumnSelectors(payload) {
  const columns = Array.isArray(payload.available_columns) ? payload.available_columns : [];
  document.querySelectorAll("[data-metadata-column]").forEach((select) => {
    const target = select.dataset.metadataColumn;
    const selected = String(payload.columns?.[target] || "");
    select.replaceChildren();
    const blank = document.createElement("option");
    blank.value = "";
    text(blank, "Not provided");
    select.appendChild(blank);
    columns.forEach((column) => {
      const option = document.createElement("option");
      option.value = column;
      text(option, column);
      option.selected = column === selected;
      select.appendChild(option);
    });
    select.value = selected;
  });
}

function renderMetadataDetectionValues(payload) {
  const section = $("#metadata-detection-contract");
  const container = $("#metadata-detection-values");
  if (!section || !container) return;
  const pcrColumn = String(payload.columns?.metadata_pcr || "");
  const values = Array.isArray(payload.detection_values) ? payload.detection_values : [];
  section.hidden = !pcrColumn;
  container.replaceChildren();
  if (!pcrColumn) return;
  if (!values.length) {
    const empty = document.createElement("p");
    empty.className = "field-help";
    text(empty, "The selected detection column contains no non-empty values.");
    container.appendChild(empty);
    return;
  }
  values.forEach((item) => {
    const row = document.createElement("div");
    row.className = "metadata-detection-row";
    const value = document.createElement("div");
    value.className = "metadata-detection-row__value";
    const label = document.createElement("strong");
    text(label, item.value);
    const count = document.createElement("span");
    text(count, `${Number(item.count || 0).toLocaleString()} row${Number(item.count || 0) === 1 ? "" : "s"}`);
    value.append(label, count);
    const select = document.createElement("select");
    select.setAttribute("aria-label", `Meaning of ${item.value}`);
    [
      ["positive", "PCR positive"],
      ["negative", "PCR negative"],
      ["ignore", "Ignore value"],
      ["review", "Needs review"]
    ].forEach(([state, labelText]) => {
      const option = document.createElement("option");
      option.value = state;
      text(option, labelText);
      select.appendChild(option);
    });
    select.value = item.state || "review";
    metadataContract.detection_value_map[item.value] = select.value;
    select.addEventListener("change", () => {
      metadataContract.detection_value_map[item.value] = select.value;
      invalidateScanReady();
      renderMetadataCalendarPreview(latestMetadataInspection);
    });
    row.append(value, select);
    container.appendChild(row);
  });
}

function renderMetadataStatusValues(payload) {
  const section = $("#metadata-status-contract");
  const container = $("#metadata-status-values");
  if (!section || !container) return;
  const statusColumn = String(payload.columns?.metadata_status || "");
  const counts = payload.value_counts?.metadata_status || {};
  const values = Object.entries(counts).sort((left, right) => {
    return Number(right[1] || 0) - Number(left[1] || 0)
      || String(left[0]).localeCompare(String(right[0]));
  });
  section.hidden = !statusColumn || values.length === 0;
  container.replaceChildren();
  if (section.hidden) return;
  const excluded = new Set(metadataContract.excluded_status_values || []);
  values.forEach(([rawValue, count]) => {
    const label = document.createElement("label");
    label.className = "metadata-status-row";
    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = excluded.has(rawValue);
    const value = document.createElement("span");
    const strong = document.createElement("strong");
    const detail = document.createElement("small");
    text(strong, rawValue);
    text(detail, `${Number(count || 0).toLocaleString()} row${Number(count || 0) === 1 ? "" : "s"}`);
    value.append(strong, detail);
    const action = document.createElement("em");
    text(action, "Exclude from calendar");
    checkbox.addEventListener("change", () => {
      const next = new Set(metadataContract.excluded_status_values || []);
      if (checkbox.checked) next.add(rawValue);
      else next.delete(rawValue);
      metadataContract.excluded_status_values = [...next].sort();
      invalidateScanReady();
    });
    label.append(checkbox, value, action);
    container.appendChild(label);
  });
}

function renderMetadataCalendarPreview(payload) {
  const container = $("#metadata-calendar-preview");
  const status = $("#metadata-contract-status");
  const title = $("#metadata-contract-title");
  if (!container || !status || !payload) return;
  const columns = metadataContract.columns || {};
  const hasIdentity = Boolean(columns.participant_id);
  const hasVisit = Boolean(columns.collection_date || columns.month);
  const hasDetection = Boolean(columns.metadata_pcr);
  metadataContractValid = hasIdentity && hasVisit && payload.ok !== false;
  const counts = {positive: 0, negative: 0, ignore: 0, review: 0};
  (payload.detection_values || []).forEach((item) => {
    const state = metadataContract.detection_value_map[item.value] || item.state || "review";
    counts[state] = (counts[state] || 0) + Number(item.count || 0);
  });
  container.replaceChildren();
  const summary = document.createElement("div");
  summary.className = "metadata-calendar-preview__summary";
  const headline = document.createElement("strong");
  const detail = document.createElement("span");
  if (!metadataContractValid) {
    text(title, "Mapping needed");
    text(headline, "Metadata mapping incomplete");
    text(detail, "Choose a participant column and a complete date or visit-month column.");
    status.className = "metadata-contract__status is-review";
    text(status, "Needs mapping");
  } else if (!hasDetection) {
    text(title, "Dates mapped");
    text(headline, `${Number(payload.records || 0).toLocaleString()} dated metadata rows`);
    text(detail, "Sequencing-only longitudinal calendar: ungenotyped visits cannot be classified as negative or missing.");
    status.className = "metadata-contract__status is-limited";
    text(status, "Calendar limited");
  } else {
    text(title, "Spreadsheet ready");
    text(headline, `${Number(payload.records || 0).toLocaleString()} dated metadata rows`);
    text(detail, `${counts.positive.toLocaleString()} positive | ${counts.negative.toLocaleString()} negative | ${counts.review.toLocaleString()} review | ${counts.ignore.toLocaleString()} ignored`);
    status.className = `metadata-contract__status${counts.review ? " is-limited" : " is-ready"}`;
    text(status, counts.review ? "Ready with review" : "Ready");
  }
  summary.append(headline, detail);
  container.appendChild(summary);
}

async function inspectMetadata() {
  const path = $("#metadata-path")?.value.trim();
  const container = $("#metadata-inspection");
  const contractSection = $("#metadata-contract");
  if (!container || !contractSection) return;
  if (!path) {
    contractSection.hidden = true;
    container.hidden = true;
    container.replaceChildren();
    latestMetadataInspection = null;
    metadataContractValid = true;
    return;
  }
  contractSection.hidden = false;
  container.hidden = false;
  container.className = "metadata-inspection is-loading";
  text(container, "Reading metadata structure...");
  try {
    const payload = await postJson("/api/metadata/inspect", {
      path,
      sheet: $("#metadata-sheet")?.value.trim() || "",
      date_order: $("#metadata-date-order")?.value || "auto",
      metadata_contract: metadataContractPayload()
    });
    latestMetadataInspection = payload;
    metadataContract = payload.metadata_contract || metadataContractPayload();
    populateMetadataSheets(payload);
    populateMetadataColumnSelectors(payload);
    renderMetadataDetectionValues(payload);
    renderMetadataStatusValues(payload);
    renderMetadataCalendarPreview(payload);
    const reviewDetails = $("#metadata-review-details");
    if (reviewDetails && !metadataContractValid) reviewDetails.open = true;
    container.className = "metadata-inspection";
    container.replaceChildren();
    const warnings = (payload.issues || []).filter((issue) => issue.severity !== "info");
    warnings.slice(0, 3).forEach((issue) => {
      const item = document.createElement("p");
      item.className = `metadata-inspection__issue is-${issue.severity}`;
      text(item, issue.message);
      container.append(item);
    });
    container.hidden = warnings.length === 0;
    const hasDetection = Boolean(metadataContract.columns?.metadata_pcr);
    setMetadataMessage(
      metadataContractValid
        ? hasDetection
          ? "Metadata mapped. Confirm detection values before scanning."
          : "Metadata mapped for dates. Add PCR/qPCR results for a complete visit calendar."
        : "Complete the metadata mapping before scanning.",
      metadataContractValid && hasDetection ? "ok" : "warn"
    );
  } catch (error) {
    latestMetadataInspection = null;
    metadataContractValid = false;
    contractSection.hidden = false;
    const reviewDetails = $("#metadata-review-details");
    if (reviewDetails) reviewDetails.open = true;
    container.hidden = false;
    container.className = "metadata-inspection has-error";
    text(container, error.message || "Metadata could not be inspected.");
    setMetadataMessage("Metadata needs attention before scanning.", "warn");
  }
}

function invalidateScanReady() {
  setScanReady(false);
}

function setSampleValidationVisible(visible) {
  const panel = $("#sample-validation-panel");
  if (panel) panel.hidden = !visible;
}

function row(cells) {
  const tr = document.createElement("tr");
  cells.forEach((value) => {
    const td = document.createElement("td");
    if (value instanceof Node) {
      td.appendChild(value);
    } else {
      text(td, value);
    }
    tr.appendChild(td);
  });
  return tr;
}

function boolish(value) {
  return value === true || String(value || "").toLowerCase() === "true";
}

function dateWithInferredParts(value, inferredYear, inferredDay) {
  const raw = String(value || "");
  const match = raw.match(/^(\d{4})-(\d{2})-(\d{2})$/);
  if (!match) return raw;
  const wrap = document.createElement("span");
  const year = document.createElement(boolish(inferredYear) ? "mark" : "span");
  const day = document.createElement(boolish(inferredDay) ? "mark" : "span");
  if (boolish(inferredYear)) year.className = "date-inferred";
  if (boolish(inferredDay)) day.className = "date-inferred";
  year.textContent = match[1];
  day.textContent = match[3];
  wrap.className = "date-value";
  wrap.append(year, document.createTextNode(`-${match[2]}-`), day);
  if (boolish(inferredYear) || boolish(inferredDay)) {
    const source = document.createElement("small");
    source.className = "date-source";
    source.textContent = "inferred";
    wrap.appendChild(source);
  }
  return wrap;
}

function sampleStatus(item) {
  const sampleType = String(item.sample_type || "").toLowerCase();
  const metadataStatus = String(item.metadata_match_status || "").toLowerCase();
  const validationStatus = String(item.validation_status || item.status || "").toLowerCase();
  const reviewStatus = String(item.review_status || "").toLowerCase();
  if (reviewStatus === "excluded") return {label: "Excluded", tone: "neutral", key: "excluded"};
  if (reviewStatus === "review") return {label: "Review", tone: "warning", key: "review"};
  if (reviewStatus === "ready") return {label: "Ready", tone: "success", key: "ready"};
  if (["control", "negative", "positive"].includes(sampleType)) {
    return {label: "Excluded", tone: "neutral", key: "excluded"};
  }
  if (item.validation_error || validationStatus === "failed" || validationStatus === "error") {
    return {label: "Failed", tone: "danger", key: "failed"};
  }
  if (
    boolish(item.inferred_year)
    || boolish(item.inferred_day)
    || metadataStatus === "ambiguous"
    || !item.collection_date
  ) {
    return {label: "Review", tone: "warning", key: "review"};
  }
  return {label: "Ready", tone: "success", key: "ready"};
}

function sampleStatusBadge(item) {
  const status = sampleStatus(item);
  const badge = document.createElement("span");
  badge.className = `sample-status is-${status.tone}`;
  text(badge, status.label);
  return badge;
}

function sampleAssessment(item) {
  const assessment = document.createElement("span");
  assessment.className = "sample-assessment";
  assessment.appendChild(sampleStatusBadge(item));
  const reason = String(item.review_reason || "").trim();
  if (reason) {
    const detail = document.createElement("small");
    text(detail, reason);
    assessment.appendChild(detail);
  }
  return assessment;
}

function emptyRow(colspan, message) {
  const tr = document.createElement("tr");
  const td = document.createElement("td");
  td.colSpan = colspan;
  td.className = "empty";
  text(td, message);
  tr.appendChild(td);
  return tr;
}

function displayMissing(value) {
  return value === null || value === undefined || value === "" ? "--" : value;
}

function formatNumber(value, digits = 2) {
  if (value === null || value === undefined || value === "") return "--";
  const number = Number(value);
  if (!Number.isFinite(number)) return String(value);
  if (Number.isInteger(number)) return String(number);
  return number.toFixed(digits).replace(/\.?0+$/, "");
}

function formatPValue(value) {
  if (value === null || value === undefined || value === "") return "--";
  const number = Number(value);
  if (!Number.isFinite(number)) return String(value);
  if (number > 0 && number < 0.001) return number.toExponential(2);
  return formatNumber(number, 3);
}

function updatePlotJump(gallery, countNode, button, count, emptyText, singularText, pluralText) {
  if (countNode) {
    text(countNode, count ? `${count} ${count === 1 ? singularText : pluralText} available` : emptyText);
  }
  if (!button) return;
  button.disabled = !count;
  button.onclick = count && gallery
    ? () => gallery.scrollIntoView({ behavior: "smooth", block: "start" })
    : null;
}

function filteredScanReviewRows() {
  const query = String($("#scan-search")?.value || "").trim().toLowerCase();
  return scanReviewRows.filter((item) => {
    const status = sampleStatus(item).key;
    if (scanReviewFilter !== "all" && status !== scanReviewFilter) return false;
    if (scanReviewLibrary !== "all" && String(item.library || "") !== scanReviewLibrary) return false;
    if (!query) return true;
    return [item.sample_id, item.biological_sample_id, item.participant_id, item.library, item.collection_date, item.replicate, item.review_reason]
      .some((value) => String(value || "").toLowerCase().includes(query));
  });
}

function scanRowsInCurrentContext() {
  const query = String($("#scan-search")?.value || "").trim().toLowerCase();
  return scanReviewRows.filter((item) => {
    if (scanReviewLibrary !== "all" && String(item.library || "") !== scanReviewLibrary) return false;
    if (!query) return true;
    return [item.sample_id, item.biological_sample_id, item.participant_id, item.library, item.collection_date, item.replicate, item.review_reason]
      .some((value) => String(value || "").toLowerCase().includes(query));
  });
}

function libraryCheckboxes() {
  return [...document.querySelectorAll('#run-library-options input[data-library]:not([data-library="all"])')];
}

function checkedRunLibraries() {
  return libraryCheckboxes()
    .filter((input) => input.checked)
    .map((input) => input.dataset.library);
}

function updateLibraryPickerSummary() {
  const checkboxes = libraryCheckboxes();
  const selected = checkedRunLibraries();
  const allCheckbox = $('#run-library-options input[data-library="all"]');
  if (allCheckbox) {
    allCheckbox.checked = checkboxes.length > 0 && selected.length === checkboxes.length;
    allCheckbox.indeterminate = selected.length > 0 && selected.length < checkboxes.length;
  }

  const selectedPairs = selected.reduce((sum, library) => sum + Number(detectedLibraryCounts[library] || 0), 0);
  const summary = $("#run-library-summary");
  const picker = $("#run-library-picker");
  if (!checkboxes.length) {
    text(summary, "Scan FASTQs to detect libraries");
    picker?.classList.remove("has-error");
    return;
  }
  if (!selected.length) {
    text(summary, "Select at least one library");
    picker?.classList.add("has-error");
    return;
  }
  picker?.classList.remove("has-error");
  if (selected.length === checkboxes.length) {
    text(summary, `All libraries (${selectedPairs} pairs)`);
  } else if (selected.length <= 2) {
    text(summary, `${selected.join(" + ")} (${selectedPairs} pairs)`);
  } else {
    text(summary, `${selected.length} libraries (${selectedPairs} pairs)`);
  }
}

function renderRunLibraryOptions(entries, total) {
  const options = $("#run-library-options");
  const trigger = $("#run-library-trigger");
  const picker = $("#run-library-picker");
  if (!options || !trigger || !picker) return;

  detectedLibraryCounts = Object.fromEntries(entries);
  options.replaceChildren();
  const restored = Array.isArray(restoredRunLibraries)
    ? new Set(restoredRunLibraries.filter((library) => detectedLibraryCounts[library]))
    : null;
  const selected = restored && restored.size ? restored : new Set(entries.map(([library]) => library));

  const addOption = (library, label, count, checked, isAll = false) => {
    const row = document.createElement("label");
    row.className = `library-option${isAll ? " library-option--all" : ""}`;
    const input = document.createElement("input");
    input.type = "checkbox";
    input.dataset.library = library;
    input.checked = checked;
    const name = document.createElement("span");
    name.className = "library-option__name";
    text(name, label);
    const countNode = document.createElement("span");
    countNode.className = "library-option__count";
    text(countNode, `${count} pairs`);
    row.append(input, name, countNode);
    options.appendChild(row);
  };

  addOption("all", "All libraries", total, selected.size === entries.length, true);
  entries.forEach(([library, count]) => addOption(library, library, count, selected.has(library)));
  restoredRunLibraries = null;
  trigger.disabled = entries.length === 0;
  picker.classList.toggle("is-disabled", entries.length === 0);
  updateLibraryPickerSummary();
}

function renderLibrarySelectors(counts = {}, total = 0) {
  const entries = Object.entries(counts)
    .filter(([library, count]) => library && Number(count) > 0)
    .sort(([left], [right]) => left.localeCompare(right, undefined, {numeric: true}));
  const scanSelect = $("#scan-library-filter");
  const populateScanFilter = (select, allLabel) => {
    if (!select) return;
    select.replaceChildren();
    const allOption = document.createElement("option");
    allOption.value = "all";
    allOption.textContent = allLabel;
    select.appendChild(allOption);
    entries.forEach(([library, count]) => {
      const option = document.createElement("option");
      option.value = library;
      option.textContent = `${library} (${count} pairs)`;
      select.appendChild(option);
    });
    select.value = "all";
    select.disabled = entries.length < 2;
  };
  populateScanFilter(scanSelect, entries.length ? `All libraries (${total})` : `All detected data (${total})`);
  renderRunLibraryOptions(entries, total);
  scanReviewLibrary = "all";
  if (scanSelect) scanSelect.value = "all";
  const help = entries.length
    ? `${entries.length} sequencing libraries detected. Open the selector to include any combination.`
    : "No library label was found in the FASTQ names, so the run will include all detected data.";
  text($("#run-library-help"), help);
}

function renderScanPreview(items = null) {
  if (Array.isArray(items)) {
    scanReviewRows = items;
    scanReviewFilter = "all";
    scanReviewLibrary = "all";
    scanReviewPage = 1;
    if ($("#scan-search")) $("#scan-search").value = "";
  }
  const tbody = $("#scan-preview");
  const tableWrap = document.querySelector(".scan-preview-wrap");
  if (tableWrap) tableWrap.scrollTop = 0;
  tbody.replaceChildren();
  const filtered = filteredScanReviewRows();
  renderSampleStatusSummary();
  if (!scanReviewRows.length) {
    tbody.appendChild(emptyRow(6, "Scan a folder to review parsed samples."));
    $("#scan-pagination").hidden = true;
    return;
  }
  if (!filtered.length) {
    tbody.appendChild(emptyRow(6, "No samples match this filter."));
  }
  const totalPages = Math.max(1, Math.ceil(filtered.length / SCAN_PAGE_SIZE));
  scanReviewPage = Math.min(Math.max(scanReviewPage, 1), totalPages);
  const start = (scanReviewPage - 1) * SCAN_PAGE_SIZE;
  const pageRows = filtered.slice(start, start + SCAN_PAGE_SIZE);
  pageRows.forEach((item) => {
    tbody.appendChild(row([
      item.biological_sample_id || item.sample_id,
      item.participant_id,
      item.library || "--",
      dateWithInferredParts(item.collection_date, item.inferred_year, item.inferred_day),
      item.replicate,
      sampleAssessment(item)
    ]));
  });
  const end = Math.min(start + SCAN_PAGE_SIZE, filtered.length);
  text($("#scan-page-status"), filtered.length ? `${start + 1}-${end} of ${filtered.length} samples` : "0 samples");
  text($("#scan-page-number"), `Page ${scanReviewPage} of ${totalPages}`);
  $("#scan-page-prev").disabled = scanReviewPage <= 1;
  $("#scan-page-next").disabled = scanReviewPage >= totalPages;
  $("#scan-pagination").hidden = false;
}

function renderSamplePreview(items) {
  const tbody = $("#sample-preview");
  tbody.replaceChildren();
  if (!items.length) {
    tbody.appendChild(emptyRow(6, "No sample rows available."));
    return;
  }
  items.slice(0, SAMPLE_PREVIEW_LIMIT).forEach((item) => {
    tbody.appendChild(row([
      item.sample_id,
      item.fastq_1,
      item.fastq_2,
      item.participant_id,
      dateWithInferredParts(item.collection_date, item.inferred_year, item.inferred_day),
      sampleStatusBadge(item)
    ]));
  });
  if (items.length > SAMPLE_PREVIEW_LIMIT) {
    const summary = document.createElement("tr");
    summary.className = "preview-limit-row";
    const cell = document.createElement("td");
    cell.colSpan = 6;
    text(cell, `Showing ${SAMPLE_PREVIEW_LIMIT} of ${items.length} generated rows.`);
    summary.appendChild(cell);
    tbody.appendChild(summary);
  }
}

function renderWarnings(scan) {
  const box = $("#scan-warnings");
  const readiness = $("#scan-readiness");
  const details = $("#scan-review-details");
  const reviewCount = $("#scan-review-count");
  box.replaceChildren();
  const notices = [];
  const pairCount = Number(scan.pair_count || 0);
  const totalBytes = Number(scan.total_fastq_bytes || 0);
  const pairingComplete = pairCount > 0 && !(scan.missing_r2 || []).length && !(scan.orphan_r2 || []).length;
  const identifiersUnique = pairCount > 0 && !(scan.duplicate_sample_ids || []).length;
  const resolvedCount = Number(scan.auto_disambiguated_sample_ids || 0);
  if (pairCount > 0) {
    readiness.hidden = false;
    readiness.className = `scan-readiness ${pairingComplete && identifiersUnique ? "is-ready" : "is-review"}`;
    const resolved = resolvedCount > 0 ? ` ${resolvedCount} ID collisions were resolved automatically.` : "";
    text(
      readiness,
      pairingComplete && identifiersUnique
        ? `Pairing complete. ${pairCount} libraries have unique sample IDs.${resolved}`
        : `${pairCount} libraries were found, but pairing or sample IDs need review.`
    );
  } else {
    readiness.hidden = true;
    text(readiness, "");
  }
  if (pairCount >= 50 || totalBytes >= 5 * 1024 * 1024 * 1024) {
    notices.push({
      className: "notice",
      text: "Large dataset detected. Full runs can require much more memory than a typical laptop. Run a small test first if this is a new setup."
    });
  }
  if (Number(scan.inferred_collection_dates || scan.collection_month_without_year || 0) > 0) {
    const yearCount = Number(scan.inferred_year_count || scan.collection_month_without_year || 0);
    const dayCount = Number(scan.inferred_day_count || scan.inferred_collection_dates || 0);
    const fallbackYear = String(scan.fallback_collection_year || "").trim();
    const fallbackDay = String(scan.fallback_collection_day || "27").replace(/^0+/, "") || "27";
    const parts = [];
    if (yearCount > 0) {
      const source = scan.fallback_year_autodetected ? " from metadata" : "";
      parts.push(`${yearCount} samples use inferred year ${fallbackYear || "from settings"}${source}`);
    }
    if (dayCount > 0) parts.push(`${dayCount} samples use inferred day ${fallbackDay}`);
    notices.push({
      className: "notice",
      text: `${parts.join(" and ")}. Review red date parts before running longitudinal analyses.`
    });
  }
  if (Number(scan.unresolved_collection_date_count || 0) > 0) {
    notices.push({
      className: "notice",
      text: `${scan.unresolved_collection_date_count} sample rows still have no collection date. Review unmatched metadata or add a fallback year before longitudinal analysis.`
    });
  }
  if (scan.metadata_path && scan.metadata_match_counts) {
    const matched = Number(scan.metadata_match_counts.matched || 0);
    const missing = Number(scan.metadata_match_counts.missing || 0);
    const ambiguous = Number(scan.metadata_match_counts.ambiguous || 0);
    notices.push({
      className: missing > 0 || ambiguous > 0 ? "notice" : "notice ok",
      text: `${scan.metadata_autodetected ? "Metadata detected automatically. " : ""}Matched ${matched} sample rows${missing > 0 ? `; ${missing} unmatched` : ""}${ambiguous > 0 ? `; ${ambiguous} ambiguous` : ""}.`
    });
  }
  if (scan.duplicate_sample_ids && scan.duplicate_sample_ids.length) {
    notices.push({
      className: "notice bad",
      text: `Duplicate sample IDs could not be resolved automatically: ${scan.duplicate_sample_ids.slice(0, 8).join(", ")}. Inspect duplicate or identically named FASTQ files.`
    });
  }
  if (scan.missing_r2 && scan.missing_r2.length) {
    notices.push({className: "notice", text: `Missing R2 files for ${scan.missing_r2.length} R1 files.`});
  }
  if (scan.orphan_r2 && scan.orphan_r2.length) {
    notices.push({className: "notice", text: `Found ${scan.orphan_r2.length} R2 files without R1.`});
  }
  notices.forEach((notice) => {
    const div = document.createElement("div");
    div.className = notice.className;
    text(div, notice.text);
    box.appendChild(div);
  });
  const hasBlockingIssues = Boolean(
    (scan.duplicate_sample_ids || []).length
    || (scan.missing_r2 || []).length
    || (scan.orphan_r2 || []).length
  );
  details.hidden = notices.length === 0;
  details.open = hasBlockingIssues;
  text(reviewCount, notices.length === 1 ? "1 note" : `${notices.length} notes`);
}

function renderSampleStatusSummary() {
  const summary = $("#sample-status-summary");
  const tools = $("#scan-table-tools");
  const scopedRows = scanRowsInCurrentContext();
  const ready = scopedRows.filter((item) => sampleStatus(item).key === "ready").length;
  const review = scopedRows.filter((item) => sampleStatus(item).key === "review").length;
  const excluded = scopedRows.filter((item) => sampleStatus(item).key === "excluded").length;
  const all = ready + review + excluded;
  text($("#status-all-count"), all);
  text($("#status-ready-count"), ready);
  text($("#status-review-count"), review);
  text($("#status-excluded-count"), excluded);
  tools.hidden = all === 0;
  summary.hidden = all === 0;
  document.querySelectorAll("[data-scan-filter]").forEach((button) => {
    const active = button.dataset.scanFilter === scanReviewFilter;
    button.classList.toggle("is-active", active);
    button.setAttribute("aria-pressed", String(active));
  });
}

async function scanFastqs() {
  if (scanInFlight) return;
  if ($("#metadata-path")?.value.trim() && !metadataContractValid) {
    const details = $(".metadata-workspace");
    if (details) details.open = true;
    const reviewDetails = $("#metadata-review-details");
    if (reviewDetails) reviewDetails.open = true;
    setMetadataMessage("Map the participant and visit columns before scanning.", "warn");
    $("#metadata-contract")?.scrollIntoView({behavior: "smooth", block: "center"});
    return;
  }
  scanInFlight = true;
  setScanReady(false);
  setSampleValidationVisible(false);
  saveSettings();
  const scanButton = $("#scan-button");
  const scanButtonLabel = scanButton.textContent;
  let scanStageIndex = 0;
  setPill($("#scan-status"), SCAN_STAGE_LABELS[scanStageIndex], "warn");
  scanButton.disabled = true;
  scanButton.classList.add("is-busy");
  text(scanButton, "Scanning...");
  const scanStageTimer = setInterval(() => {
    scanStageIndex = Math.min(scanStageIndex + 1, SCAN_STAGE_LABELS.length - 1);
    setPill($("#scan-status"), SCAN_STAGE_LABELS[scanStageIndex], "warn");
  }, 1400);
  try {
    const payload = await postJson("/api/scan", {
      fastq_dir: $("#fastq-dir").value,
      metadata_path: $("#metadata-path") ? $("#metadata-path").value : "",
      metadata_sheet: $("#metadata-sheet") ? $("#metadata-sheet").value : "",
      metadata_date_order: $("#metadata-date-order")?.value || "auto",
      metadata_contract: metadataContractPayload(),
      fallback_collection_year: $("#fallback-collection-year") ? $("#fallback-collection-year").value : "",
      fallback_collection_day: $("#fallback-collection-day") ? $("#fallback-collection-day").value : "27",
      absolute_paths: true,
      write_samples: false
    });
    if (payload.metadata_autodetected && payload.metadata_path && $("#metadata-path")) {
      $("#metadata-path").value = payload.metadata_path;
      const metadataName = payload.metadata_path.split(/[\\/]/).pop();
      setMetadataMessage(`Detected ${metadataName} beside the FASTQ files.`, "ok");
      resetMetadataContract();
      await inspectMetadata();
    }
    if (payload.fallback_year_autodetected && payload.fallback_collection_year && $("#fallback-collection-year")) {
      $("#fallback-collection-year").value = payload.fallback_collection_year;
    }
    text($("#metric-pairs"), payload.pair_count);
    text($("#metric-participants"), payload.participant_count || 0);
    text($("#metric-size"), payload.total_fastq_size);
    text($("#metric-missing"), payload.missing_pairs);
    renderLibrarySelectors(payload.library_counts || {}, Number(payload.pair_count || 0));
    text($("#sample-sheet-title"), "Generated sample records");
    renderScanPreview(payload.scan_rows || payload.sample_preview || payload.preview || []);
    renderSamplePreview(payload.sample_preview || []);
    setSampleValidationVisible(Boolean((payload.sample_preview || []).length));
    renderWarnings(payload);
    renderSampleStatusSummary();
    activeOutdir = "";
    syncGeneratedSampleSheetPath();
    updateSamplePathHelp();
    const hasDuplicates = Boolean(payload.duplicate_sample_ids && payload.duplicate_sample_ids.length);
    const successfulScan = Boolean(
      !hasDuplicates
      && Number(payload.pair_count || 0) > 0
      && (!payload.metadata_path || metadataContractValid)
    );
    if (successfulScan) {
      setPill($("#scan-status"), "No scan yet");
      setPill($("#sample-sheet-status"), "No preview");
      setStep("review", true);
      setScanReady(successfulScan);
      window.requestAnimationFrame(() => {
        const reviewSection = $("#scan-review-section");
        reviewSection?.scrollIntoView({behavior: "smooth", block: "start"});
        reviewSection?.focus({preventScroll: true});
      });
    } else if (hasDuplicates) {
      setPill($("#scan-status"), "Duplicates", "bad");
      setPill($("#sample-sheet-status"), "Preview only", "bad");
      setScanReady(false);
    } else {
      setPill($("#scan-status"), "No pairs", "warn");
      setPill($("#sample-sheet-status"), "No rows", "warn");
      setScanReady(false);
    }
    saveSettings();
  } catch (error) {
    setScanReady(false);
    setSampleValidationVisible(false);
    renderScanPreview([]);
    setPill($("#scan-status"), "Scan failed", "bad");
    renderWarnings({duplicate_sample_ids: [error.message]});
    renderSampleStatusSummary();
  } finally {
    clearInterval(scanStageTimer);
    scanInFlight = false;
    scanButton.classList.remove("is-busy");
    text(scanButton, scanButtonLabel || "Scan folder");
    updateScanButtonAvailability();
  }
}

function renderChecks(checks) {
  const tbody = $("#check-table");
  if (!tbody) return;
  tbody.replaceChildren();
  if (!checks.length) {
    tbody.appendChild(emptyRow(3, "No checks returned."));
    return;
  }
  checks.forEach((item) => {
    const tr = document.createElement("tr");
    const name = document.createElement("td");
    const status = document.createElement("td");
    const detail = document.createElement("td");
    const badge = document.createElement("span");
    badge.className = `check-badge ${item.status === "ok" ? "ok" : item.status === "warn" ? "warn" : "bad"}`;
    text(name, item.name);
    text(badge, item.status);
    text(detail, item.detail);
    status.appendChild(badge);
    tr.append(name, status, detail);
    tbody.appendChild(tr);
  });
}

async function runCheck() {
  saveSettings();
  const statusPill = $("#check-status");
  const checkButton = $("#check-button");
  setPill(statusPill, "Checking", "warn");
  if (checkButton) checkButton.disabled = true;
  try {
    const payload = await postJson("/api/check", {outdir: $("#outdir").value}, {allowAppError: true});
    renderChecks(payload.checks || []);
    await loadLog({silent: true});
    if (payload.failed) {
      setPill(statusPill, `${payload.failed} need attention`, "bad");
    } else {
      setPill(statusPill, "Ready", "ok");
      setStep("runtime", true);
      resetRunDisplay();
    }
  } catch (error) {
    setPill(statusPill, "Check failed", "bad");
    renderChecks([{name: "Runtime check", status: "missing", detail: error.message}]);
  } finally {
    if (checkButton) checkButton.disabled = false;
  }
}

function runPayload() {
  const availableLibraries = libraryCheckboxes();
  const selectedLibraries = checkedRunLibraries();
  if (availableLibraries.length && !selectedLibraries.length) {
    throw new Error("Select at least one sequencing library.");
  }
  const fastqDir = String($("#fastq-dir").value || "").trim();
  const outdir = String($("#outdir").value || "").trim();
  if (!scanReady || !fastqDir) throw new Error("Scan a FASTQ folder successfully before starting the run.");
  if ($("#metadata-path")?.value.trim() && !metadataContractValid) {
    throw new Error("Complete the metadata column mapping before starting the run.");
  }
  if (!outdir) throw new Error("Choose an output folder before starting the run.");
  if (!isAbsolutePath(outdir)) throw new Error("Output folder must be a full absolute path. Use Choose folder to select it.");
  return {
    fastq_dir: fastqDir,
    outdir,
    run_name: $("#run-name").value,
    metadata_path: $("#metadata-path") ? $("#metadata-path").value : "",
    metadata_sheet: $("#metadata-sheet") ? $("#metadata-sheet").value : "",
    metadata_date_order: $("#metadata-date-order")?.value || "auto",
    metadata_contract: metadataContractPayload(),
    kelt_barcode_map: $("#kelt-barcode-map")?.value || "",
    fallback_collection_year: $("#fallback-collection-year") ? $("#fallback-collection-year").value : "",
    fallback_collection_day: $("#fallback-collection-day") ? $("#fallback-collection-day").value : "27",
    libraries: selectedLibraries.length === availableLibraries.length ? [] : selectedLibraries,
    absolute_paths: true,
    write_samples: true,
    resume: false,
    dry_run: false,
    cpus: Number($("#cpus").value || 0),
    memory: $("#memory").value
  };
}

async function startRun() {
  saveSettings();
  $("#run-button").disabled = true;
  $("#run-message").className = "inline-message";
  text($("#run-message"), "Starting run...");
  try {
    const payload = await postJson("/api/run", runPayload());
    activeOutdir = payload.outdir;
    $("#results-outdir").value = payload.outdir;
    $("#run-samples").value = payload.samples || joinPath(payload.outdir, "samples.csv");
    text($("#sample-sheet-title"), pathBasename($("#run-samples").value, "samples.csv"));
    setPill($("#sample-sheet-status"), "Saved samples.csv", "ok");
    updateSamplePathHelp();
    $("#run-name").value = "";
    saveSettings();
    text($("#run-message"), payload.dry_run ? `Preview ready in ${payload.outdir}.` : `Run started in ${payload.outdir}.`);
    $("#run-message").classList.add("ok");
    if (payload.dry_run) {
      renderStages([], {status: "dry_run"}, {status: "dry_run"});
      await loadLog({forceScroll: true, statusPayload: {active: false, state: {status: "dry_run"}}});
    } else {
      showLaunchState(payload.outdir);
    }
    setStep("collect", true);
    startPolling();
    setTimeout(refreshAllRunState, 1200);
  } catch (error) {
    text($("#run-message"), error.message);
    $("#run-message").classList.add("bad");
  } finally {
    updateRunButtonAvailability();
  }
}

async function stopRun() {
  if (!activeOutdir) return;
  const confirmed = window.confirm("Stop the active malaria-amplicon-nf run? The output folder will stay on disk so you can resume later.");
  if (!confirmed) return;
  const button = $("#stop-button");
  const oldLabel = button ? button.textContent : "";
  if (button) {
    button.disabled = true;
    text(button, "Stopping...");
  }
  $("#run-message").className = "inline-message";
  text($("#run-message"), "Stopping run...");
  try {
    const payload = await postJson("/api/stop-run", {outdir: currentRunOutdir()});
    latestStatusPayload = payload;
    lastRunStatus = "stopped";
    renderStatus(payload);
    renderStages(payload.events || [], payload.summary || {status: "stopped"}, payload.state || {status: "stopped"});
    await loadLog({forceScroll: true, statusPayload: payload});
    text($("#run-message"), `Run stopped. To continue this output folder later, keep the same run name and enable resume.`);
    $("#run-message").classList.add("ok");
    await refreshAllRunState();
  } catch (error) {
    text($("#run-message"), error.message);
    $("#run-message").classList.add("bad");
  } finally {
    if (button) {
      text(button, oldLabel || "Stop run");
    }
  }
}

function showLaunchState(outdir) {
  followLog = true;
  lastRunStatus = "starting";
  displayedProgressPercent = 0;
  setProgressDisplay(0);
  renderStages([], {status: "starting", current_stage: "pending"}, {status: "starting"});
  renderTerminalLog(
    $("#technical-log"),
    `> Run started.\n> Output folder: ${outdir}\n> Waiting for first Nextflow update.`
  );
  scrollLogToBottom($("#technical-log"));
}

async function chooseFastqFolder() {
  const button = $("#browse-button");
  const scanButton = $("#scan-button");
  const oldLabel = button ? button.textContent : "";
  if (button) {
    button.disabled = true;
    button.classList.add("is-busy");
    text(button, "Choosing folder...");
  }
  if (scanButton && !scanInFlight) scanButton.disabled = true;
  setFolderMessage("Choosing a FASTQ folder...");
  const pickerHint = setTimeout(() => {
    setFolderMessage("Native folder picker is open. If you do not see it, check behind this browser or on the taskbar/Dock.", "ok");
  }, 1500);
  try {
    const payload = await withTimeout(
      chooseNativeFolder({
        initial: $("#fastq-dir").value,
        prompt: "Select the folder containing paired FASTQ files"
      }),
      NATIVE_PICKER_TIMEOUT_MS,
      "Folder picker timed out"
    );
    if (payload.selected && payload.path) {
      $("#fastq-dir").value = payload.path;
      invalidateScanReady();
      saveSettings();
      setFolderMessage("Folder selected. Click Scan folder when ready.", "ok");
      return;
    }
    setFolderMessage("Folder selection was cancelled.");
  } catch (error) {
    setFolderMessage("The system folder chooser could not open. Paste the full FASTQ folder path above, then scan.", "warn");
  } finally {
    clearTimeout(pickerHint);
    if (button) {
      button.disabled = false;
      button.classList.remove("is-busy");
      text(button, oldLabel || "Choose folder");
    }
    updateScanButtonAvailability();
  }
}

async function chooseMetadataFile() {
  const button = $("#choose-metadata-file");
  const oldLabel = button ? button.textContent : "";
  if (button) {
    button.disabled = true;
    text(button, "Opening...");
  }
  setMetadataMessage("Opening metadata file picker...");
  try {
    const payload = await withTimeout(
      chooseNativeFile({
        initial: $("#metadata-path").value || workspaceRoot || ".",
        prompt: "Select optional metadata file",
        kind: "metadata"
      }),
      NATIVE_PICKER_TIMEOUT_MS,
      "File picker timed out"
    );
    if (payload.selected && payload.path) {
      resetMetadataContract();
      $("#metadata-path").value = payload.path;
      invalidateScanReady();
      saveSettings();
      setMetadataMessage("Metadata file selected. Checking its structure...", "ok");
      await inspectMetadata();
      return;
    }
    setMetadataMessage("No metadata file selected. Leave this blank if metadata is not available.");
  } catch (error) {
    setMetadataMessage("Native file picker was not available. Type the metadata path manually.", "warn");
  } finally {
    if (button) {
      button.disabled = false;
      text(button, oldLabel || "Choose file");
    }
  }
}

function setKeltBarcodeMessage(message, tone = "") {
  const node = $("#kelt-barcode-help");
  if (!node) return;
  node.className = `field-help${tone ? ` ${tone}` : ""}`;
  text(node, message);
}

async function inspectKeltBarcodeMap() {
  const path = $("#kelt-barcode-map")?.value.trim() || "";
  const container = $("#kelt-barcode-inspection");
  if (!container) return;
  if (!path) {
    container.hidden = true;
    container.replaceChildren();
    setKeltBarcodeMessage("CSV or TSV with sample_id, forward_barcode, and reverse_barcode columns. Selecting a map enables exact barcode-pair matching.");
    return;
  }
  container.hidden = false;
  container.className = "metadata-inspection is-loading";
  text(container, "Checking expected KELT barcode pairs...");
  try {
    const payload = await postJson("/api/kelt/inspect", {path});
    container.className = "metadata-inspection";
    const libraryText = (payload.libraries || []).length
      ? ` across ${(payload.libraries || []).length} labelled libraries`
      : "";
    text(container, `${payload.pairs} sample barcode pairs loaded${libraryText}. Exact matching will run after adapter trimming.`);
    setKeltBarcodeMessage("KELT contamination QC will run for every selected sample.", "ok");
  } catch (error) {
    container.className = "metadata-inspection has-error";
    text(container, error.message || "The KELT barcode map could not be read.");
    setKeltBarcodeMessage("Fix or remove this map before starting the workflow.", "warn");
  }
}

async function chooseKeltBarcodeMap() {
  const button = $("#choose-kelt-barcode-map");
  const oldLabel = button?.textContent || "Choose file";
  if (button) {
    button.disabled = true;
    text(button, "Opening...");
  }
  setKeltBarcodeMessage("Opening KELT barcode map picker...");
  try {
    const payload = await withTimeout(
      chooseNativeFile({
        initial: $("#kelt-barcode-map")?.value || $("#fastq-dir")?.value || workspaceRoot || ".",
        prompt: "Select expected KELT barcode map",
        kind: "kelt"
      }),
      NATIVE_PICKER_TIMEOUT_MS,
      "File picker timed out"
    );
    if (payload.selected && payload.path) {
      $("#kelt-barcode-map").value = payload.path;
      saveSettings();
      await inspectKeltBarcodeMap();
      return;
    }
    setKeltBarcodeMessage("No KELT barcode map selected. This QC module remains off.");
  } catch (error) {
    setKeltBarcodeMessage("The system file chooser could not open. Paste the barcode map path above.", "warn");
  } finally {
    if (button) {
      button.disabled = false;
      text(button, oldLabel);
    }
  }
}

async function chooseOutputFolder() {
  const button = $("#choose-outdir-button");
  const oldLabel = button ? button.textContent : "";
  if (button) {
    button.disabled = true;
    text(button, "Opening...");
  }
  text($("#run-message"), "Opening output folder picker...");
  $("#run-message").className = "inline-message";
  try {
    const payload = await withTimeout(
      chooseNativeFolder({
        initial: $("#outdir").value,
        prompt: "Select the parent folder for run results",
        allowNewFolder: true
      }),
      NATIVE_PICKER_TIMEOUT_MS,
      "Folder picker timed out"
    );
    if (payload.selected && payload.path) {
      $("#outdir").value = payload.path;
      activeOutdir = "";
      syncGeneratedSampleSheetPath();
      updateSamplePathHelp();
      saveSettings();
      text($("#run-message"), "");
      $("#run-message").className = "inline-message";
      updateRunButtonAvailability();
      return;
    }
    text($("#run-message"), "Output folder selection was cancelled.");
  } catch (error) {
    text($("#run-message"), "The system folder chooser could not open. Paste the full output folder path above.");
    $("#run-message").classList.add("warn");
  } finally {
    if (button) {
      button.disabled = false;
      text(button, oldLabel || "Choose folder");
    }
  }
}

async function chooseResultsFolder() {
  const button = $("#choose-results-folder");
  const oldLabel = button ? button.textContent : "";
  if (button) {
    button.disabled = true;
    text(button, "Opening...");
  }
  try {
    const payload = await withTimeout(
      chooseNativeFolder({
        initial: $("#results-outdir").value || $("#outdir").value || workspaceRoot || ".",
        prompt: "Select a completed run results folder",
        allowNewFolder: false
      }),
      NATIVE_PICKER_TIMEOUT_MS,
      "Folder picker timed out"
    );
    if (payload.selected && payload.path) {
      $("#results-outdir").value = payload.path;
      saveSettings();
      await refreshResults();
    }
  } catch (error) {
    window.alert(error.message || "Folder picker was not available.");
  } finally {
    if (button) {
      button.disabled = false;
      text(button, oldLabel || "Choose folder");
    }
  }
}

async function openSelectedResultsFolder() {
  const button = $("#open-results-folder");
  const path = button?.dataset.path || $("#results-outdir")?.value || "";
  if (!path) return;
  const oldLabel = button.textContent;
  button.disabled = true;
  text(button, "Opening...");
  try {
    const payload = await postJson("/api/open-folder", {path}, {allowAppError: true});
    if (!payload.ok) throw new Error(payload.error || "Could not open folder.");
  } catch (error) {
    window.alert(error.message || "Could not open results folder.");
  } finally {
    button.disabled = false;
    text(button, oldLabel || "Open folder");
  }
}

function renderStages(events, summary, state) {
  const statusByStage = {};
  const messageByStage = {};
  const staleState = state?.status === "stale" || state?.status === "stopped";
  const runStatus = state?.status || summary.status || "pending";
  const hasStageEvents = events.some((event) => STAGES.includes(event.stage));
  STAGES.forEach((stage) => {
    statusByStage[stage] = "pending";
  });
  events.forEach((event) => {
    if (!STAGES.includes(event.stage)) return;
    if (staleState && ["started", "running"].includes(event.status)) return;
    statusByStage[event.stage] = event.status;
    if (event.message) messageByStage[event.stage] = event.message;
  });
  if (state && state.status === "failed") {
    const current = summary.current_stage;
    if (STAGES.includes(current) && statusByStage[current] !== "complete") {
      statusByStage[current] = "failed";
    } else if (!hasStageEvents && STAGES.length) {
      statusByStage[STAGES[0]] = "failed";
    }
  }
  if (isActiveStatus(runStatus) && !hasStageEvents && STAGES.length) {
    statusByStage[STAGES[0]] = "running";
    messageByStage[STAGES[0]] = "Launching workflow";
  }

  const list = $("#stage-list");
  STAGES.forEach((stage) => {
    let li = list.querySelector(`[data-stage="${stage}"]`);
    if (!li) {
      li = document.createElement("li");
      li.dataset.stage = stage;
      const dot = document.createElement("span");
      dot.className = "stage-dot";
      dot.setAttribute("aria-hidden", "true");
      const label = document.createElement("span");
      label.className = "stage-label";
      const stateNode = document.createElement("span");
      stateNode.className = "stage-status";
      li.append(dot, label, stateNode);
      list.appendChild(li);
    }
    const status = statusByStage[stage];
    li.className = stageClass(status);
    const label = li.querySelector(".stage-label");
    label.title = messageByStage[stage] || "";
    text(label, STAGE_LABELS[stage] || stage);
    const stateNode = li.querySelector(".stage-status");
    const statusLabel = stageStatusLabel(status);
    text(stateNode, statusLabel);
    stateNode.hidden = !statusLabel;
  });

  const completed = summary.completed_stages || 0;
  const total = summary.total_stages || STAGES.length;
  const percent = runStatus === "dry_run" || runStatus === "complete" ? 100 : Math.round((completed / Math.max(total, 1)) * 100);
  const indeterminate = isActiveStatus(runStatus) && percent === 0;
  const progressFill = $("#progress-fill");
  const progressBar = progressFill?.parentElement;
  const globalFill = $("#global-progress-fill");
  const globalBar = globalFill?.parentElement;
  [progressBar, globalBar].forEach((bar) => {
    if (!bar) return;
    bar.classList.toggle("is-running", isActiveStatus(runStatus));
    bar.classList.toggle("is-complete", runStatus === "complete" || runStatus === "dry_run");
    bar.classList.toggle("is-failed", runStatus === "failed");
    bar.classList.toggle("is-stopped", runStatus === "stopped");
    bar.classList.toggle("is-indeterminate", indeterminate);
  });
  if (progressBar) {
    if (indeterminate) progressBar.removeAttribute("aria-valuenow");
    else progressBar.setAttribute("aria-valuenow", String(percent));
  }
  animateProgressTo(percent);
  let pipelineLabel = checkStatusReady() ? "Ready" : "Idle";
  let pipelineClass = checkStatusReady() ? "is-ready" : "is-idle";
  let showPipelinePercent = false;
  if (isActiveStatus(runStatus)) {
    pipelineLabel = "Running";
    pipelineClass = "is-running";
    showPipelinePercent = true;
  }
  if (runStatus === "complete") pipelineLabel = "Complete";
  if (runStatus === "complete") {
    pipelineClass = "is-complete";
    showPipelinePercent = true;
  }
  if (runStatus === "dry_run") {
    pipelineLabel = "Preview ready";
    pipelineClass = "is-preview";
  }
  if (runStatus === "failed") {
    pipelineLabel = "Failed";
    pipelineClass = "is-failed";
    showPipelinePercent = true;
  }
  if (runStatus === "stopped") {
    pipelineLabel = "Stopped";
    pipelineClass = "is-unknown";
    showPipelinePercent = true;
  }
  if (!["pending", "starting", "running", "complete", "dry_run", "failed", "stopped"].includes(runStatus)) {
    pipelineLabel = "Check status";
    pipelineClass = "is-unknown";
  }
  text($("#pipeline-label"), pipelineLabel);

  const pipelineStatus = $(".pipeline-status");
  if (pipelineStatus) {
    pipelineStatus.hidden = !(
      isActiveStatus(runStatus)
      || ["complete", "dry_run", "failed", "stopped"].includes(runStatus)
    );
    pipelineStatus.classList.remove("is-idle", "is-ready", "is-running", "is-complete", "is-preview", "is-failed", "is-unknown", "has-percent");
    pipelineStatus.classList.add(pipelineClass);
    if (showPipelinePercent) pipelineStatus.classList.add("has-percent");
  }
  const current = STAGES.includes(summary.current_stage) ? summary.current_stage : "";
  let currentLabel = current ? (STAGE_LABELS[current] || current) : "No active run";
  if (isActiveStatus(runStatus) && !current) currentLabel = "Launching workflow";
  if (runStatus === "dry_run") currentLabel = "Preview ready";
  if (runStatus === "complete") currentLabel = "Run complete";
  if (runStatus === "failed") currentLabel = "Run failed";
  if (runStatus === "stopped") currentLabel = "Run stopped";
  if (runStatus === "stale") currentLabel = "Check run status";
  text($("#current-stage"), currentLabel);

  let progressCaption = "Waiting to start";
  if (isActiveStatus(runStatus)) {
    progressCaption = completed
      ? `${completed} of ${total} stages complete`
      : "Starting the workflow";
  }
  if (runStatus === "complete") progressCaption = "All workflow stages complete";
  if (runStatus === "dry_run") progressCaption = "Run preview prepared";
  if (runStatus === "failed") progressCaption = current ? `Stopped during ${STAGE_LABELS[current] || current}` : "Workflow stopped";
  if (runStatus === "stopped") progressCaption = "Stopped by user";
  text($("#progress-caption"), progressCaption);
  if ($("#progress-percent")) $("#progress-percent").hidden = true;

  const failureSummary = $("#run-failure-summary");
  if (failureSummary) {
    failureSummary.hidden = runStatus !== "failed";
    if (runStatus === "failed") {
      text($("#run-failure-detail"), state?.detail || "Reading the technical error details...");
    }
  }
}

function renderStatus(payload) {
  $("#stop-button").hidden = !Boolean(payload.active);
  $("#stop-button").disabled = !Boolean(payload.active);
  updateRunButtonAvailability(Boolean(payload.active));
}

function resetRunDisplay() {
  latestStatusPayload = null;
  lastRunStatus = "";
  completedRedirectKey = "";
  displayedProgressPercent = 0;
  setProgressDisplay(0);
  renderStatus({active: false, state: {status: "pending"}, summary: {status: "pending"}});
  renderStages([], {status: "pending", completed_stages: 0, total_stages: STAGES.length}, {status: "pending"});
  renderTerminalLog($("#technical-log"), EMPTY_LOG_TEXT);
}

async function fetchActiveRun() {
  return fetchJson("/api/active-run");
}

async function refreshStatus() {
  const out = encodeURIComponent(currentRunOutdir());
  const payload = await fetchJson(`/api/status?out=${out}`);
  renderStatus(payload);
  return payload;
}

async function refreshProgress(statusPayload = null) {
  const out = encodeURIComponent(currentRunOutdir());
  const payload = await fetchJson(`/api/progress?out=${out}`);
  const state = statusPayload ? statusPayload.state : {};
  renderStages(payload.events || [], payload.summary || {}, state || {});
  return payload;
}

async function refreshResults() {
  const out = activeOutdir || $("#results-outdir")?.value || $("#outdir").value || "results";
  if ($("#results-outdir")) $("#results-outdir").value = out;
  const payload = await fetchJson(`/api/results?out=${encodeURIComponent(out)}`);
  latestResultsPayload = payload;
  updateResultsRunSummary(payload, out);
  updateBundleButton(payload);
  renderQCDashboard(payload);
  await loadAnalysisInputStatus(encodeURIComponent(out));
  updateDinemitesRunButton();
  updateDciferRunButton();
  const list = $("#results-list");
  list.replaceChildren();
  if (!payload.files || !payload.files.length) {
    const div = document.createElement("div");
    div.className = "results-empty";
    text(div, "No completed outputs yet. Finish a run to see its report and tables here.");
    list.appendChild(div);
    return;
  }
  renderResultsDashboard(payload, list);
}

function updateResultsRunSummary(payload = {}, requestedOutdir = "") {
  const outdir = payload.outdir || requestedOutdir || $("#results-outdir")?.value || "";
  const hasCompletedOutput = Boolean(
    payload.report?.status === "ready"
    || Number(payload.ready_counts?.core || 0) > 0
  );
  const validRun = Boolean(payload.outdir_exists && hasCompletedOutput);
  const runName = validRun ? (payload.run_name || pathBasename(outdir, "Completed run")) : "No completed run";
  text($("#results-run-name"), runName);
  const openButton = $("#open-results-folder");
  if (openButton) {
    openButton.hidden = !validRun;
    openButton.disabled = !validRun;
    openButton.dataset.path = validRun ? outdir : "";
    openButton.title = validRun ? outdir : "A completed run will appear here automatically.";
  }
}

function updateBundleButton(payload) {
  const button = $("#download-bundle");
  if (!button) return;
  button.disabled = !payload.bundle_ready || !payload.bundle_url;
  button.dataset.url = payload.bundle_url || "";
  button.title = button.disabled ? "Bundle is available after result files are ready." : "Save a copy of the report and core result tables as a zip.";
}

function localSaveUrl(downloadUrl) {
  const url = new URL(downloadUrl, window.location.href);
  if (url.pathname === "/download-bundle") {
    return `/api/save-bundle?${url.searchParams.toString()}`;
  }
  const match = url.pathname.match(/^\/download\/([^/]+)$/);
  if (match) {
    return `/api/save-result/${encodeURIComponent(match[1])}?${url.searchParams.toString()}`;
  }
  return "";
}

function closeReportViewer() {
  const modal = $("#report-viewer-modal");
  if (!modal) return;
  const frame = modal.querySelector("iframe");
  if (frame) {
    frame.removeAttribute("src");
    frame.srcdoc = "";
  }
  modal.classList.remove("is-open");
  modal.setAttribute("aria-hidden", "true");
  if (!document.querySelector(".report-viewer-modal.is-open")) document.body.classList.remove("has-modal-open");
}

async function openReportViewer(url, label = "Run report") {
  let modal = $("#report-viewer-modal");
  if (!modal) {
    modal = document.createElement("div");
    modal.id = "report-viewer-modal";
    modal.className = "report-viewer-modal";
    modal.setAttribute("aria-hidden", "true");
    modal.innerHTML = `
      <div class="report-viewer-modal__scrim"></div>
      <section class="report-viewer-modal__panel" role="dialog" aria-modal="true" aria-label="Run report">
        <div class="report-viewer-modal__bar">
          <strong>Run report</strong>
          <button type="button" class="viewer-close-button report-viewer-modal__close" aria-label="Close report">&times;</button>
        </div>
        <iframe title="Run report"></iframe>
      </section>
    `;
    document.body.appendChild(modal);
    modal.querySelector(".report-viewer-modal__scrim").addEventListener("click", closeReportViewer);
    modal.querySelector(".report-viewer-modal__close").addEventListener("click", closeReportViewer);
  }
  const frame = modal.querySelector("iframe");
  text(modal.querySelector(".report-viewer-modal__bar strong"), label);
  modal.classList.add("is-open");
  modal.setAttribute("aria-hidden", "false");
  document.body.classList.add("has-modal-open");
  if (!frame) return;
  frame.removeAttribute("src");
  frame.srcdoc = "<!doctype html><html><body style=\"font-family: system-ui, sans-serif; padding: 24px; color: #005c68;\">Loading report...</body></html>";
  try {
    const response = await fetch(url, {credentials: "same-origin"});
    if (!response.ok) throw new Error(`Report returned HTTP ${response.status}`);
    frame.removeAttribute("srcdoc");
    frame.src = url;
  } catch (error) {
    frame.srcdoc = `<!doctype html><html><body style="font-family: system-ui, sans-serif; padding: 24px; color: #c94b42;"><strong>Could not open report.</strong><p>${escapeHtml(error.message || "Report failed to load.")}</p></body></html>`;
  }
}

function closeImageViewer() {
  const modal = $("#image-viewer-modal");
  if (!modal) return;
  const image = modal.querySelector("img");
  if (image) image.removeAttribute("src");
  modal.classList.remove("is-open");
  modal.setAttribute("aria-hidden", "true");
  if (!document.querySelector(".report-viewer-modal.is-open")) document.body.classList.remove("has-modal-open");
}

function openImageViewer(url, label = "Full-resolution plot") {
  let modal = $("#image-viewer-modal");
  if (!modal) {
    modal = document.createElement("div");
    modal.id = "image-viewer-modal";
    modal.className = "report-viewer-modal image-viewer-modal";
    modal.setAttribute("aria-hidden", "true");
    modal.innerHTML = `
      <div class="report-viewer-modal__scrim"></div>
      <section class="report-viewer-modal__panel image-viewer-modal__panel" role="dialog" aria-modal="true" aria-label="Full-resolution plot">
        <div class="report-viewer-modal__bar">
          <strong></strong>
          <div class="image-viewer-modal__actions">
            <button type="button" class="image-viewer-fit">Fit to window</button>
            <button type="button" class="viewer-close-button image-viewer-close" aria-label="Close plot">&times;</button>
          </div>
        </div>
        <div class="image-viewer-modal__body"><img alt=""></div>
      </section>
    `;
    document.body.appendChild(modal);
    modal.querySelector(".report-viewer-modal__scrim").addEventListener("click", closeImageViewer);
    modal.querySelector(".image-viewer-close").addEventListener("click", closeImageViewer);
    modal.querySelector(".image-viewer-fit").addEventListener("click", (event) => {
      const fit = modal.classList.toggle("is-fit");
      text(event.currentTarget, fit ? "Actual size" : "Fit to window");
    });
  }
  modal.classList.remove("is-fit");
  text(modal.querySelector(".report-viewer-modal__bar strong"), label);
  text(modal.querySelector(".image-viewer-fit"), "Fit to window");
  const image = modal.querySelector("img");
  image.src = url;
  image.alt = label;
  modal.classList.add("is-open");
  modal.setAttribute("aria-hidden", "false");
  document.body.classList.add("has-modal-open");
}

function closeTableViewer() {
  const modal = $("#table-viewer-modal");
  if (!modal) return;
  modal.classList.remove("is-open");
  modal.setAttribute("aria-hidden", "true");
  if (!document.querySelector(".report-viewer-modal.is-open")) document.body.classList.remove("has-modal-open");
  const body = modal.querySelector(".table-viewer-modal__body");
  if (body) body.replaceChildren();
}

function renderTablePreview(payload, container) {
  container.replaceChildren();
  const meta = document.createElement("p");
  meta.className = "table-viewer-meta";
  const columnCount = (payload.columns || []).length;
  text(meta, `${payload.shown_rows || 0} rows${payload.truncated ? " previewed" : ""} - ${columnCount} columns - ${payload.size || ""}`);
  const wrap = document.createElement("div");
  wrap.className = "table-viewer-table-wrap";
  const table = document.createElement("table");
  table.className = "table-viewer-table";
  const thead = document.createElement("thead");
  const headRow = document.createElement("tr");
  (payload.columns || []).forEach((column) => {
    const th = document.createElement("th");
    text(th, column);
    th.title = column;
    headRow.appendChild(th);
  });
  thead.appendChild(headRow);
  const tbody = document.createElement("tbody");
  (payload.rows || []).forEach((dataRow) => {
    const tr = document.createElement("tr");
    (payload.columns || []).forEach((column) => {
      const td = document.createElement("td");
      const value = dataRow[column] ?? "";
      text(td, value);
      td.title = value;
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });
  if (!(payload.rows || []).length) {
    const tr = document.createElement("tr");
    const td = document.createElement("td");
    td.colSpan = Math.max((payload.columns || []).length, 1);
    td.className = "empty";
    text(td, "No rows to show.");
    tr.appendChild(td);
    tbody.appendChild(tr);
  }
  table.append(thead, tbody);
  wrap.appendChild(table);
  container.append(meta, wrap);
}

async function openTableViewer(url, label) {
  let modal = $("#table-viewer-modal");
  if (!modal) {
    modal = document.createElement("div");
    modal.id = "table-viewer-modal";
    modal.className = "report-viewer-modal table-viewer-modal";
    modal.setAttribute("aria-hidden", "true");
    modal.innerHTML = `
      <div class="report-viewer-modal__scrim"></div>
      <section class="report-viewer-modal__panel table-viewer-modal__panel" role="dialog" aria-modal="true" aria-label="Result table">
        <div class="report-viewer-modal__bar">
          <strong></strong>
          <button type="button" class="viewer-close-button table-viewer-modal__close" aria-label="Close table">&times;</button>
        </div>
        <div class="table-viewer-modal__body"></div>
      </section>
    `;
    document.body.appendChild(modal);
    modal.querySelector(".report-viewer-modal__scrim").addEventListener("click", closeTableViewer);
    modal.querySelector(".table-viewer-modal__close").addEventListener("click", closeTableViewer);
  }
  text(modal.querySelector(".report-viewer-modal__bar strong"), label || "Result table");
  const body = modal.querySelector(".table-viewer-modal__body");
  modal.classList.add("is-open");
  modal.setAttribute("aria-hidden", "false");
  document.body.classList.add("has-modal-open");
  if (!body) return;
  body.innerHTML = '<p class="table-viewer-meta">Loading table...</p>';
  try {
    const payload = await fetchJson(url);
    renderTablePreview(payload, body);
  } catch (error) {
    body.innerHTML = `<p class="table-viewer-error">Could not load table: ${escapeHtml(error.message || "Unknown error")}</p>`;
  }
}

function renderASVAudit(payload, container) {
  container.replaceChildren();
  const steps = payload.steps || [];
  if (!steps.length) {
    const empty = document.createElement("p");
    empty.className = "table-viewer-error";
    text(empty, "No filtering stages were found in this audit.");
    container.appendChild(empty);
    return;
  }

  const intro = document.createElement("div");
  intro.className = "asv-audit-intro";
  const title = document.createElement("h2");
  const copy = document.createElement("p");
  text(title, "ASV retention through pipeline QC");
  text(copy, "Counts are cumulative. Each stage starts with the ASVs retained by the preceding stage.");
  intro.append(title, copy);

  const funnel = document.createElement("div");
  funnel.className = "asv-audit-funnel";
  const baseline = Number(steps[0].count || 0);
  const readDepth = steps.find((step) => /read-depth/i.test(step.step));
  const recurrence = steps.find((step) => /recurrence/i.test(step.step));
  const exactCigar = Number(payload.cigar?.unique_cigar_haplotypes_after_conversion || steps.at(-1)?.count || 0);
  const milestones = [
    {label: "Discovered ASVs", count: baseline, detail: "Starting sequence features"},
    {label: "Passed read depth", count: Number(readDepth?.count || 0), detail: readDepth?.requirement || "Read-depth requirement"},
    {label: "Passed recurrence", count: Number(recurrence?.count || 0), detail: recurrence?.requirement || "Recurrence requirement"},
    {label: "Final exact-CIGAR alleles", count: exactCigar, detail: "Primary allele result"},
  ];
  milestones.forEach((milestone, index) => {
    const stage = document.createElement("div");
    stage.className = "asv-audit-stage";
    const label = document.createElement("span");
    const count = document.createElement("strong");
    const retained = document.createElement("small");
    const track = document.createElement("div");
    const fill = document.createElement("i");
    track.className = "asv-audit-stage__track";
    fill.style.width = `${Math.max(2, baseline ? 100 * milestone.count / baseline : 0)}%`;
    text(label, milestone.label);
    text(count, milestone.count.toLocaleString());
    text(retained, index === 0 ? milestone.detail : `${(baseline ? 100 * milestone.count / baseline : 0).toFixed(1)}% of start`);
    stage.title = milestone.detail;
    track.appendChild(fill);
    stage.append(label, count, retained, track);
    funnel.appendChild(stage);
  });

  const heading = document.createElement("h3");
  heading.className = "asv-audit-table-title";
  text(heading, "Filtering requirements");
  const tableWrap = document.createElement("div");
  tableWrap.className = "asv-audit-table-wrap";
  const table = document.createElement("table");
  table.className = "asv-audit-table";
  table.innerHTML = "<thead><tr><th>Stage</th><th>Requirement</th><th>Retained</th><th>Removed here</th><th>From previous</th></tr></thead>";
  const tbody = document.createElement("tbody");
  steps.forEach((step, index) => {
    const row = document.createElement("tr");
    [
      step.step,
      step.requirement || "-",
      Number(step.count || 0).toLocaleString(),
      index === 0 ? "-" : Number(step.removed || 0).toLocaleString(),
      index === 0 ? "Baseline" : `${Number(step.retained_previous_pct || 0).toFixed(1)}%`,
    ].forEach((value) => {
      const cell = document.createElement("td");
      text(cell, value);
      row.appendChild(cell);
    });
    tbody.appendChild(row);
  });
  table.appendChild(tbody);
  tableWrap.appendChild(table);

  const sensitivity = document.createElement("section");
  sensitivity.className = "asv-audit-sensitivity";
  const sensitivityHeading = document.createElement("div");
  const sensitivityTitle = document.createElement("h3");
  const sensitivityCopy = document.createElement("p");
  text(sensitivityTitle, "Optional CD-HIT comparison");
  text(sensitivityCopy, "A separate sensitivity result. The primary exact-CIGAR table is not reduced or overwritten.");
  sensitivityHeading.append(sensitivityTitle, sensitivityCopy);
  const sensitivityFlow = document.createElement("div");
  sensitivityFlow.className = "asv-audit-sensitivity__flow";
  const cdhit = payload.cdhit || {};
  const clusters = Number(cdhit.output_clusters || 0);
  const consolidated = Number(cdhit.sequence_features_grouped_away_by_clustering || Math.max(0, exactCigar - clusters));
  const exactBlock = document.createElement("div");
  const arrow = document.createElement("span");
  const clusterBlock = document.createElement("div");
  const consolidation = document.createElement("p");
  exactBlock.innerHTML = `<strong>${exactCigar.toLocaleString()}</strong><span>Exact ASVs</span>`;
  arrow.className = "asv-audit-sensitivity__arrow";
  arrow.setAttribute("aria-hidden", "true");
  text(arrow, "→");
  clusterBlock.innerHTML = `<strong>${clusters.toLocaleString()}</strong><span>Similarity clusters</span>`;
  text(consolidation, `${consolidated.toLocaleString()} ASVs consolidated into existing clusters.`);
  sensitivityFlow.append(exactBlock, arrow, clusterBlock);
  sensitivity.append(sensitivityHeading, sensitivityFlow, consolidation);
  container.append(intro, funnel, heading, tableWrap, sensitivity);
}

async function openASVAuditViewer(url) {
  let modal = $("#asv-audit-modal");
  if (!modal) {
    modal = document.createElement("div");
    modal.id = "asv-audit-modal";
    modal.className = "report-viewer-modal asv-audit-modal";
    modal.setAttribute("aria-hidden", "true");
    modal.innerHTML = `
      <div class="report-viewer-modal__scrim"></div>
      <section class="report-viewer-modal__panel asv-audit-modal__panel" role="dialog" aria-modal="true" aria-label="ASV filtering audit">
        <div class="report-viewer-modal__bar">
          <strong>ASV filtering audit</strong>
          <button type="button" class="viewer-close-button asv-audit-modal__close" aria-label="Close audit">&times;</button>
        </div>
        <div class="asv-audit-modal__body"></div>
      </section>`;
    document.body.appendChild(modal);
    const close = () => {
      modal.classList.remove("is-open");
      modal.setAttribute("aria-hidden", "true");
      document.body.classList.remove("has-modal-open");
    };
    modal.querySelector(".report-viewer-modal__scrim").addEventListener("click", close);
    modal.querySelector(".asv-audit-modal__close").addEventListener("click", close);
  }
  modal.classList.add("is-open");
  modal.setAttribute("aria-hidden", "false");
  document.body.classList.add("has-modal-open");
  const body = modal.querySelector(".asv-audit-modal__body");
  body.innerHTML = '<p class="table-viewer-meta">Loading filtering audit...</p>';
  try {
    renderASVAudit(await fetchJson(url), body);
  } catch (error) {
    body.innerHTML = `<p class="table-viewer-error">Could not load audit: ${escapeHtml(error.message || "Unknown error")}</p>`;
  }
}

async function downloadFile(url, button) {
  const oldText = button ? button.textContent : "";
  const oldTitle = button ? button.title : "";
  if (button) {
    button.disabled = true;
    text(button, "Saving...");
  }
  try {
    const saveUrl = localSaveUrl(url);
    if (!saveUrl) {
      window.location.href = url;
      return;
    }
    const response = await fetch(saveUrl, {method: "POST"});
    const payload = await response.json().catch(() => ({}));
    if (!response.ok || payload.ok === false) {
      throw new Error(payload.error || `Save failed: ${response.status}`);
    }
    if (button) {
      text(button, "Saved");
      button.title = payload.path ? `Saved to ${payload.path}` : "Saved to Downloads";
      window.setTimeout(() => {
        text(button, oldText || "Save");
        button.title = oldTitle;
      }, 1800);
    }
  } catch (error) {
    window.alert(error.message || "Save failed.");
  } finally {
    if (button) {
      button.disabled = false;
      if (button.textContent !== "Saved") text(button, oldText || "Save");
    }
  }
}

function resultAction(url, label, primary = false, action = "save", meta = {}) {
  const button = document.createElement("button");
  button.type = "button";
  if (primary) button.className = "primary";
  text(button, label);
  button.addEventListener("click", () => {
    if (action === "report") {
      openReportViewer(url, meta.label || label);
      return;
    }
    if (action === "table") {
      openTableViewer(url, meta.label || label);
      return;
    }
    downloadFile(url, button);
  });
  return button;
}

function resultTableRow(item) {
  const row = document.createElement("tr");
  const fileCell = document.createElement("td");
  fileCell.className = "result-file-cell";
  const title = document.createElement("strong");
  text(title, item.label);
  fileCell.appendChild(title);

  const sizeCell = document.createElement("td");
  text(sizeCell, item.size);

  const actionsCell = document.createElement("td");
  const actions = document.createElement("div");
  actions.className = "result-actions";
  if (item.table_url) actions.appendChild(resultAction(item.table_url, "View table", true, "table", {label: item.label}));
  if (item.view_url) actions.appendChild(resultAction(item.view_url, "View", true, "report"));
  if (item.download_url) actions.appendChild(resultAction(item.download_url, "Download", false, "download"));
  actionsCell.appendChild(actions);

  row.append(fileCell, sizeCell, actionsCell);
  return row;
}

function resultTable(files, emptyText) {
  const section = document.createElement("section");
  section.className = "results-table-section";
  const wrap = document.createElement("div");
  wrap.className = "results-table-wrap";
  const table = document.createElement("table");
  table.className = "results-table";
  const thead = document.createElement("thead");
  thead.innerHTML = "<tr><th>Table</th><th>Size</th><th></th></tr>";
  const tbody = document.createElement("tbody");
  if (!files.length) {
    const row = document.createElement("tr");
    const cell = document.createElement("td");
    cell.colSpan = 3;
    cell.className = "empty";
    text(cell, emptyText);
    row.appendChild(cell);
    tbody.appendChild(row);
  } else {
    files.forEach((item) => tbody.appendChild(resultTableRow(item)));
  }
  table.append(thead, tbody);
  wrap.appendChild(table);
  section.appendChild(wrap);
  return section;
}

function reportMetric(metric) {
  const card = document.createElement("article");
  card.className = "results-summary-metric";
  const label = document.createElement("span");
  const value = document.createElement("strong");
  const detail = document.createElement("small");
  const spans = metric.querySelectorAll("span");
  text(label, spans[0]?.textContent?.trim() || "Metric");
  text(value, metric.querySelector("b")?.textContent?.trim() || "-");
  text(detail, spans[spans.length - 1]?.textContent?.trim() || "");
  card.append(label, value, detail);
  return card;
}

function reportChart(chart) {
  const panel = document.createElement("article");
  panel.className = "results-summary-chart";
  const heading = document.createElement("h4");
  text(heading, chart.querySelector("h3")?.textContent?.trim() || "Quality summary");
  panel.appendChild(heading);

  chart.querySelectorAll(".chart-row").forEach((sourceRow) => {
    const row = document.createElement("div");
    row.className = "results-summary-chart-row";
    const label = document.createElement("span");
    label.className = "results-summary-chart-label";
    text(label, sourceRow.querySelector(".chart-label")?.textContent?.trim() || "-");

    const track = document.createElement("span");
    track.className = "results-summary-chart-track";
    const fill = document.createElement("i");
    const sourceFill = sourceRow.querySelector(".chart-track span");
    const width = Number.parseFloat(sourceFill?.style?.width || "0");
    fill.style.width = `${Math.min(100, Math.max(0, Number.isFinite(width) ? width : 0))}%`;
    track.appendChild(fill);

    const value = document.createElement("strong");
    text(value, sourceRow.querySelector(".chart-value")?.textContent?.trim() || "-");
    row.append(label, track, value);
    panel.appendChild(row);
  });
  return panel;
}

async function renderInlineRunSummary(url, container) {
  container.replaceChildren();
  const loading = document.createElement("p");
  loading.className = "results-summary-state";
  text(loading, "Loading run summary...");
  container.appendChild(loading);

  try {
    const response = await fetch(url, {credentials: "same-origin"});
    if (!response.ok) throw new Error(`Run report returned HTTP ${response.status}`);
    const report = new DOMParser().parseFromString(await response.text(), "text/html");
    const metricNodes = Array.from(report.querySelectorAll(".metrics .metric"));
    const chartNodes = Array.from(report.querySelectorAll(".viz-grid .chart"));
    if (!metricNodes.length && !chartNodes.length) throw new Error("Run summary data was not found in the report.");

    container.replaceChildren();
    if (metricNodes.length) {
      const metrics = document.createElement("div");
      metrics.className = "results-summary-metrics";
      metricNodes.forEach((metric) => metrics.appendChild(reportMetric(metric)));
      container.appendChild(metrics);
    }
    if (chartNodes.length) {
      const heading = document.createElement("div");
      heading.className = "results-summary-subheading";
      const title = document.createElement("h3");
      const copy = document.createElement("span");
      text(title, "Quality overview");
      text(copy, "Read depth, retained variants, and locus recovery");
      heading.append(title, copy);
      const charts = document.createElement("div");
      charts.className = "results-summary-charts";
      chartNodes.forEach((chart) => charts.appendChild(reportChart(chart)));
      container.append(heading, charts);
    }
  } catch (error) {
    container.replaceChildren();
    const unavailable = document.createElement("p");
    unavailable.className = "results-summary-state is-error";
    text(unavailable, `The run summary could not be displayed here. ${error.message || ""}`.trim());
    container.appendChild(unavailable);
  }
}

function renderResultsDashboard(payload, list) {
  const report = payload.report;
  const coreFiles = (payload.core_files || []).filter((item) => item.status === "ready" && item.label !== ASV_FILTERING_SUMMARY_LABEL);
  const keltReport = payload.kelt_report;
  const keltFiles = (payload.kelt_files || []).filter((item) => item.status === "ready");
  const reportReady = Boolean(report && report.status === "ready" && report.view_url);
  const keltReportReady = Boolean(keltReport && keltReport.status === "ready" && keltReport.view_url);

  const dashboard = document.createElement("div");
  dashboard.className = "results-dashboard";

  if (!reportReady && !coreFiles.length && !keltReportReady && !keltFiles.length) {
    const empty = document.createElement("div");
    empty.className = "results-empty";
    text(empty, "No completed outputs yet. Report and table links appear here as the workflow finishes.");
    dashboard.appendChild(empty);
    list.appendChild(dashboard);
    return;
  }

  if (reportReady) {
    const summary = document.createElement("section");
    summary.className = "results-inline-summary";
    const heading = document.createElement("div");
    heading.className = "results-inline-summary__heading";
    const title = document.createElement("h3");
    const detail = document.createElement("span");
    text(title, "Run summary");
    text(detail, "Key statistics and quality-control patterns from this run");
    heading.append(title, detail);
    const body = document.createElement("div");
    body.className = "results-inline-summary__body";
    summary.append(heading, body);
    dashboard.appendChild(summary);
    renderInlineRunSummary(report.view_url, body);
  }

  if (coreFiles.length) {
    const tableHeader = document.createElement("div");
    tableHeader.className = "results-section-heading";
    const heading = document.createElement("h3");
    text(heading, "Output tables");
    const count = document.createElement("span");
    text(count, `${coreFiles.length} ${coreFiles.length === 1 ? "table" : "tables"}`);
    tableHeader.append(heading, count);
    dashboard.append(tableHeader, resultTable(coreFiles, "No completed output tables yet."));
  }

  if (keltReportReady || keltFiles.length) {
    const keltHeader = document.createElement("div");
    keltHeader.className = "results-section-heading results-section-heading--kelt";
    const heading = document.createElement("h3");
    text(heading, "KELT contamination QC");
    const detail = document.createElement("span");
    text(detail, "Well-specific inline-barcode checks");
    keltHeader.append(heading, detail);
    dashboard.appendChild(keltHeader);
    if (keltReportReady) {
      const actions = document.createElement("div");
      actions.className = "results-kelt-report";
      const copy = document.createElement("span");
      text(copy, "Contamination summary and samples requiring review");
      actions.append(copy, resultAction(keltReport.view_url, "View report", true, "report", {label: "KELT contamination QC"}));
      dashboard.appendChild(actions);
    }
    if (keltFiles.length) dashboard.appendChild(resultTable(keltFiles, "No KELT QC tables are available."));
  }

  list.appendChild(dashboard);
}

function updateQCWorkspaceState(hasRun = hasSelectedRunOutputs()) {
  const emptyState = $("#qc-empty-state");
  const workspace = $("#qc-workspace");
  if (emptyState) emptyState.hidden = hasRun;
  if (workspace) workspace.hidden = !hasRun;
  ["dinemites", "dcifer"].forEach((name) => {
    const button = document.querySelector(`.downstream-subtab[data-downstream-tab="${name}"]`);
    if (!button) return;
    button.disabled = !hasRun;
    button.title = hasRun ? "" : "Select a completed run before starting an analysis.";
  });
}

function renderQCDashboard(payload = {}) {
  const list = $("#qc-results-list");
  if (!list) return;
  list.replaceChildren();
  updateQCWorkspaceState();

  const asvSummary = (payload.core_files || []).find((item) => {
    return item.label === ASV_FILTERING_SUMMARY_LABEL && item.status === "ready";
  });
  const cdhitFiles = (payload.cdhit_files || []).filter((item) => item.status === "ready");
  const cdhitSummary = payload.cdhit_summary || {};

  const dashboard = document.createElement("div");
  dashboard.className = "results-dashboard";

  if (!asvSummary && !cdhitFiles.length) {
    const empty = document.createElement("div");
    empty.className = "results-empty";
    text(empty, "QC files are still being prepared.");
    dashboard.appendChild(empty);
    list.appendChild(dashboard);
    return;
  }

  if (asvSummary) {
    const header = document.createElement("div");
    header.className = "results-section-heading";
    const heading = document.createElement("h3");
    text(heading, "ASV filtering");
    const detail = document.createElement("span");
    text(detail, "Cumulative audit");
    header.append(heading, detail);
    const auditRow = document.createElement("div");
    auditRow.className = "qc-audit-row";
    const auditCopy = document.createElement("div");
    const auditTitle = document.createElement("strong");
    const auditDetail = document.createElement("span");
    text(auditTitle, "Filtering audit");
    text(auditDetail, "Starting ASVs, each requirement, and retained counts");
    auditCopy.append(auditTitle, auditDetail);
    const actions = document.createElement("div");
    actions.className = "result-actions";
    const auditButton = document.createElement("button");
    auditButton.type = "button";
    auditButton.className = "primary";
    text(auditButton, "View audit");
    auditButton.addEventListener("click", () => {
      openASVAuditViewer(`/api/asv-filtering-audit?out=${encodeURIComponent(dinemitesOutdir())}`);
    });
    actions.appendChild(auditButton);
    if (asvSummary.download_url) actions.appendChild(resultAction(asvSummary.download_url, "Save copy", false, "download"));
    auditRow.append(auditCopy, actions);
    dashboard.append(header, auditRow);
  }

  if (cdhitFiles.length) {
    const identity = Number(cdhitSummary.identity_threshold || 0.989);
    const header = document.createElement("div");
    header.className = "results-section-heading results-section-heading--cdhit";
    const heading = document.createElement("h3");
    text(heading, "CD-HIT clustering output");
    const detail = document.createElement("span");
    text(detail, "Created by the pipeline");
    header.append(heading, detail);
    dashboard.appendChild(header);

    const sourceNote = document.createElement("div");
    sourceNote.className = "qc-cdhit-source-note";
    const sourceTitle = document.createElement("strong");
    const sourceCopy = document.createElement("span");
    text(sourceTitle, `${(identity * 100).toLocaleString(undefined, {maximumFractionDigits: 1})}% sequence identity`);
    text(sourceCopy, "These are raw ASV-to-cluster assignments and counts for each sequencing library. Participant-level abundance and replicate rules have not been applied yet.");
    sourceNote.append(sourceTitle, sourceCopy);
    dashboard.appendChild(sourceNote);

    const metrics = document.createElement("div");
    metrics.className = "results-cdhit-summary";
    [
      [cdhitSummary.input_asvs, "ASVs clustered"],
      [cdhitSummary.clusters, "Clusters"],
      [cdhitSummary.singleton_clusters, "Singleton clusters"],
      [cdhitSummary.largest_cluster_size, "Largest cluster"],
    ].forEach(([value, label]) => {
      const metric = document.createElement("div");
      const strong = document.createElement("strong");
      const caption = document.createElement("span");
      text(strong, Number.isFinite(Number(value)) ? Number(value).toLocaleString() : "-");
      text(caption, label);
      metric.append(strong, caption);
      metrics.appendChild(metric);
    });
    dashboard.append(metrics, resultTable(cdhitFiles, "No CD-HIT cluster tables are available."));
  }

  list.appendChild(dashboard);
}

function analysisInputMatchesSettings(input = {}) {
  if (!input.available) return false;
  const filter = input.summary?.abundance_filter || {};
  return Math.abs(Number(filter.threshold_percent) - analysisMinAbundancePct()) < 0.000001
    && String(filter.denominator || "") === analysisAbundanceDenominator()
    && String(input.summary?.replicate_policy || "") === "strict";
}

async function buildAnalysisTable(mode = "cdhit98") {
  const msg = $("#analysis-input-message");
  if (!msg) return;
  msg.className = "inline-message";
  saveSettings();
  const primary = mode === "primary";
  const button = primary ? $("#analysis-build-primary") : $("#analysis-build-table");
  if (button) button.disabled = true;
  text(msg, primary
    ? "Confirming alleles across technical replicates and applying the biological call rules..."
    : "Preparing a sensitivity input from the existing 98.9% CD-HIT clusters...");
  try {
    const payload = await postJson("/api/analysis-table/build", {
      outdir: dinemitesOutdir(),
      mode,
      min_abundance_pct: analysisMinAbundancePct(),
      abundance_denominator: analysisAbundanceDenominator()
    });
    updateAnalysisInputAvailability(payload.inputs || {});
    msg.classList.add("ok");
    const summary = payload.summary || {};
    const output = summary.output || {};
    if (primary) {
      text(msg, `Exact-allele input ready: ${Number(output.participant_visits || 0).toLocaleString()} participant visits and ${Number(output.alleles || 0).toLocaleString()} alleles.`);
    } else {
      const clusters = Number(summary.clusters || output.alleles || 0).toLocaleString();
      text(msg, `CD-HIT input ready: ${clusters} locus-specific clusters after analysis filtering and replicate merging. CD-HIT was not rerun.`);
    }
  } catch (error) {
    msg.classList.add("bad");
    text(msg, userMessage(error.message));
  } finally {
    if (button) button.disabled = false;
  }
}

function updateAnalysisInputAvailability(inputs = {}) {
  latestAnalysisInputStatus = inputs;
  const primary = inputs.primary || {};
  const sensitivity = inputs.cdhit98 || {};
  const primaryReady = analysisInputMatchesSettings(primary);
  const available = analysisInputMatchesSettings(sensitivity);
  const hasRun = hasSelectedRunOutputs();
  updateQCWorkspaceState(hasRun);
  ["#dinemites-analysis-mode", "#dcifer-analysis-mode"].forEach((selector) => {
    const select = $(selector);
    const option = select?.querySelector('option[value="summed"]');
    if (option) option.disabled = !available;
  });

  const primaryStatus = $("#analysis-primary-status");
  if (primaryStatus) {
    primaryStatus.className = `analysis-input-state${primaryReady ? " is-ready" : ""}`;
    if (primaryReady) {
      const output = primary.summary?.output || {};
      text(primaryStatus, `Ready - ${Number(output.participant_visits || primary.samples || 0).toLocaleString()} visits, ${Number(output.alleles || primary.alleles || 0).toLocaleString()} alleles`);
    } else {
      text(primaryStatus, primary.available ? "Settings changed - rebuild required" : "Not built");
    }
  }
  const primaryBuild = $("#analysis-build-primary");
  if (primaryBuild) {
    primaryBuild.disabled = !hasRun;
    text(primaryBuild, primaryReady ? "Update exact-allele input" : "Prepare exact-allele input");
  }
  const primaryView = $("#analysis-view-primary");
  const primaryDownload = $("#analysis-download-primary");
  if (primaryView) {
    primaryView.hidden = !primaryReady;
    primaryView.disabled = !primaryReady;
    primaryView.onclick = primaryReady && primary.table_url
      ? () => openTableViewer(primary.table_url, "Shared exact-CIGAR analysis table")
      : null;
  }
  if (primaryDownload) {
    primaryDownload.hidden = !primaryReady;
    primaryDownload.disabled = !primaryReady;
    primaryDownload.onclick = primaryReady && primary.download_url
      ? () => downloadFile(primary.download_url, primaryDownload)
      : null;
  }

  const status = $("#analysis-cdhit-status");
  if (status) {
    status.className = `analysis-input-state${available ? " is-ready" : ""}`;
    if (available) {
      const alleles = Number(sensitivity.alleles || sensitivity.summary?.clusters || 0).toLocaleString();
      const samples = Number(sensitivity.samples || 0).toLocaleString();
      text(status, `Ready - ${alleles} clusters, ${samples} samples`);
    } else {
      text(status, "Not built");
    }
  }
  const buildButton = $("#analysis-build-table");
  if (buildButton) {
    buildButton.disabled = !hasRun;
    text(buildButton, available ? "Update CD-HIT input" : "Prepare CD-HIT input");
  }
  const viewButton = $("#analysis-view-table");
  const downloadButton = $("#analysis-download-table");
  if (viewButton) {
    viewButton.disabled = !available;
    viewButton.hidden = !available;
    viewButton.onclick = available && sensitivity.table_url
      ? () => openTableViewer(sensitivity.table_url, "CD-HIT 98.9% sensitivity input")
      : null;
  }
  if (downloadButton) {
    downloadButton.disabled = !available;
    downloadButton.hidden = !available;
    downloadButton.onclick = available && sensitivity.download_url
      ? () => downloadFile(sensitivity.download_url, downloadButton)
      : null;
  }
  ["#analysis-open-cdhit-dinemites", "#analysis-open-cdhit-dcifer"].forEach((selector) => {
    const button = $(selector);
    if (!button) return;
    button.disabled = !available;
    button.hidden = !available;
    button.title = available ? "" : "Prepare the CD-HIT 98.9% sensitivity input first.";
  });
  const cdhitLaunchActions = $("#analysis-cdhit-launch-actions");
  if (cdhitLaunchActions) cdhitLaunchActions.hidden = !available;

  ["#qc-continue-dinemites", "#qc-continue-dcifer"].forEach((selector) => {
    const button = $(selector);
    if (button) {
      button.disabled = !primaryReady;
      button.hidden = !primaryReady;
      button.title = primaryReady ? "" : "Prepare the exact-allele input first.";
    }
  });
  const primaryLaunchActions = $("#analysis-primary-launch-actions");
  if (primaryLaunchActions) primaryLaunchActions.hidden = !primaryReady;
  updateDinemitesRunButton();
  updateDciferRunButton();
}

async function loadAnalysisInputStatus(encodedOutdir = "") {
  const out = encodedOutdir || encodeURIComponent(dinemitesOutdir());
  try {
    const payload = await fetchJson(`/api/analysis-table/status?out=${out}`);
    updateAnalysisInputAvailability(payload.inputs || {});
  } catch (_error) {
    updateAnalysisInputAvailability({});
  }
}

async function refreshAllRunState() {
  if (!activeOutdir) {
    resetRunDisplay();
    return {active: false};
  }
  const previousStatus = lastRunStatus;
  try {
    const status = await refreshStatus();
    latestStatusPayload = status;
    try {
      await refreshProgress(status);
    } catch (_error) {
      // Keep polling even if one progress read races the file writer.
    }
    try {
      await refreshResults();
    } catch (_error) {
      // Results may not exist yet while the pipeline is still active.
    }
    await loadLog({silent: true, statusPayload: status});
    const currentStatus = payloadStatus(status);
    const wasRunning = isActiveStatus(previousStatus);
    if (isActiveStatus(currentStatus) && !pollTimer) {
      startPolling();
    }
    if ((currentStatus === "complete" || currentStatus === "dry_run") && wasRunning) {
      const state = status.state || {};
      const key = `${currentRunOutdir()}:${state.completed_at || currentStatus}`;
      if (completedRedirectKey !== key) {
        completedRedirectKey = key;
        $("#results-outdir").value = currentRunOutdir();
        saveSettings();
        await refreshResults();
        selectTab("results");
        saveSettings();
      }
    }
    lastRunStatus = currentStatus;
    const active = Boolean(status.active);
    if (!active && !isActiveStatus(currentStatus) && pollTimer) {
      clearInterval(pollTimer);
      pollTimer = null;
    }
    return status;
  } catch (_error) {
    if (!isActiveStatus(lastRunStatus) && pollTimer) {
      clearInterval(pollTimer);
      pollTimer = null;
    }
    return null;
  }
}

function startPolling() {
  if (pollTimer) return;
  pollTimer = setInterval(refreshAllRunState, 3000);
}

async function loadLog() {
  const options = arguments[0] || {};
  if (logInFlight) return;
  const logNode = $("#technical-log");
  if (!activeOutdir) {
    renderTerminalLog(logNode, EMPTY_LOG_TEXT);
    return;
  }
  logInFlight = true;
  const out = encodeURIComponent(currentRunOutdir());
  const statusPayload = options.statusPayload || latestStatusPayload;
  const shouldStick =
    options.forceScroll ||
    followLog ||
    Boolean(statusPayload?.active);
  try {
    const status = payloadStatus(statusPayload);
    const hasRunState = Boolean(statusPayload?.state?.status);
    const canFetchLog =
      statusPayload?.active ||
      hasRunState ||
      isActiveStatus(status) ||
      status === "complete" ||
      status === "failed" ||
      status === "dry_run";
    if (!canFetchLog) {
      renderTerminalLog(logNode, EMPTY_LOG_TEXT);
      return;
    }
    const payload = await fetchJson(`/api/logs?out=${out}&max_bytes=120000`);
    renderTerminalLog(logNode, compactTerminalLog(payload.text, statusPayload));
    if (status === "failed") updateRunFailureDetail(payload.text);
    if (shouldStick && logNode) {
      scrollLogToBottom(logNode);
    }
  } catch (error) {
    if (!options.silent) {
      renderTerminalLog(logNode, error.message);
    }
  } finally {
    logInFlight = false;
  }
}

function renderCommonPaths(paths) {
  const box = $("#common-paths");
  if (!box) return;
  box.replaceChildren();
  paths.forEach((item) => {
    const button = document.createElement("button");
    text(button, item.label);
    button.title = item.path;
    button.addEventListener("click", () => {
      $("#fastq-dir").value = item.path;
      invalidateScanReady();
      saveSettings();
      setFolderMessage("Folder selected. Click Scan folder when ready.", "ok");
    });
    box.appendChild(button);
  });
}

function resetInstalledPath(selector, appRoot, fallback) {
  const node = $(selector);
  if (!node || !appRoot || !node.value) return;
  const normalizedValue = node.value.replaceAll("\\", "/").toLowerCase();
  const normalizedRoot = appRoot.replaceAll("\\", "/").toLowerCase();
  if (normalizedValue.startsWith(normalizedRoot)) {
    node.value = fallback;
  }
}

function currentRunOutdir() {
  return activeOutdir || $("#outdir").value || "";
}

async function loadHealth() {
  const payload = await fetchJson("/api/health");
  pathStyle = payload.path_style || pathStyle;
  workspaceRoot = payload.workspace_root || workspaceRoot;
  renderCommonPaths(payload.common_paths || []);
  if (payload.workspace_root && $("#workspace-root")) {
    text($("#workspace-root"), payload.workspace_root);
  }
  resetInstalledPath("#fastq-dir", payload.app_root, "");
  resetInstalledPath("#run-samples", payload.app_root, "samples.csv");
  resetInstalledPath("#outdir", payload.app_root, "");
  resetInstalledPath("#results-outdir", payload.app_root, "");
  saveSettings();
}

function closeLibraryPicker() {
  const trigger = $("#run-library-trigger");
  const options = $("#run-library-options");
  if (!trigger || !options) return;
  options.hidden = true;
  trigger.setAttribute("aria-expanded", "false");
}

function bindLibraryPicker() {
  const picker = $("#run-library-picker");
  const trigger = $("#run-library-trigger");
  const options = $("#run-library-options");
  if (!picker || !trigger || !options) return;
  trigger.addEventListener("click", () => {
    if (trigger.disabled) return;
    const willOpen = options.hidden;
    options.hidden = !willOpen;
    trigger.setAttribute("aria-expanded", String(willOpen));
  });
  options.addEventListener("change", (event) => {
    const checkbox = event.target.closest('input[type="checkbox"][data-library]');
    if (!checkbox) return;
    const entries = libraryCheckboxes();
    if (checkbox.dataset.library === "all") {
      entries.forEach((entry) => { entry.checked = checkbox.checked; });
    }
    updateLibraryPickerSummary();
    updateRunButtonAvailability();
    saveSettings();
  });
  document.addEventListener("click", (event) => {
    if (!picker.contains(event.target)) closeLibraryPicker();
  });
}

function bindHelpControls() {
  const appButton = $(".app-help-button");
  const appPopover = $("#app-help-popover");
  const parameterPopover = document.createElement("div");
  parameterPopover.className = "parameter-help-popover";
  parameterPopover.setAttribute("role", "tooltip");
  parameterPopover.hidden = true;
  document.body.appendChild(parameterPopover);
  let activeParameterButton = null;

  const closeAppHelp = () => {
    if (!appButton || !appPopover) return;
    appPopover.hidden = true;
    appButton.setAttribute("aria-expanded", "false");
  };
  const closeParameterHelp = () => {
    parameterPopover.hidden = true;
    if (activeParameterButton) activeParameterButton.setAttribute("aria-expanded", "false");
    activeParameterButton = null;
  };
  const openParameterHelp = (button) => {
    const message = button.dataset.tooltip || "";
    if (!message) return;
    closeAppHelp();
    if (activeParameterButton === button && !parameterPopover.hidden) {
      closeParameterHelp();
      return;
    }
    activeParameterButton = button;
    text(parameterPopover, message);
    parameterPopover.hidden = false;
    button.setAttribute("aria-expanded", "true");
    const rect = button.getBoundingClientRect();
    const popoverRect = parameterPopover.getBoundingClientRect();
    const left = Math.min(window.innerWidth - popoverRect.width - 12, Math.max(12, rect.left - popoverRect.width / 2 + rect.width / 2));
    const below = rect.bottom + 8;
    const top = below + popoverRect.height <= window.innerHeight - 12
      ? below
      : Math.max(12, rect.top - popoverRect.height - 8);
    parameterPopover.style.left = `${left}px`;
    parameterPopover.style.top = `${top}px`;
  };

  if (appButton && appPopover) {
    appButton.addEventListener("click", () => {
      const willOpen = appPopover.hidden;
      closeParameterHelp();
      appPopover.hidden = !willOpen;
      appButton.setAttribute("aria-expanded", String(willOpen));
    });
  }
  document.querySelectorAll(".param-help[data-tooltip]").forEach((button) => {
    button.setAttribute("aria-expanded", "false");
    button.addEventListener("click", (event) => {
      event.stopPropagation();
      openParameterHelp(button);
    });
  });
  document.addEventListener("click", (event) => {
    if (appButton && appPopover && !appButton.contains(event.target) && !appPopover.contains(event.target)) closeAppHelp();
    if (!parameterPopover.contains(event.target) && !event.target.closest(".param-help")) closeParameterHelp();
  });
  document.addEventListener("keydown", (event) => {
    if (event.key !== "Escape") return;
    closeAppHelp();
    closeParameterHelp();
    closeLibraryPicker();
    closeReportViewer();
    closeTableViewer();
  });
}

function bindEvents() {
  bindHelpControls();
  bindLibraryPicker();
  document.querySelectorAll(".tab").forEach((tab) => {
    tab.addEventListener("click", () => {
      selectTab(tab.dataset.tab);
      saveSettings();
    });
  });
  document.querySelectorAll(".downstream-subtab").forEach((tab) => {
    tab.addEventListener("click", () => {
      selectTab(tab.dataset.downstreamTab);
      saveSettings();
    });
  });
  $("#qc-open-run-files")?.addEventListener("click", () => {
    selectTab("results");
    saveSettings();
  });
  $("#qc-continue-dinemites")?.addEventListener("click", () => {
    setAnalysisMode("off");
    updateDinemitesRunButton();
    selectTab("dinemites");
    saveSettings();
    window.scrollTo({top: 0, behavior: "smooth"});
  });
  $("#qc-continue-dcifer")?.addEventListener("click", () => {
    setAnalysisMode("off");
    updateDciferRunButton();
    selectTab("dcifer");
    saveSettings();
    window.scrollTo({top: 0, behavior: "smooth"});
  });
  $("#analysis-open-cdhit-dinemites")?.addEventListener("click", async () => {
    setAnalysisMode("summed");
    selectTab("dinemites");
    updateDinemitesRunButton();
    await loadDinemitesResults();
    saveSettings();
    window.scrollTo({top: 0, behavior: "smooth"});
  });
  $("#analysis-open-cdhit-dcifer")?.addEventListener("click", async () => {
    setAnalysisMode("summed");
    selectTab("dcifer");
    updateDciferRunButton();
    await loadDciferResults();
    saveSettings();
    window.scrollTo({top: 0, behavior: "smooth"});
  });
  const proceedRunButton = $("#proceed-run-button");
  if (proceedRunButton) {
    proceedRunButton.addEventListener("click", () => {
      if (proceedRunButton.disabled || !scanReady) return;
      syncGeneratedSampleSheetPath();
      proceedRunButton.classList.add("is-committing");
      window.setTimeout(() => {
        proceedRunButton.classList.remove("is-committing");
        selectTab("run");
        saveSettings();
        window.scrollTo({top: 0, behavior: "smooth"});
      }, 120);
    });
  }
  ["#fastq-dir", "#metadata-path", "#metadata-sheet", "#fallback-collection-year", "#fallback-collection-day"].forEach((selector) => {
    const input = $(selector);
    if (!input) return;
    input.addEventListener("input", () => {
      invalidateScanReady();
      updateSamplePathHelp();
      if (selector === "#fastq-dir") updateScanButtonAvailability();
    });
    input.addEventListener("change", () => {
      invalidateScanReady();
      updateSamplePathHelp();
      if (selector === "#fastq-dir") updateScanButtonAvailability();
    });
  });
  document.querySelectorAll("input").forEach((input) => {
    input.addEventListener("change", saveSettings);
  });
  ["#analysis-abundance-denominator"].forEach((selector) => {
    const select = $(selector);
    if (select) select.addEventListener("change", () => {
      updateAnalysisFilterSummary();
      saveSettings();
      updateAnalysisInputAvailability(latestAnalysisInputStatus || {});
    });
  });
  $("#analysis-min-abundance-pct")?.addEventListener("input", () => {
    updateAnalysisFilterSummary();
    updateAnalysisInputAvailability(latestAnalysisInputStatus || {});
  });
  updateAnalysisFilterSummary();
  ["#dinemites-analysis-mode", "#dcifer-analysis-mode"].forEach((selector) => {
    const select = $(selector);
    if (!select) return;
    select.addEventListener("change", async () => {
      setAnalysisMode(select.value);
      const dmResults = $("#dinemites-results");
      const dcResults = $("#dcifer-results");
      if (dmResults) dmResults.hidden = true;
      if (dcResults) dcResults.hidden = true;
      updateDinemitesRunButton();
      updateDciferRunButton();
      await Promise.all([loadDinemitesResults(), loadDciferResults()]);
    });
  });
  document.querySelectorAll('input[name="analysis-input-choice"]').forEach((radio) => {
    radio.addEventListener("change", () => {
      if (!radio.checked) return;
      setAnalysisMode(radio.value);
      updateDinemitesRunButton();
      updateDciferRunButton();
    });
  });
  document.querySelectorAll(".analysis-option").forEach((option) => {
    option.addEventListener("click", (event) => {
      if (event.target.closest("button, input, label, a")) return;
      const radio = option.querySelector('input[name="analysis-input-choice"]');
      if (radio && !radio.checked) radio.click();
    });
  });
  const buildAnalysisButton = $("#analysis-build-table");
  if (buildAnalysisButton) buildAnalysisButton.addEventListener("click", () => {
    setAnalysisMode("summed");
    buildAnalysisTable("cdhit98");
  });
  const buildPrimaryButton = $("#analysis-build-primary");
  if (buildPrimaryButton) buildPrimaryButton.addEventListener("click", () => {
    setAnalysisMode("off");
    buildAnalysisTable("primary");
  });
  $("#scan-button").addEventListener("click", scanFastqs);
  document.querySelectorAll("[data-scan-filter]").forEach((button) => {
    button.addEventListener("click", () => {
      scanReviewFilter = button.dataset.scanFilter || "all";
      scanReviewPage = 1;
      document.querySelectorAll("[data-scan-filter]").forEach((item) => {
        const active = item === button;
        item.classList.toggle("is-active", active);
        item.setAttribute("aria-pressed", String(active));
      });
      renderScanPreview();
    });
  });
  $("#scan-search")?.addEventListener("input", () => {
    scanReviewPage = 1;
    renderScanPreview();
  });
  $("#scan-library-filter")?.addEventListener("change", () => {
    scanReviewLibrary = $("#scan-library-filter").value || "all";
    scanReviewPage = 1;
    renderScanPreview();
    saveSettings();
  });
  $("#scan-page-prev")?.addEventListener("click", () => {
    scanReviewPage = Math.max(1, scanReviewPage - 1);
    renderScanPreview();
  });
  $("#scan-page-next")?.addEventListener("click", () => {
    scanReviewPage += 1;
    renderScanPreview();
  });
  const checkButton = $("#check-button");
  if (checkButton) checkButton.addEventListener("click", runCheck);
  $("#run-button").addEventListener("click", startRun);
  $("#stop-button").addEventListener("click", stopRun);
  const refreshResultsButton = $("#refresh-results");
  if (refreshResultsButton) {
    refreshResultsButton.addEventListener("click", () => {
      withButtonFeedback(refreshResultsButton, "Refreshing...", async () => {
        if (!$("#results-outdir").value) $("#results-outdir").value = currentRunOutdir();
        saveSettings();
        await refreshResults();
      });
    });
  }
  const bundleButton = $("#download-bundle");
  if (bundleButton) {
    bundleButton.addEventListener("click", () => {
      const url = bundleButton.dataset.url;
      if (url) downloadFile(url, bundleButton);
    });
  }
  const openResultsButton = $("#open-results-folder");
  if (openResultsButton) openResultsButton.addEventListener("click", openSelectedResultsFolder);
  const chooseResultsButton = $("#choose-results-folder");
  if (chooseResultsButton) chooseResultsButton.addEventListener("click", chooseResultsFolder);
  const logNode = $("#technical-log");
  if (logNode) {
    logNode.addEventListener("scroll", () => {
      followLog = logNode.scrollTop + logNode.clientHeight >= logNode.scrollHeight - 36;
    });
  }
  $("#show-run-log")?.addEventListener("click", () => {
    const details = $("#technical-log-details");
    if (!details) return;
    details.open = true;
    details.scrollIntoView({behavior: "smooth", block: "start"});
  });
  $("#browse-button").addEventListener("click", chooseFastqFolder);
  $("#choose-metadata-file").addEventListener("click", chooseMetadataFile);
  $("#choose-kelt-barcode-map")?.addEventListener("click", chooseKeltBarcodeMap);
  $("#kelt-barcode-map")?.addEventListener("change", () => {
    saveSettings();
    inspectKeltBarcodeMap();
  });
  $("#metadata-sheet")?.addEventListener("change", () => {
    metadataContract.columns = {};
    metadataContract.detection_value_map = {};
    invalidateScanReady();
    inspectMetadata();
  });
  $("#metadata-path")?.addEventListener("change", () => {
    resetMetadataContract();
    invalidateScanReady();
    inspectMetadata();
  });
  document.querySelectorAll("[data-metadata-column]").forEach((select) => {
    select.addEventListener("change", () => {
      metadataContract.columns[select.dataset.metadataColumn] = select.value;
      if (select.dataset.metadataColumn === "metadata_pcr") {
        metadataContract.detection_value_map = {};
      }
      if (select.dataset.metadataColumn === "metadata_status") {
        metadataContract.excluded_status_values = [];
      }
      invalidateScanReady();
      inspectMetadata();
    });
  });
  $("#metadata-date-order")?.addEventListener("change", () => {
    invalidateScanReady();
    saveSettings();
    inspectMetadata();
  });
  $("#choose-outdir-button").addEventListener("click", chooseOutputFolder);
  $("#dry-run").addEventListener("change", () => {
    $("#dry-run").checked = false;
    text($("#run-button"), "Start run");
    saveSettings();
  });
  $("#outdir").addEventListener("input", updateRunButtonAvailability);
  $("#outdir").addEventListener("change", () => {
    activeOutdir = "";
    syncGeneratedSampleSheetPath();
    updateSamplePathHelp();
    saveSettings();
    refreshAllRunState();
    updateRunButtonAvailability();
  });
  $("#run-samples").addEventListener("input", updateSamplePathHelp);
  $("#run-samples").addEventListener("change", () => {
    updateSamplePathHelp();
    saveSettings();
  });
}


// ---------------------------------------------------------------------------
// DINEMITES analysis
// ---------------------------------------------------------------------------

let dinemitesPollTimer = null;
let dinemitesPlotItems = [];
let dinemitesPlotIndex = 0;
let dinemitesPlotSwapToken = 0;
let dinemitesPlotZoom = 1.75;

function setDinemitesResultsView(name) {
  const viewName = name === "tables" ? "tables" : "plots";
  document.querySelectorAll("[data-dinemites-results-view]").forEach((button) => {
    const active = button.dataset.dinemitesResultsView === viewName;
    button.classList.toggle("is-active", active);
    button.setAttribute("aria-selected", active ? "true" : "false");
  });
  const plots = $("#dinemites-results-plots");
  const tables = $("#dinemites-results-tables");
  if (plots) plots.hidden = viewName !== "plots";
  if (tables) tables.hidden = viewName !== "tables";
}

function setDinemitesPlotZoom(value) {
  dinemitesPlotZoom = Math.max(1, Math.min(2.5, Number(value) || 1));
  const image = document.querySelector(".dinemites-plot-viewport img");
  const label = $("#dm-plot-zoom-label");
  if (image) image.style.width = `${Math.round(dinemitesPlotZoom * 100)}%`;
  if (label) text(label, `${Math.round(dinemitesPlotZoom * 100)}%`);
  const zoomOut = $("#dm-plot-zoom-out");
  const zoomIn = $("#dm-plot-zoom-in");
  if (zoomOut) zoomOut.disabled = dinemitesPlotZoom <= 1;
  if (zoomIn) zoomIn.disabled = dinemitesPlotZoom >= 2.5;
}

function dinemitesOutdir() {
  return activeOutdir || $("#results-outdir")?.value || "results";
}

function hasSelectedRunOutputs() {
  const outdir = activeOutdir || $("#results-outdir")?.value || "";
  if (!outdir) return false;
  if (latestResultsPayload && latestResultsPayload.outdir === outdir) {
    return Number(latestResultsPayload.ready_counts?.core || 0) > 0;
  }
  return true;
}

function selectedAnalysisInputState() {
  const mode = analysisApiMode();
  const input = latestAnalysisInputStatus?.[mode] || {};
  return {
    mode,
    ready: analysisInputMatchesSettings(input),
    missingMessage: mode === "cdhit98"
      ? "Prepare the CD-HIT 98.9% input in Quality control first."
      : "Prepare the exact-allele input in Quality control first."
  };
}

function analysisResultLabel(mode) {
  if (mode === "cdhit98") return "Sensitivity result - CD-HIT 98.9% locus-specific clusters";
  return "Primary result - exact CIGAR alleles";
}

function updateDinemitesRunButton() {
  const btn = $("#dinemites-run");
  if (!btn) return;
  const selected = selectedAnalysisInputState();
  const ready = hasSelectedRunOutputs()
    && selected.ready
    && latestDinemitesReadiness?.dinemites_ready !== false;
  btn.disabled = !ready;
  btn.title = ready
    ? ""
    : (!selected.ready
      ? selected.missingMessage
      : (hasSelectedRunOutputs() ? "DINEMITES needs complete dates and repeated participant visits." : "Run the main pipeline first."));
}

function readinessItem(label, value, state = "ok") {
  const item = document.createElement("div");
  item.className = `readiness-item is-${state}`;
  const marker = document.createElement("span");
  marker.className = "readiness-marker";
  marker.setAttribute("aria-hidden", "true");
  const name = document.createElement("span");
  text(name, label);
  const detail = document.createElement("strong");
  text(detail, value);
  item.append(marker, name, detail);
  return item;
}

function renderDinemitesReadiness(readiness) {
  latestDinemitesReadiness = readiness || null;
  const container = $("#dinemites-readiness");
  if (!container) return;
  container.replaceChildren();
  if (!readiness || !readiness.samples) {
    container.appendChild(readinessItem("Completed run", "Not available", "bad"));
    updateDinemitesRunButton();
    return;
  }
  container.append(
    readinessItem("Collection dates", readiness.missing_dates ? `${readiness.missing_dates} missing` : `${readiness.samples} complete`, readiness.missing_dates ? "bad" : "ok"),
    readinessItem("Repeated participants", `${readiness.repeated_subjects} of ${readiness.participants}`, readiness.repeated_subjects ? "ok" : "bad"),
    readinessItem("Detected loci", String(readiness.loci), readiness.loci ? "ok" : "bad"),
    readinessItem("Stan runtime", readiness.cmdstan_available ? "Available" : "Checked when model starts", readiness.cmdstan_available ? "ok" : "neutral")
  );
  updateDinemitesRunButton();
}

function dinemitesTLagValue() {
  const noCutoff = $("#dinemites-no-day-cutoff")?.checked;
  if (noCutoff) return "Inf";
  return $("#dinemites-t-lag")?.value || "Inf";
}

function selectedDinemitesCovariates() {
  const covariates = [];
  if ($("#dinemites-covariate-season")?.checked) {
    covariates.push("covariate_season", "covariate_season_missing");
  }
  if ($("#dinemites-covariate-age")?.checked) {
    covariates.push("covariate_age", "covariate_age_missing");
  }
  if ($("#dinemites-covariate-gender")?.checked) {
    covariates.push("covariate_gender", "covariate_gender_missing");
  }
  const custom = $("#dinemites-infection-covariates")?.value || "";
  custom.split(",").map((item) => item.trim()).filter(Boolean).forEach((item) => {
    if (!covariates.includes(item)) covariates.push(item);
  });
  return covariates.length > 0 ? covariates.join(",") : "none";
}

function updateDinemitesModelSettingsVisibility() {
  const model = $("#dinemites-model")?.value || "simple";
  const simpleSettings = $("#dinemites-simple-settings");
  const bayesianSettings = $("#dinemites-bayesian-settings");
  if (simpleSettings) simpleSettings.hidden = model !== "simple";
  if (bayesianSettings) bayesianSettings.hidden = model !== "bayesian";
  document.querySelectorAll("[data-dinemites-model]").forEach((button) => {
    const selected = button.dataset.dinemitesModel === model;
    button.classList.toggle("is-active", selected);
    button.setAttribute("aria-pressed", selected ? "true" : "false");
  });
  const notes = {
    simple: "Uses recent visits to classify persistence. Fast and deterministic; useful as a transparent first pass.",
    clustering: "Learns persistence patterns without individual covariates. The first run may compile a cached Stan model.",
    bayesian: "Models persistence and uncertainty with optional season covariates. Compare it with the other models when conclusions differ."
  };
  text($("#dinemites-model-note"), notes[model]);
}

function handleDinemitesDayCutoffToggle() {
  const noCutoff = $("#dinemites-no-day-cutoff");
  const tLag = $("#dinemites-t-lag");
  if (tLag && noCutoff) {
    tLag.disabled = noCutoff.checked;
  }
  saveSettings();
}

function handleDinemitesToggle() {
  const controls = $("#dinemites-controls");
  const pill = $("#dinemites-status");
  if ($("#dinemites-enable")) $("#dinemites-enable").checked = true;
  if (controls) controls.hidden = false;
  updateDinemitesModelSettingsVisibility();
  setPill(pill, "Ready", "ok");
  updateDinemitesRunButton();
  saveSettings();
}

async function runDinemites() {
  const btn = $("#dinemites-run");
  const msg = $("#dinemites-message");
  btn.disabled = true;
  setAnalysisModeControlsDisabled(true);
  msg.className = "inline-message";
  text(msg, "Starting DINEMITES analysis...");
  try {
    const payload = await postJson("/api/dinemites/run", {
      model_type: $("#dinemites-model").value,
      outdir: dinemitesOutdir(),
      samples: $("#run-samples").value,
      n_lags: Number($("#dinemites-n-lags").value || 3),
      t_lag: dinemitesTLagValue(),
      min_abundance_pct: analysisMinAbundancePct(),
      abundance_denominator: analysisAbundanceDenominator(),
      analysis_mode: analysisApiMode(),
      no_day_cutoff: $("#dinemites-no-day-cutoff").checked,
      seed: Number($("#dinemites-seed").value || 1),
      refresh: Number($("#dinemites-refresh-interval").value || 100),
      bayesian_lag_days: Number($("#dinemites-bayesian-lag-days").value || 30),
      bayesian_chains: Number($("#dinemites-bayesian-chains").value || 4),
      bayesian_parallel_chains: Number($("#dinemites-bayesian-parallel-chains").value || 2),
      bayesian_iter_warmup: Number($("#dinemites-bayesian-warmup").value || 500),
      bayesian_iter_sampling: Number($("#dinemites-bayesian-sampling").value || 500),
      bayesian_adapt_delta: Number($("#dinemites-bayesian-adapt-delta").value || 0.99),
      infection_general_covariates: selectedDinemitesCovariates(),
      bayesian_drop_out: $("#dinemites-bayesian-drop-out").checked
    });
    text(msg, "DINEMITES is running. Results will appear automatically when complete.");
    setPill($("#dinemites-status"), "Running", "warn");
    startDinemitesPolling();
  } catch (error) {
    text(msg, userMessage(error.message));
    msg.classList.add("bad");
    setPill($("#dinemites-status"), "Failed", "bad");
    btn.disabled = false;
    setAnalysisModeControlsDisabled(false);
  }
}

function startDinemitesPolling() {
  if (dinemitesPollTimer) return;
  dinemitesPollTimer = setInterval(pollDinemitesStatus, 3000);
}

function stopDinemitesPolling() {
  if (dinemitesPollTimer) {
    clearInterval(dinemitesPollTimer);
    dinemitesPollTimer = null;
  }
}

async function pollDinemitesStatus() {
  const out = encodeURIComponent(dinemitesOutdir());
  const mode = encodeURIComponent(analysisApiMode());
  try {
    const payload = await fetchJson(`/api/dinemites/status?out=${out}&mode=${mode}`);
    const status = payload.status || "idle";
    const msg = $("#dinemites-message");
    if (status === "running") {
      setPill($("#dinemites-status"), "Running", "warn");
      const btn = $("#dinemites-run");
      if (btn) btn.disabled = true;
      if (msg && !msg.textContent.trim()) {
        msg.className = "inline-message";
        text(msg, "DINEMITES is running. Results will appear automatically when complete.");
      }
    } else if (status === "complete") {
      setPill($("#dinemites-status"), "Complete", "ok");
      stopDinemitesPolling();
      await loadDinemitesResults();
      if (msg) {
        msg.className = "inline-message ok";
        text(msg, "DINEMITES complete. Results are shown below.");
      }
      updateDinemitesRunButton();
      setAnalysisModeControlsDisabled(false);
    } else if (status === "failed") {
      setPill($("#dinemites-status"), "Failed", "bad");
      const detail = payload.state?.detail || "DINEMITES analysis failed.";
      msg.className = "inline-message bad";
      text(msg, userMessage(detail));
      stopDinemitesPolling();
      updateDinemitesRunButton();
      setAnalysisModeControlsDisabled(false);
    } else {
      stopDinemitesPolling();
      updateDinemitesRunButton();
      setAnalysisModeControlsDisabled(false);
    }
  } catch (_error) {
    // Silently retry on transient errors.
  }
}

function renderDinemitesResultFiles(files = {}) {
  const container = $("#dm-output-tables");
  const count = $("#dm-output-table-count");
  if (!container) return;
  const definitions = [
    ["input", "Model input"],
    ["allele_probabilities", "Allele probabilities"],
    ["allele_key", "Allele key"],
    ["molfoi", "Molecular force of infection"],
    ["new_infections", "New infection events"],
  ];
  const items = definitions
    .filter(([key]) => files[key]?.exists)
    .map(([key, label]) => ({label, ...files[key]}));
  text(count, items.length ? `${items.length} table${items.length === 1 ? "" : "s"}` : "");
  container.replaceChildren(resultTable(items, "DINEMITES tables will appear after the analysis completes."));
}

function renderDinemitesPlots(plots) {
  const gallery = $("#dinemites-plot-gallery");
  if (!gallery) return;
  const items = Array.isArray(plots) ? plots.filter((plot) => plot && plot.exists && plot.view_url) : [];
  const previousFilename = dinemitesPlotItems[dinemitesPlotIndex]?.filename || $("#dm-plot-selector")?.value || "";
  dinemitesPlotItems = items;
  dinemitesPlotIndex = Math.max(0, items.findIndex((plot) => plot.filename === previousFilename));
  if (dinemitesPlotIndex < 0) dinemitesPlotIndex = 0;

  updatePlotJump(
    gallery,
    $("#dm-plot-count"),
    null,
    items.length,
    "No DINEMITES plots available yet.",
    "DINEMITES plot",
    "DINEMITES plots"
  );

  updateDinemitesPlotBrowser();
  if (!items.length) {
    gallery.replaceChildren();
    gallery.hidden = true;
    return;
  }

  gallery.hidden = false;
  renderSelectedDinemitesPlot();
}

function dinemitesPlotLabel(plot) {
  return plot?.subject || plot?.filename || "Subject";
}

function updateDinemitesPlotBrowser() {
  const selector = $("#dm-plot-selector");
  const prev = $("#dm-prev-plot");
  const next = $("#dm-next-plot");
  const hasPlots = dinemitesPlotItems.length > 0;

  if (selector) {
    selector.replaceChildren();
    if (!hasPlots) {
      const option = document.createElement("option");
      option.value = "";
      text(option, "No plots");
      selector.appendChild(option);
    } else {
      dinemitesPlotItems.forEach((plot, index) => {
        const option = document.createElement("option");
        option.value = plot.filename || String(index);
        text(option, dinemitesPlotLabel(plot));
        selector.appendChild(option);
      });
      selector.selectedIndex = dinemitesPlotIndex;
    }
    selector.disabled = !hasPlots;
  }

  if (prev) prev.disabled = !hasPlots || dinemitesPlotIndex <= 0;
  if (next) next.disabled = !hasPlots || dinemitesPlotIndex >= dinemitesPlotItems.length - 1;
}

function setDinemitesPlotIndex(index) {
  if (!dinemitesPlotItems.length) return;
  const gallery = $("#dinemites-plot-gallery");
  const scrollPosition = window.scrollY;
  const previousHeight = gallery?.getBoundingClientRect().height || 0;
  if (gallery && previousHeight > 0) gallery.style.minHeight = `${Math.ceil(previousHeight)}px`;
  const swapToken = ++dinemitesPlotSwapToken;
  dinemitesPlotIndex = Math.min(dinemitesPlotItems.length - 1, Math.max(0, index));
  updateDinemitesPlotBrowser();
  const image = renderSelectedDinemitesPlot();
  window.scrollTo({top: scrollPosition, left: window.scrollX, behavior: "auto"});

  const finishSwap = () => {
    if (swapToken !== dinemitesPlotSwapToken) return;
    if (gallery) gallery.style.minHeight = "";
    window.scrollTo({top: scrollPosition, left: window.scrollX, behavior: "auto"});
  };
  if (image?.complete) {
    requestAnimationFrame(finishSwap);
  } else if (image) {
    image.addEventListener("load", finishSwap, {once: true});
    image.addEventListener("error", finishSwap, {once: true});
  } else {
    finishSwap();
  }
}

function renderSelectedDinemitesPlot() {
  const gallery = $("#dinemites-plot-gallery");
  if (!gallery) return;
  gallery.replaceChildren();
  const plot = dinemitesPlotItems[dinemitesPlotIndex];
  if (!plot) {
    gallery.hidden = true;
    return;
  }
  gallery.hidden = false;

  const figure = document.createElement("figure");
  figure.className = "dinemites-plot-card";

  const toolbar = document.createElement("div");
  toolbar.className = "dinemites-plot-toolbar";
  const toolbarTitle = document.createElement("strong");
  text(toolbarTitle, dinemitesPlotLabel(plot));
  toolbar.appendChild(toolbarTitle);

  const zoomControls = document.createElement("div");
  zoomControls.className = "dinemites-zoom-controls";
  const zoomOut = document.createElement("button");
  zoomOut.id = "dm-plot-zoom-out";
  zoomOut.type = "button";
  zoomOut.title = "Zoom out";
  zoomOut.setAttribute("aria-label", "Zoom DINEMITES plot out");
  text(zoomOut, "−");
  zoomOut.addEventListener("click", () => setDinemitesPlotZoom(dinemitesPlotZoom - 0.25));
  const zoomLabel = document.createElement("span");
  zoomLabel.id = "dm-plot-zoom-label";
  zoomLabel.className = "dinemites-zoom-label";
  const zoomIn = document.createElement("button");
  zoomIn.id = "dm-plot-zoom-in";
  zoomIn.type = "button";
  zoomIn.title = "Zoom in";
  zoomIn.setAttribute("aria-label", "Zoom DINEMITES plot in");
  text(zoomIn, "+");
  zoomIn.addEventListener("click", () => setDinemitesPlotZoom(dinemitesPlotZoom + 0.25));
  const fit = document.createElement("button");
  fit.type = "button";
  fit.className = "dinemites-zoom-fit";
  fit.title = "Fit plot to available width";
  text(fit, "Fit");
  fit.addEventListener("click", () => setDinemitesPlotZoom(1));
  zoomControls.append(zoomOut, zoomLabel, zoomIn, fit);
  toolbar.appendChild(zoomControls);
  figure.appendChild(toolbar);

  const img = document.createElement("img");
  img.src = plot.view_url;
  img.alt = `DINEMITES plot for ${dinemitesPlotLabel(plot)}`;
  img.loading = "eager";
  img.decoding = "async";
  if (Number(plot.width) > 0 && Number(plot.height) > 0) {
    img.width = Number(plot.width);
    img.height = Number(plot.height);
  }
  const viewport = document.createElement("div");
  viewport.className = "dinemites-plot-viewport";
  viewport.appendChild(img);
  figure.appendChild(viewport);
  requestAnimationFrame(() => setDinemitesPlotZoom(dinemitesPlotZoom));

  const caption = document.createElement("figcaption");
  const captionNote = document.createElement("span");
  text(captionNote, "Use zoom for labels or open the original image for detailed inspection.");
  caption.appendChild(captionNote);

  if (plot.download_url) {
    const actions = document.createElement("div");
    actions.className = "dinemites-plot-actions";
    const viewButton = document.createElement("button");
    viewButton.type = "button";
    text(viewButton, "View full resolution");
    viewButton.addEventListener("click", () => {
      openImageViewer(plot.view_url, `DINEMITES plot - ${dinemitesPlotLabel(plot)}`);
    });
    actions.appendChild(viewButton);

    const button = document.createElement("button");
    button.type = "button";
    text(button, "Download plot");
    button.addEventListener("click", () => {
      window.location.href = plot.download_url;
    });
    actions.appendChild(button);
    caption.appendChild(actions);
  }

  figure.appendChild(caption);
  gallery.appendChild(figure);
  return img;
}

async function loadDinemitesResults() {
  const out = encodeURIComponent(dinemitesOutdir());
  const mode = encodeURIComponent(analysisApiMode());
  try {
    const payload = await fetchJson(`/api/dinemites/results?out=${out}&mode=${mode}`);
    const state = payload.state || {};
    renderDinemitesReadiness(payload.readiness || null);
    const status = state.status || "idle";
    const resultsPanel = $("#dinemites-results");

    if (status === "complete") {
      if (resultsPanel) resultsPanel.hidden = false;
      setPill($("#dinemites-status"), "Complete", "ok");
      text(
        $("#dm-results-input-mode"),
        analysisResultLabel(payload.analysis_mode)
      );

      text($("#dm-new-infections"), formatNumber(payload.summary?.new_infections, 0));
      text($("#dm-molfoi"), formatNumber(payload.summary?.molfoi, 2));
      text($("#dm-subjects"), formatNumber(payload.summary?.subjects, 0));
      text($("#dm-model"), String(payload.model_summary?.model_type || "--"));

      const diagnosticsNode = $("#dm-diagnostics");
      const diagnostics = payload.model_summary?.diagnostics || {};
      const diagnosticMessages = [];
      if (Number(diagnostics.divergent_transitions || 0) > 0) diagnosticMessages.push(`${diagnostics.divergent_transitions} divergent transition(s)`);
      if (Number(diagnostics.max_treedepth_transitions || 0) > 0) diagnosticMessages.push(`${diagnostics.max_treedepth_transitions} transition(s) reached maximum tree depth`);
      if (Number(diagnostics.max_rhat || 0) > 1.05) diagnosticMessages.push(`maximum R-hat ${formatNumber(diagnostics.max_rhat, 3)}`);
      if (diagnosticsNode) {
        diagnosticsNode.hidden = !diagnosticMessages.length;
        text(diagnosticsNode, diagnosticMessages.length ? `Bayesian diagnostics need review: ${diagnosticMessages.join("; ")}.` : "");
      }
      renderDinemitesPlots(payload.plots || []);
      renderDinemitesResultFiles(payload.files || {});
    } else if (status === "running") {
      setPill($("#dinemites-status"), "Running", "warn");
      setAnalysisModeControlsDisabled(true);
      const runBtn = $("#dinemites-run");
      const msg = $("#dinemites-message");
      if (runBtn) runBtn.disabled = true;
      if (msg && !msg.textContent.trim()) {
        msg.className = "inline-message";
        text(msg, "DINEMITES is running. Results will appear automatically when complete.");
      }
      startDinemitesPolling();
    } else if (status === "failed") {
      if (resultsPanel) resultsPanel.hidden = true;
      renderDinemitesPlots([]);
      renderDinemitesResultFiles({});
      setPill($("#dinemites-status"), "Failed", "bad");
      setAnalysisModeControlsDisabled(false);
    }
  } catch (_error) {
    // Results may not exist yet.
    renderDinemitesPlots([]);
    renderDinemitesResultFiles({});
  }
}

function bindDinemitesEvents() {
  const toggle = $("#dinemites-enable");
  if (toggle) toggle.addEventListener("change", handleDinemitesToggle);
  const runBtn = $("#dinemites-run");
  if (runBtn) runBtn.addEventListener("click", runDinemites);
  const plotSelector = $("#dm-plot-selector");
  if (plotSelector) {
    plotSelector.addEventListener("change", () => {
      setDinemitesPlotIndex(plotSelector.selectedIndex);
    });
  }
  const prevPlot = $("#dm-prev-plot");
  if (prevPlot) prevPlot.addEventListener("click", () => setDinemitesPlotIndex(dinemitesPlotIndex - 1));
  const nextPlot = $("#dm-next-plot");
  if (nextPlot) nextPlot.addEventListener("click", () => setDinemitesPlotIndex(dinemitesPlotIndex + 1));
  document.querySelectorAll("[data-dinemites-results-view]").forEach((button) => {
    button.addEventListener("click", () => setDinemitesResultsView(button.dataset.dinemitesResultsView));
  });
  const modelSelect = $("#dinemites-model");
  if (modelSelect) {
    modelSelect.addEventListener("change", () => {
      updateDinemitesModelSettingsVisibility();
      saveSettings();
    });
  }
  document.querySelectorAll("[data-dinemites-model]").forEach((button) => {
    button.addEventListener("click", () => {
      if (!modelSelect) return;
      modelSelect.value = button.dataset.dinemitesModel;
      modelSelect.dispatchEvent(new Event("change", { bubbles: true }));
    });
  });
  const nLagsInput = $("#dinemites-n-lags");
  if (nLagsInput) nLagsInput.addEventListener("change", saveSettings);
  const tLagInput = $("#dinemites-t-lag");
  if (tLagInput) tLagInput.addEventListener("change", saveSettings);
  const noDayCutoff = $("#dinemites-no-day-cutoff");
  if (noDayCutoff) noDayCutoff.addEventListener("change", handleDinemitesDayCutoffToggle);
  const seedInput = $("#dinemites-seed");
  if (seedInput) seedInput.addEventListener("change", saveSettings);
  const refreshInput = $("#dinemites-refresh-interval");
  if (refreshInput) refreshInput.addEventListener("change", saveSettings);
  [
    "#dinemites-bayesian-lag-days",
    "#dinemites-bayesian-chains",
    "#dinemites-bayesian-parallel-chains",
    "#dinemites-bayesian-warmup",
    "#dinemites-bayesian-sampling",
    "#dinemites-bayesian-adapt-delta",
    "#dinemites-covariate-season",
    "#dinemites-covariate-age",
    "#dinemites-covariate-gender",
    "#dinemites-infection-covariates",
    "#dinemites-bayesian-drop-out"
  ].forEach((selector) => {
    const control = $(selector);
    if (control) control.addEventListener("change", saveSettings);
  });
}

// ---------------------------------------------------------------------------
// dcifer analysis
// ---------------------------------------------------------------------------

let dciferPollTimer = null;

function dciferOutdir() {
  return activeOutdir || $("#results-outdir")?.value || "results";
}

function updateDciferRunButton() {
  const btn = $("#dcifer-run");
  if (!btn) return;
  const selected = selectedAnalysisInputState();
  const ready = hasSelectedRunOutputs() && selected.ready;
  btn.disabled = !ready;
  btn.title = ready
    ? ""
    : (!selected.ready
      ? selected.missingMessage
      : "Run the main pipeline first.");
}

function renderDciferReadiness(readiness) {
  latestDciferReadiness = readiness || null;
  const container = $("#dcifer-readiness");
  if (!container) return;
  container.replaceChildren();
  if (!readiness || !readiness.samples) {
    container.appendChild(readinessItem("Completed run", "Not available", "bad"));
    return;
  }
  container.append(
    readinessItem("Samples", String(readiness.samples), readiness.samples >= 2 ? "ok" : "bad"),
    readinessItem("Detected loci", String(readiness.loci), readiness.loci >= 2 ? "ok" : (readiness.loci === 1 ? "warn" : "bad")),
    readinessItem("Interpretation", readiness.dcifer_single_locus ? "Exploratory: single locus" : "Multi-locus", readiness.dcifer_single_locus ? "warn" : "ok")
  );
}

function handleDciferToggle() {
  const toggle = $("#dcifer-enable");
  if (toggle) toggle.checked = true;
  const controls = $("#dcifer-controls");
  const pill = $("#dcifer-status");
  if (controls) controls.hidden = false;
  setPill(pill, "Ready", "ok");
  updateDciferRunButton();
  saveSettings();
}

async function runDcifer() {
  const btn = $("#dcifer-run");
  const msg = $("#dcifer-message");
  if (btn) btn.disabled = true;
  setAnalysisModeControlsDisabled(true);
  if (msg) {
    msg.className = "inline-message";
    text(msg, "Starting Dcifer analysis...");
  }
  try {
    await postJson("/api/dcifer/run", {
      outdir: dciferOutdir(),
      samples: $("#run-samples").value,
      min_abundance_pct: analysisMinAbundancePct(),
      abundance_denominator: analysisAbundanceDenominator(),
      analysis_mode: analysisApiMode(),
      coi_lrank: Number($("#dcifer-coi-lrank").value || 2),
      ibd_grid_nr: Number($("#dcifer-ibd-grid-nr").value || 1000),
      alpha: Number($("#dcifer-alpha").value || 0.05),
      afreq_mode: "current_run"
    });
    if (msg) {
      text(msg, "Dcifer is running. Results will appear automatically when complete.");
    }
    setPill($("#dcifer-status"), "Running", "warn");
    startDciferPolling();
  } catch (error) {
    if (msg) {
      text(msg, userMessage(error.message));
      msg.classList.add("bad");
    }
    setPill($("#dcifer-status"), "Failed", "bad");
    updateDciferRunButton();
    setAnalysisModeControlsDisabled(false);
  }
}

function startDciferPolling() {
  if (dciferPollTimer) return;
  dciferPollTimer = setInterval(pollDciferStatus, 3000);
}

function stopDciferPolling() {
  if (dciferPollTimer) {
    clearInterval(dciferPollTimer);
    dciferPollTimer = null;
  }
}

async function pollDciferStatus() {
  const out = encodeURIComponent(dciferOutdir());
  const mode = encodeURIComponent(analysisApiMode());
  try {
    const payload = await fetchJson(`/api/dcifer/status?out=${out}&mode=${mode}`);
    const status = payload.status || "idle";
    const msg = $("#dcifer-message");
    if (status === "running") {
      setPill($("#dcifer-status"), "Running", "warn");
      const btn = $("#dcifer-run");
      if (btn) btn.disabled = true;
      if (msg && !msg.textContent.trim()) {
        msg.className = "inline-message";
        text(msg, "Dcifer is running. Results will appear automatically when complete.");
      }
    } else if (status === "complete") {
      setPill($("#dcifer-status"), "Complete", "ok");
      stopDciferPolling();
      await loadDciferResults();
      if (msg) {
        msg.className = "inline-message ok";
        text(msg, "Dcifer complete. Results are shown below.");
      }
      updateDciferRunButton();
      setAnalysisModeControlsDisabled(false);
    } else if (status === "failed") {
      setPill($("#dcifer-status"), "Failed", "bad");
      const detail = payload.state?.detail || "Dcifer analysis failed.";
      if (msg) {
        msg.className = "inline-message bad";
        text(msg, userMessage(detail));
      }
      stopDciferPolling();
      updateDciferRunButton();
      setAnalysisModeControlsDisabled(false);
    } else {
      stopDciferPolling();
      updateDciferRunButton();
      setAnalysisModeControlsDisabled(false);
    }
  } catch (_error) {
    // Leave polling active for transient local server hiccups.
  }
}

function matrixHasValues(matrix) {
  return Boolean(
    matrix &&
    Array.isArray(matrix.labels) &&
    matrix.labels.length &&
    Array.isArray(matrix.rows) &&
    matrix.rows.length
  );
}

function colorFromStops(stops, value) {
  const clamped = Math.max(0, Math.min(1, value));
  const scaled = clamped * (stops.length - 1);
  const index = Math.min(stops.length - 2, Math.floor(scaled));
  const local = scaled - index;
  const left = stops[index];
  const right = stops[index + 1];
  const mix = left.map((channel, offset) => Math.round(channel + (right[offset] - channel) * local));
  return `rgb(${mix[0]}, ${mix[1]}, ${mix[2]})`;
}

function dciferHeatmapColor(kind, value) {
  if (value === null || value === undefined || Number.isNaN(Number(value))) return "#dedede";
  const number = Number(value);
  if (kind === "pvalue") {
    const intensity = 1 - Math.max(0, Math.min(1, number / 0.5));
    return colorFromStops([
      [247, 247, 252],
      [218, 218, 235],
      [158, 154, 200],
      [106, 81, 163],
      [63, 0, 125]
    ], intensity);
  }
  return colorFromStops([
    [247, 251, 255],
    [198, 219, 239],
    [107, 174, 214],
    [33, 113, 181],
    [8, 48, 107]
  ], Math.max(0, Math.min(1, number)));
}

function matrixValueText(kind, value) {
  if (value === null || value === undefined || Number.isNaN(Number(value))) return "NA";
  return kind === "pvalue" ? formatPValue(value) : formatNumber(value, 3);
}

function dciferPlotForKind(plots, kind) {
  const needle = kind === "pvalue" ? "pvalue" : "relatedness";
  return (Array.isArray(plots) ? plots : []).find((plot) => {
    const name = `${plot?.filename || ""} ${plot?.title || ""}`.toLowerCase();
    return name.includes(needle);
  });
}

function dciferSampleDisplayMap(labels, rows) {
  const ordered = [];
  const seen = new Set();
  const addSample = (value) => {
    const sampleId = String(value || "").trim();
    if (!sampleId || seen.has(sampleId)) return;
    seen.add(sampleId);
    ordered.push(sampleId);
  };
  labels.forEach(addSample);
  rows.forEach((rowItem) => addSample(rowItem?.sample_id));
  const digits = Math.max(2, String(ordered.length || 1).length);
  const entries = ordered.map((sampleId, index) => ({
    sampleId,
    displayId: `S${String(index + 1).padStart(digits, "0")}`,
  }));
  return {
    entries,
    bySample: new Map(entries.map((entry) => [entry.sampleId, entry.displayId])),
  };
}

function dciferDisplayId(labelMap, sampleId, fallbackIndex) {
  const full = String(sampleId || "").trim();
  if (labelMap.bySample.has(full)) return labelMap.bySample.get(full);
  const digits = Math.max(2, String(Math.max(labelMap.entries.length, fallbackIndex + 1)).length);
  return `S${String(fallbackIndex + 1).padStart(digits, "0")}`;
}

function renderDciferHeatmapKey(entries) {
  const details = document.createElement("details");
  details.className = "dcifer-heatmap-key";

  const summary = document.createElement("summary");
  text(summary, "Sample key");
  details.appendChild(summary);

  const list = document.createElement("div");
  list.className = "dcifer-heatmap-key-list";
  entries.forEach((entry) => {
    const item = document.createElement("div");
    item.className = "dcifer-heatmap-key-item";

    const displayId = document.createElement("strong");
    text(displayId, entry.displayId);
    item.appendChild(displayId);

    const sample = document.createElement("span");
    text(sample, entry.sampleId);
    sample.title = entry.sampleId;
    item.appendChild(sample);

    list.appendChild(item);
  });
  details.appendChild(list);
  return details;
}

function renderDciferHeatmapCard(kind, matrix, plotInfo) {
  const namespace = "http://www.w3.org/2000/svg";
  const labels = matrix.labels || [];
  const rows = matrix.rows || [];
  const labelMap = dciferSampleDisplayMap(labels, rows);
  const dimension = Math.max(labels.length, rows.length);
  const plotSpan = dimension > 50 ? 900 : dimension > 30 ? 820 : 720;
  const left = 70;
  const top = 34;
  const right = 18;
  const bottom = 28;
  const cell = plotSpan / Math.max(1, dimension);
  const labelStep = Math.max(1, Math.ceil(dimension / 14));
  const labelFont = dimension > 45 ? 12 : 13;
  const width = left + labels.length * cell + right;
  const height = top + rows.length * cell + bottom;
  const titleText = kind === "pvalue"
    ? "p-value heatmap"
    : "Relatedness heatmap";

  const figure = document.createElement("figure");
  figure.className = "dcifer-heatmap-card";

  const header = document.createElement("div");
  header.className = "dcifer-heatmap-header";
  const heading = document.createElement("h3");
  text(heading, titleText);
  header.appendChild(heading);
  const headerNote = document.createElement("p");
  const previewNote = matrix.truncated
    ? `Showing an ${labels.length} x ${rows.length} preview of the full ${matrix.total_columns || dimension}-sample matrix.`
    : `${dimension} samples. Axis IDs map to full sample names in the key.`;
  text(headerNote, `${previewNote} Hover over the matrix to inspect a pair.`);
  header.appendChild(headerNote);
  figure.appendChild(header);

  const svg = document.createElementNS(namespace, "svg");
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.setAttribute("preserveAspectRatio", "xMidYMid meet");
  svg.setAttribute("role", "img");
  svg.setAttribute("aria-label", `${titleText}; hover cells for sample pair values`);

  labels.forEach((label, index) => {
    if (index % labelStep !== 0 && index !== labels.length - 1) return;
    const x = left + index * cell + cell / 2;
    const xLabel = document.createElementNS(namespace, "text");
    xLabel.setAttribute("x", String(x));
    xLabel.setAttribute("y", String(top - 10));
    xLabel.setAttribute("fill", "#494b44");
    xLabel.setAttribute("font-family", "Inter, Arial, sans-serif");
    xLabel.setAttribute("font-size", String(labelFont));
    xLabel.setAttribute("font-weight", "700");
    xLabel.setAttribute("text-anchor", "middle");
    xLabel.textContent = dciferDisplayId(labelMap, label, index);
    svg.appendChild(xLabel);
  });

  const inspector = document.createElement("div");
  inspector.className = "dcifer-heatmap-inspector is-idle";
  inspector.setAttribute("aria-live", "polite");
  const inspectorPrompt = document.createElement("span");
  inspectorPrompt.className = "dcifer-heatmap-inspector__prompt";
  text(inspectorPrompt, "Hover a square");
  const inspectorPair = document.createElement("strong");
  inspectorPair.className = "dcifer-heatmap-inspector__pair";
  text(inspectorPair, "Inspect any sample pair");
  const inspectorValue = document.createElement("span");
  inspectorValue.className = "dcifer-heatmap-inspector__value";
  text(inspectorValue, `${matrix.value_label || titleText} will appear here.`);
  inspector.append(inspectorPrompt, inspectorPair, inspectorValue);

  const showPair = (rowDisplayId, rowLabel, columnDisplayId, columnFullLabel, value) => {
    inspector.classList.remove("is-idle");
    text(inspectorPrompt, `${rowDisplayId} compared with ${columnDisplayId}`);
    text(inspectorPair, `${rowLabel}  vs  ${columnFullLabel}`);
    text(inspectorValue, `${matrix.value_label || titleText}: ${matrixValueText(kind, value)}`);
  };

  rows.forEach((rowItem, rowIndex) => {
    const rowLabel = String(rowItem.sample_id || "");
    const rowDisplayId = dciferDisplayId(labelMap, rowLabel, rowIndex);
    const y = top + rowIndex * cell;
    if (rowIndex % labelStep === 0 || rowIndex === rows.length - 1) {
      const yLabel = document.createElementNS(namespace, "text");
      yLabel.setAttribute("x", String(left - 8));
      yLabel.setAttribute("y", String(y + cell * 0.68));
      yLabel.setAttribute("fill", "#494b44");
      yLabel.setAttribute("font-family", "Inter, Arial, sans-serif");
      yLabel.setAttribute("font-size", String(labelFont));
      yLabel.setAttribute("font-weight", "700");
      yLabel.setAttribute("text-anchor", "end");
      yLabel.textContent = rowDisplayId;
      svg.appendChild(yLabel);
    }

    labels.forEach((columnLabel, columnIndex) => {
      const columnFullLabel = String(columnLabel || "");
      const columnDisplayId = dciferDisplayId(labelMap, columnFullLabel, columnIndex);
      const value = Array.isArray(rowItem.values) ? rowItem.values[columnIndex] : null;
      const rect = document.createElementNS(namespace, "rect");
      rect.classList.add("dcifer-heatmap-cell");
      rect.setAttribute("x", String(left + columnIndex * cell));
      rect.setAttribute("y", String(y));
      rect.setAttribute("width", String(cell));
      rect.setAttribute("height", String(cell));
      rect.setAttribute("fill", dciferHeatmapColor(kind, value));

      const hoverText = `${rowDisplayId} (${rowLabel}) vs ${columnDisplayId} (${columnFullLabel}): ${matrix.value_label || titleText} = ${matrixValueText(kind, value)}`;
      const nativeTitle = document.createElementNS(namespace, "title");
      nativeTitle.textContent = hoverText;
      rect.appendChild(nativeTitle);
      rect.addEventListener("mouseenter", () => {
        showPair(rowDisplayId, rowLabel, columnDisplayId, columnFullLabel, value);
      });
      svg.appendChild(rect);
    });
  });

  svg.addEventListener("mousemove", (event) => {
    const bounds = svg.getBoundingClientRect();
    if (!bounds.width || !bounds.height) return;
    const svgX = (event.clientX - bounds.left) * (width / bounds.width);
    const svgY = (event.clientY - bounds.top) * (height / bounds.height);
    const columnIndex = Math.floor((svgX - left) / cell);
    const rowIndex = Math.floor((svgY - top) / cell);
    if (rowIndex < 0 || rowIndex >= rows.length || columnIndex < 0 || columnIndex >= labels.length) return;
    const rowItem = rows[rowIndex];
    const rowLabel = String(rowItem.sample_id || "");
    const columnLabel = String(labels[columnIndex] || "");
    const value = Array.isArray(rowItem.values) ? rowItem.values[columnIndex] : null;
    showPair(
      dciferDisplayId(labelMap, rowLabel, rowIndex),
      rowLabel,
      dciferDisplayId(labelMap, columnLabel, columnIndex),
      columnLabel,
      value
    );
  });

  const viewport = document.createElement("div");
  viewport.className = "dcifer-heatmap-viewport";
  viewport.appendChild(svg);
  figure.appendChild(viewport);
  figure.appendChild(inspector);

  const caption = document.createElement("figcaption");
  const captionTitle = document.createElement("span");
  text(captionTitle, "Axis labels use short IDs; expand the sample key for full names.");
  caption.appendChild(captionTitle);

  if (plotInfo?.download_url) {
    const button = document.createElement("button");
    button.type = "button";
    text(button, "Download original PNG");
    button.addEventListener("click", () => {
      window.location.href = plotInfo.download_url;
    });
    caption.appendChild(button);
  }

  figure.appendChild(caption);
  figure.appendChild(renderDciferHeatmapKey(labelMap.entries));
  return figure;
}

function renderDciferPlots(plots, matrices = {}) {
  const gallery = $("#dcifer-plot-gallery");
  if (!gallery) return;
  gallery.replaceChildren();
  const items = Array.isArray(plots) ? plots.filter((plot) => plot && plot.exists && plot.view_url) : [];
  const matrixItems = [
    ["relatedness", matrices?.relatedness],
    ["pvalue", matrices?.pvalue],
  ].filter(([, matrix]) => matrixHasValues(matrix));
  const visibleCount = matrixItems.length || items.length;
  updatePlotJump(
    gallery,
    $("#dcifer-plot-count"),
    $("#dcifer-view-plots"),
    visibleCount,
    "No Dcifer heatmaps available yet.",
    "Dcifer heatmap",
    "Dcifer heatmaps"
  );
  if (matrixItems.length) {
    gallery.hidden = false;
    const grid = document.createElement("div");
    grid.className = "dcifer-interactive-grid";
    matrixItems.forEach(([kind, matrix]) => {
      grid.appendChild(renderDciferHeatmapCard(kind, matrix, dciferPlotForKind(items, kind)));
    });
    gallery.appendChild(grid);
    return;
  }

  if (!items.length) {
    gallery.hidden = true;
    return;
  }
  gallery.hidden = false;
  items.forEach((plot) => {
    const figure = document.createElement("figure");
    figure.className = "dcifer-plot-card";

    const img = document.createElement("img");
    img.src = plot.view_url;
    img.alt = plot.title || plot.filename || "Dcifer plot";
    img.loading = "lazy";
    figure.appendChild(img);

    const caption = document.createElement("figcaption");
    const title = document.createElement("span");
    text(title, plot.title || plot.filename || "Dcifer plot");
    caption.appendChild(title);

    if (plot.download_url) {
      const button = document.createElement("button");
      button.type = "button";
      text(button, "Download original PNG");
      button.addEventListener("click", () => {
        window.location.href = plot.download_url;
      });
      caption.appendChild(button);
    }

    figure.appendChild(caption);
    gallery.appendChild(figure);
  });
}

function renderDciferPairs(pairs) {
  const tbody = $("#dcifer-pairs-table");
  if (!tbody) return;
  tbody.replaceChildren();
  const rows = Array.isArray(pairs) ? pairs : [];
  if (!rows.length) {
    tbody.appendChild(emptyRow(6, "No pairwise relatedness rows available."));
    return;
  }
  rows.forEach((item) => {
    tbody.appendChild(row([
      displayMissing(item.sample_a),
      displayMissing(item.sample_b),
      formatNumber(item.estimate, 3),
      formatPValue(item.p_value),
      formatPValue(item.q_value),
      displayMissing(item.comparison_type)
    ]));
  });
}

async function loadDciferResults() {
  const out = encodeURIComponent(dciferOutdir());
  const mode = encodeURIComponent(analysisApiMode());
  try {
    const payload = await fetchJson(`/api/dcifer/results?out=${out}&mode=${mode}`);
    const state = payload.state || {};
    renderDciferReadiness(payload.readiness || null);
    const status = state.status || "idle";
    const resultsPanel = $("#dcifer-results");

    if (status === "complete") {
      if (resultsPanel) resultsPanel.hidden = false;
      setPill($("#dcifer-status"), "Complete", "ok");
      text(
        $("#dcifer-results-input-mode"),
        analysisResultLabel(payload.analysis_mode)
      );
      const summary = payload.summary || {};
      text($("#dcifer-samples"), formatNumber(summary.samples, 0));
      text($("#dcifer-pairs"), formatNumber(summary.pairs, 0));
      text($("#dcifer-max-relatedness"), formatNumber(summary.max_relatedness, 3));
      text($("#dcifer-q-le-alpha"), formatNumber(summary.q_le_alpha, 0));
      text(
        $("#dcifer-summary-note"),
        summary.caveat || "False-discovery-rate adjusted q-values are exploratory unless allele frequencies come from an adequate representative background population."
      );
      renderDciferPlots(payload.plots || [], payload.matrices || {});
      renderDciferPairs(payload.pairs || []);

      const files = payload.files || {};
      enableDciferDownload("#dcifer-dl-pairs", files.pairwise_relatedness);
      enableDciferDownload("#dcifer-dl-coi", files.coi);
      enableDciferDownload("#dcifer-dl-input", files.input);
      enableDciferDownload("#dcifer-dl-matrix", files.relatedness_matrix);
    } else if (status === "running") {
      setPill($("#dcifer-status"), "Running", "warn");
      setAnalysisModeControlsDisabled(true);
      const runBtn = $("#dcifer-run");
      const msg = $("#dcifer-message");
      if (runBtn) runBtn.disabled = true;
      if (msg && !msg.textContent.trim()) {
        msg.className = "inline-message";
        text(msg, "Dcifer is running. Results will appear automatically when complete.");
      }
      startDciferPolling();
    } else if (status === "failed") {
      if (resultsPanel) resultsPanel.hidden = true;
      text($("#dcifer-summary-note"), "");
      renderDciferPairs([]);
      renderDciferPlots([]);
      setPill($("#dcifer-status"), "Failed", "bad");
      setAnalysisModeControlsDisabled(false);
    }
  } catch (_error) {
    // Results may not exist yet.
    renderDciferPlots([]);
  }
}

function enableDciferDownload(selector, fileInfo) {
  const btn = $(selector);
  if (!btn) return;
  if (fileInfo && fileInfo.exists && fileInfo.download_url) {
    btn.disabled = false;
    btn.onclick = () => { window.location.href = fileInfo.download_url; };
  } else {
    btn.disabled = true;
    btn.onclick = null;
  }
}

function bindDciferEvents() {
  const toggle = $("#dcifer-enable");
  if (toggle) toggle.addEventListener("change", handleDciferToggle);
  const runBtn = $("#dcifer-run");
  if (runBtn) runBtn.addEventListener("click", runDcifer);
  [
    "#dcifer-coi-lrank",
    "#dcifer-ibd-grid-nr",
    "#dcifer-alpha"
  ].forEach((selector) => {
    const control = $(selector);
    if (control) control.addEventListener("change", saveSettings);
  });
}

// ---------------------------------------------------------------------------
// Init
// ---------------------------------------------------------------------------

async function init() {
  const startupParams = new URLSearchParams(window.location.search);
  const startupOutdir = String(startupParams.get("out") || "").trim();
  const startupTab = String(startupParams.get("tab") || "").trim().toLowerCase();
  restoreSettings(await loadPersistedSettings());
  bindEvents();
  bindPlotWheelScrolling();
  bindDinemitesEvents();
  bindDciferEvents();
  selectTab("inputs");
  setScanReady(false);
  updateScanButtonAvailability();
  setSampleValidationVisible(false);
  text($("#run-button"), "Start run");
  try {
    await loadHealth();
  } catch (_error) {
    renderCommonPaths([]);
  }
  if ($("#kelt-barcode-map")?.value.trim()) inspectKeltBarcodeMap();
  resetRunDisplay();
  const activeRun = await fetchActiveRun();
  if (activeRun?.active && activeRun.outdir) {
    activeOutdir = activeRun.outdir;
    $("#results-outdir").value = activeRun.outdir;
    syncGeneratedSampleSheetPath();
    updateSamplePathHelp();
    saveSettings();
  }
  if (startupOutdir) {
    activeOutdir = startupOutdir;
    $("#results-outdir").value = startupOutdir;
    syncGeneratedSampleSheetPath();
    updateSamplePathHelp();
  }
  const status = await refreshAllRunState();
  if (status?.active) {
    selectTab("run");
    startPolling();
  } else if (["results", "qc", "dinemites", "dcifer"].includes(startupTab)) {
    selectTab(startupTab);
  }
  // Restore DINEMITES toggle state and check for existing results
  handleDinemitesDayCutoffToggle();
  updateDinemitesModelSettingsVisibility();
  handleDinemitesToggle();
  loadDinemitesResults();
  handleDciferToggle();
  loadDciferResults();
}

init();
