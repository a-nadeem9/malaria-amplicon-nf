import { invoke, isTauri } from "@tauri-apps/api/core";
import { listen } from "@tauri-apps/api/event";
import { getCurrentWindow } from "@tauri-apps/api/window";
import "./styles.css";

const hasTauriRuntime = isTauri();
const appWindow = hasTauriRuntime ? getCurrentWindow() : null;
const previewMode = new URLSearchParams(window.location.search).get("preview");

const loadingView = document.querySelector("#loading-view");
const setupView = document.querySelector("#setup-view");
const runtimeBadge = document.querySelector("#runtime-badge");
const modeCopy = document.querySelector("#mode-copy");
const runtimeLine = document.querySelector(".runtime-line");
const loadingTitleEl = document.querySelector("#loading-title");
const loadingDetailEl = document.querySelector("#loading-detail");
const progressBarEl = document.querySelector("#launcher-progress");
const progressLabelEl = document.querySelector("#progress-label");
const progressPercentEl = document.querySelector("#progress-percent");
const progressEtaEl = document.querySelector("#progress-eta");
const setupStatusEl = document.querySelector("#setup-status");
const setupTitleEl = document.querySelector("#setup-title");
const setupCopyEl = document.querySelector("#setup-copy");
const setupBadgeEl = document.querySelector("#setup-badge");
const installButton = document.querySelector("#install");
const retryButton = document.querySelector("#retry");
const shell = document.querySelector(".shell");
const installConsoles = [...document.querySelectorAll(".install-console")];
const installConsoleOutputs = [...document.querySelectorAll(".install-console-output")];
const runtimeLogLines = [];

const PROGRESS_STAGES = Object.freeze({
  preview: { start: 12, end: 76, duration: 12, after: 0, label: "Loading preview", showEta: false },
  "runtime-check": { start: 6, end: 34, duration: 8, after: 0, label: "Checking runtime", showEta: false },
  "environment-detect": { start: 4, end: 8, duration: 10, after: 960, label: "Checking Ubuntu", showEta: true },
  "ubuntu-install": { start: 6, end: 14, duration: 180, after: 945, label: "Installing Ubuntu", showEta: true },
  "ubuntu-ready": { start: 14, end: 16, duration: 4, after: 940, label: "Ubuntu ready", showEta: true },
  "packages-install": { start: 4, end: 7, duration: 3, after: 940, label: "Preparing installation", showEta: true },
  "install-files": { start: 7, end: 11, duration: 4, after: 925, label: "Installing workflow files", showEta: true },
  "download-runtime": { start: 7, end: 10, duration: 60, after: 925, label: "Downloading workflow files", showEta: true },
  "verify-download": { start: 10, end: 11, duration: 10, after: 925, label: "Verifying workflow files", showEta: true },
  "install-micromamba": { start: 11, end: 18, duration: 45, after: 880, label: "Preparing package manager", showEta: true },
  "create-runtime": { start: 18, end: 74, duration: 600, after: 280, label: "Installing analysis environment", showEta: true },
  "install-r-packages": { start: 74, end: 91, duration: 240, after: 45, label: "Installing analysis modules", showEta: true },
  "create-launcher": { start: 91, end: 94, duration: 5, after: 40, label: "Creating local launcher", showEta: true },
  "check-path": { start: 94, end: 96, duration: 3, after: 37, label: "Finalizing launcher", showEta: true },
  "verify-runtime": { start: 96, end: 99, duration: 40, after: 0, label: "Verifying installation", showEta: true },
  "packages-ready": { start: 100, end: 100, duration: 0, after: 0, label: "Installation complete", showEta: false },
  "backend-start": { start: 72, end: 96, duration: 30, after: 0, label: "Starting local interface", showEta: false },
});

let progressValue = 0;
let activeProgressStage = PROGRESS_STAGES["runtime-check"];
let progressStageStartedAt = Date.now();
let progressTicker = null;

shell.addEventListener("mousedown", (event) => {
  if (!appWindow || event.button !== 0 || event.target.closest("button")) return;
  appWindow.startDragging().catch(() => {});
});

function showView(view) {
  loadingView.hidden = view !== "loading";
  setupView.hidden = view !== "setup";
}

function setBadge(label, mode) {
  if (!runtimeBadge) return;
  runtimeBadge.textContent = label;
  runtimeBadge.dataset.mode = mode;
}

function setRuntimeLineVisible(visible) {
  if (runtimeLine) runtimeLine.hidden = !visible;
}

function setLoadingStatus(message) {
  if (!loadingDetailEl) return;
  loadingDetailEl.textContent = message;
  loadingDetailEl.hidden = !String(message || "").trim();
}

function formatRemainingTime(seconds) {
  if (seconds <= 45) return "Estimated time remaining: less than a minute";
  const minutes = Math.max(1, Math.ceil(seconds / 60));
  return `Estimated time remaining: about ${minutes} min`;
}

function renderProgress() {
  const rounded = Math.max(0, Math.min(100, Math.round(progressValue)));
  progressBarEl?.style.setProperty("--progress", `${rounded}%`);
  progressBarEl?.setAttribute("aria-valuenow", String(rounded));
  if (progressPercentEl) progressPercentEl.textContent = `${rounded}%`;
}

function stopProgressTicker() {
  if (progressTicker !== null) {
    window.clearInterval(progressTicker);
    progressTicker = null;
  }
}

function updateProgressEstimate() {
  const stage = activeProgressStage;
  const elapsed = Math.max(0, (Date.now() - progressStageStartedAt) / 1000);
  if (stage.duration > 0) {
    const ratio = Math.min(0.92, elapsed / stage.duration);
    const estimatedProgress = stage.start + ((stage.end - stage.start) * ratio);
    progressValue = Math.max(progressValue, estimatedProgress);
  }
  renderProgress();

  if (!progressEtaEl) return;
  progressEtaEl.hidden = !stage.showEta;
  if (!stage.showEta) return;
  if (stage.duration > 0 && elapsed > stage.duration * 1.25) {
    progressEtaEl.textContent = "Taking longer than estimated; installation is still active.";
    return;
  }
  const remaining = Math.max(0, stage.duration - elapsed) + stage.after;
  progressEtaEl.textContent = formatRemainingTime(remaining);
}

function setProgressStage(phase) {
  const stage = PROGRESS_STAGES[phase] || PROGRESS_STAGES["packages-install"];
  activeProgressStage = stage;
  progressStageStartedAt = Date.now();
  progressValue = Math.max(progressValue, stage.start);
  if (progressLabelEl) progressLabelEl.textContent = stage.label;
  stopProgressTicker();
  updateProgressEstimate();
  if (stage.end < 100) {
    progressTicker = window.setInterval(updateProgressEstimate, 1000);
  }
}

function resetProgress(phase) {
  progressValue = 0;
  setProgressStage(phase);
}

function setLoadingPhase(title, detail, phase) {
  if (loadingTitleEl) loadingTitleEl.textContent = title;
  setLoadingStatus(detail);
  if (phase) setProgressStage(phase);
}

function renderRuntimeLog() {
  const hasOutput = runtimeLogLines.length > 0;
  const content = runtimeLogLines.join("\n");
  installConsoles.forEach((consoleEl) => {
    consoleEl.hidden = !hasOutput;
  });
  installConsoleOutputs.forEach((outputEl) => {
    outputEl.textContent = content;
    outputEl.scrollTop = outputEl.scrollHeight;
  });
}

function resetRuntimeLog() {
  runtimeLogLines.length = 0;
  renderRuntimeLog();
}

function appendRuntimeLog(payload) {
  const line = String(payload?.line || "").trim();
  if (!line) return;
  runtimeLogLines.push(line);
  if (runtimeLogLines.length > 120) {
    runtimeLogLines.splice(0, runtimeLogLines.length - 120);
  }
  renderRuntimeLog();
}

function setSetupStatus(message) {
  setupStatusEl.textContent = message;
}

function isWslPlatformMissing(runtime) {
  return runtime.mode === "wsl_missing";
}

function isLinuxEnvironmentMissing(runtime) {
  return runtime.mode === "wsl_distribution_missing";
}

function humanizeSetupError(error) {
  const message = String(error || "");
  if (/rate limit|403/i.test(message)) {
    return "The DINEMITES download was temporarily limited. Retry setup; no GitHub account or token is required.";
  }
  if (/bzip2|Cannot exec/i.test(message)) {
    return "The runtime archive could not be unpacked. Retry setup to use the built-in compatibility fallback.";
  }
  if (/Could not resolve|Could not connect|connection|timed out|curl:\s*\([56728]/i.test(message)) {
    return "The runtime download was interrupted. Check your internet connection, then retry setup.";
  }
  if (/restart Windows|reboot/i.test(message)) {
    return "Windows needs a restart to finish enabling WSL. Restart the computer, then open malaria-amplicon-nf again.";
  }
  if (/Linux environment|distribution|Ubuntu/i.test(message)) {
    return "Ubuntu could not be prepared in WSL. Restart Windows if requested, then retry setup.";
  }
  return "Setup stopped before the analysis runtime finished installing. Retry setup to continue.";
}

function configureSetup(runtime) {
  stopProgressTicker();
  const wslMissing = isWslPlatformMissing(runtime);
  const environmentMissing = isLinuxEnvironmentMissing(runtime);
  setupTitleEl.textContent = wslMissing
    ? "Windows setup required"
    : environmentMissing
      ? "Ubuntu setup is required"
      : "Runtime setup was interrupted";
  setupCopyEl.textContent = wslMissing
    ? "This app needs Windows Subsystem for Linux. Install WSL once, restart Windows if requested, then check again."
    : environmentMissing
      ? "WSL is installed but has no Linux distribution. malaria-amplicon-nf can add standard Ubuntu, then install its workflow packages there."
      : "No additional software or account is required. Retry setup to continue the workflow package installation.";
  setupBadgeEl.textContent = wslMissing ? "WSL required" : "Needs retry";
  installButton.textContent = wslMissing
    ? "Install WSL"
    : environmentMissing
      ? "Install Ubuntu"
      : "Retry Setup";
  installButton.dataset.action = wslMissing
    ? "wsl-platform"
    : environmentMissing
      ? "wsl-environment"
      : "runtime";
  setSetupStatus(runtime.detail || runtime.message || "Setup is required.");
}

function configureBackendStartupError(error) {
  stopProgressTicker();
  showView("setup");
  setupTitleEl.textContent = "Local interface could not start";
  setupCopyEl.textContent = "The analysis runtime is installed, but its local interface did not start. Try again or review the status below.";
  setupBadgeEl.textContent = "Startup issue";
  installButton.textContent = "Try Again";
  installButton.dataset.action = "backend";
  installButton.disabled = false;
  retryButton.disabled = false;
  setSetupStatus(`The local interface could not start. ${String(error || "Unknown startup error.")}`);
}

async function openSimplseq() {
  showView("loading");
  if (runtimeBadge) runtimeBadge.hidden = true;
  setRuntimeLineVisible(false);
  setBadge("Opening", "running");
  modeCopy.textContent = "malaria-amplicon-nf runtime is ready.";
  setLoadingPhase("Opening malaria-amplicon-nf...", "Starting the local workflow interface...", "backend-start");

  const result = await invoke("start_backend");
  setLoadingStatus("");
}

async function openSimplseqWithRecovery() {
  try {
    await openSimplseq();
    return true;
  } catch (error) {
    configureBackendStartupError(error);
    return false;
  }
}

async function installRuntimeAndOpen(runtime, resetLog = true) {
  showView("loading");
  if (runtimeBadge) runtimeBadge.hidden = true;
  setRuntimeLineVisible(false);
  modeCopy.textContent = runtime.mode === "outdated"
    ? "Updating malaria-amplicon-nf runtime..."
    : "Installing malaria-amplicon-nf runtime...";
  setLoadingPhase(
    runtime.mode === "outdated" ? "Updating analysis runtime..." : "Installing workflow packages...",
    "Preparing the pinned Nextflow, Python, and R packages. This is a one-time setup step.",
    "packages-install"
  );
  if (resetLog) resetRuntimeLog();

  try {
    await invoke("install_runtime");
    const updatedRuntime = await invoke("detect_runtime");
    if (!["wsl", "native"].includes(updatedRuntime.mode)) {
      throw new Error(updatedRuntime.detail || "Runtime setup did not complete.");
    }
  } catch (error) {
    showView("setup");
    configureSetup(runtime);
    installButton.disabled = false;
    retryButton.disabled = false;
    setSetupStatus(humanizeSetupError(error));
    return;
  }

  await openSimplseqWithRecovery();
}

async function installWslEnvironmentAndOpen(runtime) {
  showView("loading");
  if (runtimeBadge) runtimeBadge.hidden = true;
  setRuntimeLineVisible(false);
  modeCopy.textContent = "Preparing Ubuntu in WSL...";
  setLoadingPhase(
    "Checking Ubuntu...",
    "Using the existing WSL Linux distribution when one is available...",
    "environment-detect"
  );
  resetProgress("environment-detect");
  resetRuntimeLog();

  try {
    await invoke("install_wsl");
    const updatedRuntime = await invoke("detect_runtime");
    if (isLinuxEnvironmentMissing(updatedRuntime) || isWslPlatformMissing(updatedRuntime)) {
      throw new Error(updatedRuntime.detail || "The WSL Linux environment is not ready.");
    }
    await installRuntimeAndOpen(updatedRuntime, false);
  } catch (error) {
    showView("setup");
    configureSetup(runtime);
    installButton.disabled = false;
    retryButton.disabled = false;
    setSetupStatus(humanizeSetupError(error));
  }
}

async function checkAndOpen() {
  if (!hasTauriRuntime) {
    showView(previewMode === "loading" ? "loading" : "setup");
    if (runtimeBadge) runtimeBadge.hidden = true;
    setRuntimeLineVisible(false);
    resetProgress("preview");
    setLoadingPhase("Loading...", "Previewing the desktop launcher.", "preview");
    setSetupStatus("Preview only. In the packaged app, this button installs the runtime for WSL, Linux, or macOS.");
    return;
  }

  showView("loading");
  resetProgress("runtime-check");
  if (runtimeBadge) runtimeBadge.hidden = true;
  setRuntimeLineVisible(false);
  setBadge("Checking", "starting");
  modeCopy.textContent = "Checking local malaria-amplicon-nf runtime...";
  setLoadingPhase("Checking runtime...", "Looking for the local analysis environment...", "runtime-check");

  const runtime = await invoke("detect_runtime");
  if (isWslPlatformMissing(runtime)) {
    showView("setup");
    configureSetup(runtime);
    return;
  }
  if (isLinuxEnvironmentMissing(runtime)) {
    await installWslEnvironmentAndOpen(runtime);
    return;
  }
  if (runtime.mode === "missing" || runtime.mode === "outdated") {
    await installRuntimeAndOpen(runtime);
    return;
  }

  await openSimplseqWithRecovery();
}

installButton.addEventListener("click", async () => {
  if (!hasTauriRuntime) {
    setSetupStatus("Preview only. The packaged desktop app can run the runtime installer for the user.");
    return;
  }

  installButton.disabled = true;
  retryButton.disabled = true;
  const action = installButton.dataset.action;
  if (action === "backend") {
    await openSimplseqWithRecovery();
    return;
  }
  const installingWsl = action === "wsl-platform" || action === "wsl-environment";
  resetProgress(installingWsl ? "environment-detect" : "packages-install");
  setSetupStatus(action === "wsl-platform"
    ? "Windows will ask for administrator approval. A restart may be required after WSL is installed."
    : action === "wsl-environment"
      ? "Installing standard Ubuntu in WSL. This is needed only when WSL has no Linux distribution."
      : "Installing malaria-amplicon-nf runtime. This can take several minutes on first run...");
  resetRuntimeLog();

  try {
    await invoke(installingWsl ? "install_wsl" : "install_runtime");
    await checkAndOpen();
  } catch (error) {
    installButton.disabled = false;
    retryButton.disabled = false;
    setSetupStatus(humanizeSetupError(error));
  }
});

retryButton.addEventListener("click", () => {
  checkAndOpen().catch((error) => {
    showView("setup");
    setSetupStatus(humanizeSetupError(error));
  });
});

if (hasTauriRuntime) {
  listen("runtime-progress", ({ payload }) => {
    if (!payload) return;
    showView("loading");
    if (runtimeBadge) runtimeBadge.hidden = true;
    setRuntimeLineVisible(false);
    setLoadingPhase(
      payload.title || "Preparing analysis runtime...",
      payload.detail || "Working...",
      payload.phase || "packages-install"
    );
  }).catch(() => {});
  listen("runtime-log", ({ payload }) => {
    appendRuntimeLog(payload);
  }).catch(() => {});
}

checkAndOpen().catch((error) => {
  showView("setup");
  configureSetup({mode: "error", detail: String(error)});
  setSetupStatus(humanizeSetupError(error));
});
