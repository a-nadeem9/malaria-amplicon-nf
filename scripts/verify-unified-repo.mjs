import fs from "node:fs";
import path from "node:path";

const root = process.cwd();
const failures = [];

function read(file) {
  return fs.readFileSync(path.join(root, file), "utf8");
}

function exists(file) {
  return fs.existsSync(path.join(root, file));
}

function expect(condition, message) {
  if (!condition) failures.push(message);
}

const packageJson = JSON.parse(read("package.json"));
const tauriConfig = JSON.parse(read("src-tauri/tauri.conf.json"));
const libRs = read("src-tauri/src/lib.rs");
const indexHtml = read("index.html");
const launcherJs = read("src/main.js");
const readme = read("README.md");
const styles = read("src/styles.css");
const releaseWorkflow = read(".github/workflows/release.yml");
const postProcDada2 = read("runtime/workflow/scripts/postProc_dada2.R");
const cargoToml = read("src-tauri/Cargo.toml");
const runtimePyproject = read("runtime/pyproject.toml");
const runtimeInit = read("runtime/src/simplseq/__init__.py");
const installer = read("runtime/install-simplseq.sh");
const desktopCapability = read("src-tauri/capabilities/default.json");
const backendJs = read("runtime/gui/static/js/app.js");
const flaskApp = read("runtime/gui/flask_app.py");
const runtimeEnvironment = read("runtime/environment.yml");
const linuxRuntimeLock = read("runtime/locks/linux-64-explicit.txt");
const macRuntimeLock = read("runtime/locks/osx-64-explicit.txt");
const desktopRuntimePrep = read("scripts/prepare-desktop-runtime.mjs");
const runDinemites = read("runtime/workflow/scripts/run_dinemites.R");

expect(packageJson.name === "malaria-amplicon-nf", "package.json name must be malaria-amplicon-nf");
expect(tauriConfig.productName === "malaria-amplicon-nf", "Tauri productName must be malaria-amplicon-nf");
const version = packageJson.version;
expect(tauriConfig.version === version, "Tauri config version must match package.json");
expect(
  tauriConfig.bundle?.macOS?.signingIdentity === "-",
  "macOS release bundles must use an explicit ad-hoc signature"
);
expect(
  tauriConfig.build?.beforeBuildCommand ===
    "node scripts/prepare-desktop-runtime.mjs && node node_modules/vite/bin/vite.js build",
  "Tauri builds must prepare the bundled workflow before invoking Vite"
);
expect(cargoToml.includes(`version = "${version}"`), "Cargo version must match package.json");
expect(runtimePyproject.includes(`version = "${version}"`), "Runtime package version must match package.json");
expect(runtimeInit.includes(`__version__ = "${version}"`), "Runtime Python version must match package.json");
expect(libRs.includes(`const RUNTIME_VERSION: &str = "v${version}"`), "Launcher runtime version must match package.json");
expect(
  libRs.includes('include_str!("../../runtime/install-simplseq.sh")'),
  "Desktop launcher must embed the versioned runtime installer"
);
expect(installer.includes(`SIMPLSEQ_VERSION:-v${version}`), "Runtime installer version must match package.json");
expect(
  readme.includes("https://github.com/a-nadeem9/malaria-amplicon-nf/releases"),
  "README must link to the desktop releases page"
);
expect(indexHtml.includes("malaria-amplicon-nf"), "index.html must use the new app name");
expect(indexHtml.includes("/assets/launcher-logo.png"), "index.html must use the optimized launcher logo");
expect(!indexHtml.includes("Local analysis app"), "launcher header must not show the local-analysis eyebrow");
expect(!indexHtml.includes("Desktop app"), "loading panel must not show the desktop-app eyebrow");
expect(!indexHtml.includes("First-time setup"), "setup panel must not show the first-time-setup eyebrow");
expect(!indexHtml.includes("progress-rail"), "launcher must not render the old progress rail");
expect(indexHtml.includes("runtime-spinner"), "launcher must render a minimal runtime spinner");
expect(readme.includes("malaria-amplicon-nf"), "README must describe the unified repo");
expect(
  readme.includes("assets/malaria-amplicon-nf-readme-banner.png"),
  "README must use the branded logo + wordmark banner"
);
expect(styles.includes("#005c68"), "styles.css must include the sampled dark teal text color");
expect(styles.includes("#20988c"), "styles.css must include the sampled section teal color");
expect(styles.includes("#ec1850"), "styles.css must include the sampled accent red color");
expect(styles.includes("Aptos"), "styles.css must use the selected Aptos-like font stack");
expect(styles.includes("--wordmark-font"), "styles.css must define the standalone wordmark font stack");
expect(styles.includes("font-family: var(--wordmark-font)"), "brand name must use the standalone wordmark font stack");
expect(!styles.includes(".progress-rail"), "styles.css must not keep the old progress rail styles");
expect(styles.includes(".runtime-spinner"), "styles.css must define the minimal runtime spinner");
expect(
  postProcDada2.includes("write.fasta(lapply(seqs, s2c)"),
  "postProc_dada2.R must write ASV FASTA sequences with lapply, not header-only sapply output"
);
expect(
  runDinemites.includes('na.strings = c("", "NA")') &&
    runDinemites.includes('allele = na_if(trimws(as.character(.data$allele)), "")') &&
    runDinemites.includes("the filled dataset contains a blank allele"),
  "DINEMITES must preserve blank allele rows as visit metadata, never as alleles"
);
expect(
  runDinemites.includes("dataset_filled <- add_qpcr_times") &&
    runDinemites.includes("imputed_datasets <- impute_dataset") &&
    runDinemites.includes("has_qpcr_imputation <- TRUE") &&
    !runDinemites.includes("n_unknown_genotypes"),
  "DINEMITES must impute PCR-positive visits that lack genotype calls"
);
expect(
  runDinemites.includes('format(dates, "%d %b %Y")') &&
    runDinemites.includes("guide = guide_axis(n.dodge = 2)"),
  "DINEMITES plots must show full collection dates on non-overlapping axis rows"
);

expect(!libRs.includes("curl -fsSL {INSTALL_URL} | bash"), "desktop launcher must not pipe a remote script into bash");
expect(
  Array.isArray(tauriConfig.bundle?.resources) &&
    tauriConfig.bundle.resources.includes("resources/runtime/"),
  "desktop bundles must include the managed workflow source"
);
expect(
  libRs.includes('resolve("resources/runtime", BaseDirectory::Resource)') &&
    libRs.includes("SIMPLSEQ_BUNDLED_RUNTIME_DIR") &&
    installer.includes('BUNDLED_RUNTIME_DIR="${SIMPLSEQ_BUNDLED_RUNTIME_DIR:-}"'),
  "desktop setup must install its bundled workflow instead of downloading matching release assets"
);
expect(
  desktopRuntimePrep.includes('path.join(repoRoot, "runtime")') &&
    desktopRuntimePrep.includes('path.join(repoRoot, "src-tauri", "resources", "runtime")'),
  "desktop runtime preparation must copy the repository runtime into Tauri resources"
);
expect(installer.includes("command -v python3"), "fresh Ubuntu must be able to unpack micromamba without bzip2");
expect(!installer.includes("api.github.com"), "runtime setup must not consume GitHub API quota");
expect(!installer.includes("install_github"), "runtime setup must not install DINEMITES through the GitHub API");
expect(!runtimeEnvironment.includes("r-devtools"), "runtime must not carry the obsolete devtools installer stack");
expect(flaskApp.includes("from waitress import serve"), "desktop backend must use the embedded Waitress server");
expect(
  linuxRuntimeLock.includes("waitress-3.0.2") && macRuntimeLock.includes("waitress-3.0.2"),
  "pinned desktop environments must include Waitress"
);
expect(installer.includes("import flask, waitress"), "runtime setup must verify the local GUI server dependencies");
expect(
  launcherJs.includes("Local interface could not start") && launcherJs.includes('dataset.action = "backend"'),
  "launcher must distinguish backend startup failures from installation failures"
);
expect(flaskApp.includes("launch_job("), "downstream analyses must use durable worker processes");
expect(!flaskApp.includes("DINEMITES_THREADS"), "DINEMITES must not be owned by a daemon thread");
expect(!flaskApp.includes("DCIFER_THREADS"), "Dcifer must not be owned by a daemon thread");
expect(!desktopCapability.includes('"remote"'), "remote backend pages must not receive Tauri capabilities");
expect(!tauriConfig.app?.withGlobalTauri, "the global Tauri bridge must stay disabled");
expect(
  tauriConfig.app?.windows?.[0]?.theme === "Light" &&
    tauriConfig.app.windows[0].decorations !== false &&
    libRs.includes("DwmSetWindowAttribute") &&
    libRs.includes("style_native_titlebar(&window)"),
  "Windows must retain native controls while styling the title bar with the app palette"
);
expect(!backendJs.includes("window.__TAURI__"), "the Flask backend must not depend on a remote Tauri bridge");
expect(libRs.includes("GET /api/health HTTP/1.1"), "desktop launcher must wait for the backend health endpoint");
expect(libRs.includes('mode: "wsl_missing"'), "desktop launcher must detect missing WSL separately");
expect(
  libRs.includes('const DEFAULT_WSL_DISTRO_NAME: &str = "Ubuntu"') &&
    libRs.includes("const LEGACY_MANAGED_WSL_DISTRO_NAMES: [&str; 3]") &&
    libRs.includes('"malaria-amplicon-nf-runtime-v1"') &&
    libRs.includes("find_runtime_wsl_distro"),
  "Windows must reuse existing WSL distributions and retain legacy runtime compatibility"
);
expect(
  !libRs.includes("WSL_ROOTFS_URL") && !libRs.includes('"--import"'),
  "Windows must not download and import a second private Ubuntu image"
);
expect(
  libRs.includes('command.args(["--distribution", distro_name, "--user", "root", "--exec"])') &&
    libRs.includes("runtime_wsl_command()"),
  "WSL commands must avoid interactive Ubuntu account setup"
);
expect(
  libRs.includes('command.args(["--list", "--quiet"])') &&
    libRs.includes('name.eq_ignore_ascii_case(DEFAULT_WSL_DISTRO_NAME)') &&
    libRs.includes('starts_with("ubuntu")') &&
    !libRs.includes(".or_else(|| registered.first().cloned())"),
  "Windows setup must prefer Ubuntu and must not install into an unrelated WSL distribution"
);
expect(
  libRs.includes('"--install"') &&
    libRs.includes('"--distribution"') &&
    libRs.includes("DEFAULT_WSL_DISTRO_NAME") &&
    libRs.includes('"--no-launch"') &&
    libRs.includes("Start-Process -FilePath 'wsl.exe'") &&
    libRs.includes("-Verb RunAs") &&
    !libRs.includes('return Err("Windows Subsystem for Linux is not installed.'),
  "Windows setup must install standard Ubuntu and elevate when the WSL platform is absent"
);
expect(
  flaskApp.includes("is_windows_network_path") &&
    flaskApp.includes("mount the share inside WSL"),
  "WSL folder pickers must reject inaccessible Windows network-share paths clearly"
);
expect(
  libRs.includes('"runtime-progress"') && libRs.includes('phase: "packages-install"'),
  "first-run setup must expose real progress stages"
);
expect(
  libRs.includes('"runtime-log"') &&
    libRs.includes("spawn_runtime_log_reader") &&
    launcherJs.includes('listen("runtime-log"'),
  "first-run setup must stream live installer output into the launcher"
);
expect(
  libRs.includes("installer_progress_from_log_line") &&
    launcherJs.includes("const PROGRESS_STAGES") &&
    launcherJs.includes("Estimated time remaining:") &&
    indexHtml.includes('id="progress-percent"') &&
    indexHtml.includes('id="progress-eta"'),
  "first-run setup must expose checkpoint-driven percentage progress and an estimated time remaining"
);

expect(exists("runtime/main.nf"), "runtime/main.nf must exist");
expect(exists("runtime/nextflow.config"), "runtime/nextflow.config must exist");
expect(exists("runtime/src/simplseq/cli.py"), "runtime Python CLI must exist");
expect(exists("runtime/workflow/scripts/AmpliconPipeline.py"), "runtime workflow scripts must exist");
expect(exists("runtime/install-simplseq.sh"), "runtime installer must exist");
expect(
  installer.includes('DINEMITES_LDFLAGS="${LDFLAGS:-}"') &&
    installer.includes("-mlinker-version=0") &&
    installer.includes('LDFLAGS="$DINEMITES_LDFLAGS"'),
  "Apple Silicon setup must avoid the macOS 26 versioned libLTO linker failure"
);

expect(exists("assets/logo.png"), "new logo must exist at assets/logo.png");
expect(exists("assets/launcher-logo.png"), "optimized launcher logo must exist");
expect(
  exists("assets/malaria-amplicon-nf-readme-banner.png"),
  "README banner must exist at assets/malaria-amplicon-nf-readme-banner.png"
);

expect(releaseWorkflow.includes("npm run verify:repo"), "release workflow must verify the unified repo before packaging");
expect(
  releaseWorkflow.includes("Verify macOS ad-hoc bundle signature") &&
    releaseWorkflow.includes("codesign --verify --deep --strict"),
  "release workflow must reject malformed macOS ad-hoc signatures"
);
expect(
    releaseWorkflow.includes("malaria-amplicon-nf-${version}-Setup.exe") &&
    releaseWorkflow.includes("malaria-amplicon-nf-${version}-arm64.dmg"),
  "release workflow must publish two versioned desktop installer assets"
);
expect(
  !releaseWorkflow.includes("macOS-Intel.dmg") &&
    !releaseWorkflow.includes("x86_64-apple-darwin") &&
    !releaseWorkflow.includes("macos-15-intel") &&
    !releaseWorkflow.includes("Linux.AppImage") &&
    !releaseWorkflow.includes("runtime.tar.gz") &&
    !releaseWorkflow.includes("SHA256SUMS.txt"),
  "public releases must contain only Windows and Apple Silicon installers"
);
expect(
  releaseWorkflow.includes("Verify public release downloads"),
  "release workflow must verify assets through their public download URLs"
);

if (failures.length) {
  console.error("Unified repo verification failed:");
  for (const failure of failures) console.error(`- ${failure}`);
  process.exit(1);
}

console.log("Unified repo verification passed.");
