use serde::{Deserialize, Serialize};
use std::{
    env,
    fs::{self, File, OpenOptions},
    io::{BufRead, BufReader, Read, Write},
    net::{TcpListener, TcpStream},
    path::{Path, PathBuf},
    process::{Child, Command, Stdio},
    sync::Mutex,
    thread,
    time::{Duration, Instant},
};
use tauri::{
    path::BaseDirectory, AppHandle, Emitter, LogicalSize, Manager, State, Url, WebviewWindow,
    WindowEvent,
};
use tauri_plugin_dialog::DialogExt;

const RUNTIME_VERSION: &str = "v1.0.3";
const DEFAULT_WSL_DISTRO_NAME: &str = "Ubuntu";
const LEGACY_MANAGED_WSL_DISTRO_NAMES: [&str; 3] = [
    "malaria-amplicon-nf-runtime-v3",
    "malaria-amplicon-nf-runtime-v2",
    "malaria-amplicon-nf-runtime-v1",
];
const INSTALL_SCRIPT: &str = include_str!("../../runtime/install-simplseq.sh");
const APP_WINDOW_WIDTH: f64 = 1180.0;
const APP_WINDOW_HEIGHT: f64 = 820.0;
const APP_WINDOW_MIN_WIDTH: f64 = 920.0;
const APP_WINDOW_MIN_HEIGHT: f64 = 680.0;

#[cfg(windows)]
const WINDOWS_CAPTION_BACKGROUND: u32 = 0x00FF_FFFF;
#[cfg(windows)]
const WINDOWS_CAPTION_TEXT: u32 = 0x0068_5C00;
#[cfg(windows)]
const WINDOWS_CAPTION_BORDER: u32 = 0x00F0_EFE4;

#[cfg(windows)]
use std::os::windows::process::CommandExt;

#[cfg(windows)]
const CREATE_NO_WINDOW: u32 = 0x08000000;

#[derive(Default)]
struct BackendState {
    child: Mutex<Option<Child>>,
    url: Mutex<Option<String>>,
    runtime_verified: Mutex<bool>,
    picker_bridge_dir: Mutex<Option<PathBuf>>,
}

#[derive(Serialize)]
struct RuntimeInfo {
    mode: String,
    message: String,
    detail: String,
}

#[derive(Clone, Serialize)]
struct RuntimeProgress {
    phase: String,
    title: String,
    detail: String,
}

#[derive(Clone, Serialize)]
struct RuntimeLog {
    line: String,
    stream: String,
}

#[derive(Serialize)]
struct BackendLaunch {
    url: String,
}

#[derive(Deserialize)]
struct PickerRequest {
    picker_type: String,
    initial: Option<String>,
    prompt: Option<String>,
    allow_new_folder: Option<bool>,
    kind: Option<String>,
}

#[derive(Serialize)]
struct PickerResponse {
    ok: bool,
    selected: bool,
    path: String,
    error: String,
}

fn selected_dialog_path(
    path: Option<tauri_plugin_dialog::FilePath>,
) -> Result<Option<String>, String> {
    path.map(|value| {
        value
            .into_path()
            .map(|path| path.to_string_lossy().into_owned())
            .map_err(|_| "The selected item is not a local filesystem path.".to_string())
    })
    .transpose()
}

fn configure_dialog_directory<R: tauri::Runtime>(
    mut dialog: tauri_plugin_dialog::FileDialogBuilder<R>,
    initial: Option<String>,
) -> tauri_plugin_dialog::FileDialogBuilder<R> {
    let Some(initial) = initial.filter(|value| !value.trim().is_empty()) else {
        return dialog;
    };
    let path = PathBuf::from(initial);
    let directory = if path.is_dir() {
        Some(path)
    } else {
        path.parent()
            .filter(|parent| parent.is_dir())
            .map(Path::to_path_buf)
    };
    if let Some(directory) = directory {
        dialog = dialog.set_directory(directory);
    }
    dialog
}

fn picker_response(app: &AppHandle, request: PickerRequest) -> PickerResponse {
    let mut dialog = configure_dialog_directory(app.dialog().file(), request.initial);
    if let Some(prompt) = request.prompt.filter(|value| !value.trim().is_empty()) {
        dialog = dialog.set_title(prompt);
    }
    let selection = if request.picker_type == "folder" {
        dialog
            .set_can_create_directories(request.allow_new_folder.unwrap_or(false))
            .blocking_pick_folder()
    } else {
        let dialog = if request.kind.as_deref() == Some("kelt") {
            dialog.add_filter("KELT barcode maps", &["csv", "tsv"])
        } else {
            dialog.add_filter("Metadata files", &["csv", "tsv", "xlsx", "xlsm"])
        };
        dialog.blocking_pick_file()
    };
    match selected_dialog_path(selection) {
        Ok(path) => PickerResponse {
            ok: true,
            selected: path.is_some(),
            path: path.unwrap_or_default(),
            error: String::new(),
        },
        Err(error) => PickerResponse {
            ok: false,
            selected: false,
            path: String::new(),
            error,
        },
    }
}

fn start_picker_bridge(app: AppHandle) -> Result<PathBuf, String> {
    let base = env::var_os("LOCALAPPDATA")
        .or_else(|| env::var_os("XDG_RUNTIME_DIR"))
        .map(PathBuf::from)
        .unwrap_or_else(env::temp_dir);
    let directory = base
        .join("malaria-amplicon-nf")
        .join(format!("picker-{}", std::process::id()));
    fs::create_dir_all(&directory)
        .map_err(|err| format!("Could not create desktop picker bridge: {err}"))?;
    let worker_directory = directory.clone();
    thread::spawn(move || loop {
        if let Ok(entries) = fs::read_dir(&worker_directory) {
            for entry in entries.flatten() {
                let path = entry.path();
                let name = path
                    .file_name()
                    .and_then(|value| value.to_str())
                    .unwrap_or("");
                if !name.ends_with(".request.json") {
                    continue;
                }
                let response_path =
                    path.with_file_name(name.replace(".request.json", ".response.json"));
                let temporary_path = response_path.with_extension("json.tmp");
                let response = fs::read_to_string(&path)
                    .map_err(|err| err.to_string())
                    .and_then(|value| {
                        serde_json::from_str::<PickerRequest>(&value).map_err(|err| err.to_string())
                    })
                    .map(|request| picker_response(&app, request))
                    .unwrap_or_else(|error| PickerResponse {
                        ok: false,
                        selected: false,
                        path: String::new(),
                        error,
                    });
                if let Ok(encoded) = serde_json::to_vec(&response) {
                    if fs::write(&temporary_path, encoded).is_ok() {
                        let _ = fs::rename(&temporary_path, &response_path);
                    }
                }
                let _ = fs::remove_file(&path);
            }
        }
        thread::sleep(Duration::from_millis(100));
    });
    Ok(directory)
}

fn command_available(command: &str, args: &[&str]) -> bool {
    let mut command = Command::new(command);
    command.args(args);
    hide_command_window(&mut command);
    let Ok(mut child) = command
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
    else {
        return false;
    };

    let started = Instant::now();
    while started.elapsed() < Duration::from_secs(20) {
        match child.try_wait() {
            Ok(Some(status)) => return status.success(),
            Ok(None) => thread::sleep(Duration::from_millis(100)),
            Err(_) => return false,
        }
    }
    let _ = child.kill();
    let _ = child.wait();
    false
}

fn hide_command_window(command: &mut Command) {
    #[cfg(windows)]
    {
        command.creation_flags(CREATE_NO_WINDOW);
    }
    #[cfg(not(windows))]
    {
        let _ = command;
    }
}

fn decode_command_output(bytes: &[u8]) -> String {
    let looks_utf16 = bytes
        .chunks(2)
        .take(16)
        .any(|chunk| chunk.len() == 2 && chunk[1] == 0);
    if looks_utf16 {
        let units = bytes
            .chunks_exact(2)
            .map(|chunk| u16::from_le_bytes([chunk[0], chunk[1]]))
            .collect::<Vec<_>>();
        String::from_utf16_lossy(&units)
    } else {
        String::from_utf8_lossy(bytes).into_owned()
    }
}

fn capture_command(
    command: &mut Command,
    timeout: Duration,
) -> Result<(bool, String, String), String> {
    hide_command_window(command);
    let mut child = command
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|err| err.to_string())?;
    let started = Instant::now();
    loop {
        match child.try_wait() {
            Ok(Some(status)) => {
                let mut stdout = Vec::new();
                let mut stderr = Vec::new();
                if let Some(mut stream) = child.stdout.take() {
                    let _ = stream.read_to_end(&mut stdout);
                }
                if let Some(mut stream) = child.stderr.take() {
                    let _ = stream.read_to_end(&mut stderr);
                }
                return Ok((
                    status.success(),
                    decode_command_output(&stdout),
                    decode_command_output(&stderr),
                ));
            }
            Ok(None) if started.elapsed() < timeout => {
                thread::sleep(Duration::from_millis(50));
            }
            Ok(None) => {
                let _ = child.kill();
                let _ = child.wait();
                return Err(format!(
                    "command timed out after {} seconds",
                    timeout.as_secs()
                ));
            }
            Err(error) => return Err(error.to_string()),
        }
    }
}

fn registered_wsl_distributions() -> Result<Vec<String>, String> {
    let mut command = Command::new("wsl.exe");
    command.args(["--list", "--quiet"]);
    let (success, stdout, stderr) = capture_command(&mut command, Duration::from_secs(5))?;
    if !success {
        let detail = format!("{stdout}\n{stderr}");
        if detail
            .to_ascii_lowercase()
            .contains("has no installed distributions")
        {
            return Ok(Vec::new());
        }
        let detail = detail.trim();
        return Err(if detail.is_empty() {
            "Windows Subsystem for Linux is not available.".into()
        } else {
            detail.into()
        });
    }
    Ok(stdout
        .lines()
        .map(|line| line.trim_matches(|value: char| value == '\0' || value.is_whitespace()))
        .filter(|line| !line.is_empty())
        .map(String::from)
        .collect())
}

fn find_runtime_wsl_distro(registered: &[String]) -> Option<String> {
    registered
        .iter()
        .find(|name| name.eq_ignore_ascii_case(DEFAULT_WSL_DISTRO_NAME))
        .cloned()
        .or_else(|| {
            registered
                .iter()
                .find(|name| name.to_ascii_lowercase().starts_with("ubuntu"))
                .cloned()
        })
        .or_else(|| {
            LEGACY_MANAGED_WSL_DISTRO_NAMES
                .iter()
                .find_map(|candidate| {
                    registered
                        .iter()
                        .find(|name| name.eq_ignore_ascii_case(candidate))
                        .cloned()
                })
        })
}

fn wsl_command_for(distro_name: &str) -> Command {
    let mut command = Command::new("wsl.exe");
    command.args(["--distribution", distro_name, "--user", "root", "--exec"]);
    command
}

fn wsl_distro_ready(distro_name: &str) -> bool {
    let mut command = wsl_command_for(distro_name);
    command.args(["bash", "-lc", "true"]);
    hide_command_window(&mut command);
    let Ok(mut child) = command
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
    else {
        return false;
    };
    let started = Instant::now();
    while started.elapsed() < Duration::from_secs(20) {
        match child.try_wait() {
            Ok(Some(status)) => return status.success(),
            Ok(None) => thread::sleep(Duration::from_millis(100)),
            Err(_) => return false,
        }
    }
    let _ = child.kill();
    let _ = child.wait();
    false
}

fn runtime_wsl_distro_name() -> String {
    registered_wsl_distributions()
        .ok()
        .and_then(|registered| find_runtime_wsl_distro(&registered))
        .unwrap_or_else(|| DEFAULT_WSL_DISTRO_NAME.into())
}

fn runtime_wsl_command() -> Command {
    wsl_command_for(&runtime_wsl_distro_name())
}

fn windows_path_for_wsl(path: &Path) -> Result<String, String> {
    let mut raw = path.to_string_lossy().replace('\\', "/");
    if let Some(unc_path) = raw.strip_prefix("//?/UNC/") {
        raw = format!("//{unc_path}");
    } else if let Some(drive_path) = raw
        .strip_prefix("//?/")
        .or_else(|| raw.strip_prefix("//./"))
    {
        raw = drive_path.to_string();
    }
    let bytes = raw.as_bytes();
    if bytes.len() >= 2 && bytes[1] == b':' && bytes[0].is_ascii_alphabetic() {
        let drive = (bytes[0] as char).to_ascii_lowercase();
        let remainder = raw[2..].trim_start_matches('/');
        return Ok(format!("/mnt/{drive}/{remainder}"));
    }
    if raw.starts_with('/') {
        return Ok(raw);
    }
    Err(format!(
        "Could not translate the desktop picker path for WSL: {raw}"
    ))
}

#[cfg(test)]
mod path_tests {
    use super::windows_path_for_wsl;
    use std::path::Path;

    #[test]
    fn translates_regular_windows_drive_path_for_wsl() {
        let translated = windows_path_for_wsl(Path::new(
            r"C:\Users\Adina rajan\AppData\Local\malaria-amplicon-nf\resources\runtime",
        ))
        .unwrap();

        assert_eq!(
            translated,
            "/mnt/c/Users/Adina rajan/AppData/Local/malaria-amplicon-nf/resources/runtime"
        );
    }

    #[test]
    fn translates_extended_windows_drive_path_for_wsl() {
        let translated = windows_path_for_wsl(Path::new(
            r"\\?\C:\Users\Adina rajan\AppData\Local\malaria-amplicon-nf\resources\runtime",
        ))
        .unwrap();

        assert_eq!(
            translated,
            "/mnt/c/Users/Adina rajan/AppData/Local/malaria-amplicon-nf/resources/runtime"
        );
        assert!(!translated.contains("//?/"));
    }
}

fn bundled_runtime_directory(app: &AppHandle) -> Result<PathBuf, String> {
    let runtime_dir = app
        .path()
        .resolve("resources/runtime", BaseDirectory::Resource)
        .map_err(|err| format!("Could not locate the bundled workflow files: {err}"))?;

    for required in ["main.nf", "environment.yml"] {
        if !runtime_dir.join(required).is_file() {
            return Err(format!(
                "The desktop installer is missing bundled workflow file {required}. Reinstall malaria-amplicon-nf."
            ));
        }
    }

    Ok(runtime_dir)
}

fn shell_quote(value: &str) -> String {
    format!("'{}'", value.replace('\'', "'\"'\"'"))
}

#[cfg(windows)]
fn style_native_titlebar(window: &WebviewWindow) -> Result<(), String> {
    use std::{ffi::c_void, mem::size_of};
    use windows::Win32::Graphics::Dwm::{
        DwmSetWindowAttribute, DWMWA_BORDER_COLOR, DWMWA_CAPTION_COLOR, DWMWA_TEXT_COLOR,
    };

    let hwnd = window
        .hwnd()
        .map_err(|err| format!("Could not access the Windows title bar: {err}"))?;
    for (attribute, color) in [
        (DWMWA_CAPTION_COLOR, WINDOWS_CAPTION_BACKGROUND),
        (DWMWA_TEXT_COLOR, WINDOWS_CAPTION_TEXT),
        (DWMWA_BORDER_COLOR, WINDOWS_CAPTION_BORDER),
    ] {
        unsafe {
            DwmSetWindowAttribute(
                hwnd,
                attribute,
                (&color as *const u32).cast::<c_void>(),
                size_of::<u32>() as u32,
            )
        }
        .map_err(|err| format!("Could not style the Windows title bar: {err}"))?;
    }
    Ok(())
}

#[cfg(not(windows))]
fn style_native_titlebar(_window: &WebviewWindow) -> Result<(), String> {
    Ok(())
}

fn wsl_command_available_for(distro_name: &str, args: &[&str]) -> bool {
    let mut command = wsl_command_for(distro_name);
    command.args(args);
    hide_command_window(&mut command);
    let Ok(mut child) = command
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
    else {
        return false;
    };

    let started = Instant::now();
    while started.elapsed() < Duration::from_secs(20) {
        match child.try_wait() {
            Ok(Some(status)) => return status.success(),
            Ok(None) => thread::sleep(Duration::from_millis(100)),
            Err(_) => return false,
        }
    }
    let _ = child.kill();
    let _ = child.wait();
    false
}

fn native_simplseq_path() -> Option<String> {
    let home = env::var_os("HOME").or_else(|| env::var_os("USERPROFILE"))?;
    let mut path = PathBuf::from(home);
    path.push(".local");
    path.push("bin");
    path.push("simplseq");
    path.to_str().map(String::from)
}

fn native_runtime_ready_marker() -> Option<PathBuf> {
    let home = env::var_os("HOME").or_else(|| env::var_os("USERPROFILE"))?;
    let mut path = PathBuf::from(home);
    path.push(".local");
    path.push("share");
    path.push("simplseq");
    path.push("versions");
    path.push(RUNTIME_VERSION);
    path.push(".install-ready");
    Some(path)
}

fn wsl_runtime_installed_for(distro_name: &str) -> bool {
    wsl_command_available_for(
        distro_name,
        &[
            "bash",
            "-lc",
            "test -x ~/.local/bin/simplseq || command -v simplseq >/dev/null 2>&1",
        ],
    )
}

fn wsl_platform_installed() -> bool {
    let mut command = Command::new("wsl.exe");
    command.arg("--status");
    capture_command(&mut command, Duration::from_secs(5)).is_ok_and(|(success, _, _)| success)
}

fn wsl_distribution_ready() -> bool {
    registered_wsl_distributions()
        .ok()
        .and_then(|registered| find_runtime_wsl_distro(&registered))
        .is_some_and(|distro| wsl_distro_ready(&distro))
}

fn wsl_runtime_matches_expected_version_for(distro_name: &str) -> bool {
    let script = format!(
        "test -x ~/.local/bin/simplseq && test \"$(cat ~/.local/share/simplseq/versions/{}/.install-ready 2>/dev/null)\" = \"{}\"",
        RUNTIME_VERSION,
        RUNTIME_VERSION
    );
    wsl_command_available_for(distro_name, &["bash", "-lc", &script])
}

fn cleanup_stale_wsl_backends() {
    if !cfg!(target_os = "windows") {
        return;
    }
    let script = r#"
for pid_file in /tmp/malaria-amplicon-nf-backend-*.pid; do
  test -e "$pid_file" || continue
  pid="$(cat "$pid_file" 2>/dev/null || true)"
  case "$pid" in
    ''|*[!0-9]*) rm -f "$pid_file"; continue ;;
  esac
  if kill -0 "$pid" 2>/dev/null; then
    command_line="$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null || true)"
    case "$command_line" in
      *"python -m simplseq run"*)
        kill "$pid" 2>/dev/null || true
        for _ in 1 2 3 4 5; do
          kill -0 "$pid" 2>/dev/null || break
          sleep 0.1
        done
        kill -KILL "$pid" 2>/dev/null || true
        ;;
    esac
  fi
  rm -f "$pid_file"
done
"#;
    let mut command = runtime_wsl_command();
    command.args(["bash", "-lc", script]);
    hide_command_window(&mut command);
    let _ = command
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();
}

fn native_runtime_installed() -> bool {
    native_simplseq_path()
        .as_deref()
        .is_some_and(|path| command_available(path, &["--help"]))
        || command_available("simplseq", &["--help"])
}

fn native_runtime_matches_expected_version() -> bool {
    native_simplseq_path().is_some_and(|path| std::path::Path::new(&path).is_file())
        && native_runtime_ready_marker()
            .and_then(|path| std::fs::read_to_string(path).ok())
            .is_some_and(|version| version.trim() == RUNTIME_VERSION)
}

fn find_port() -> Result<u16, String> {
    let listener = TcpListener::bind("127.0.0.1:0").map_err(|err| err.to_string())?;
    let port = listener.local_addr().map_err(|err| err.to_string())?.port();
    drop(listener);
    Ok(port)
}

fn backend_healthy(port: u16) -> bool {
    let Ok(mut stream) = TcpStream::connect(("127.0.0.1", port)) else {
        return false;
    };
    let _ = stream.set_read_timeout(Some(Duration::from_secs(2)));
    let _ = stream.set_write_timeout(Some(Duration::from_secs(2)));
    let request =
        format!("GET /api/health HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\nConnection: close\r\n\r\n");
    if stream.write_all(request.as_bytes()).is_err() {
        return false;
    }
    let mut response = [0_u8; 8192];
    let Ok(size) = stream.read(&mut response) else {
        return false;
    };
    let text = String::from_utf8_lossy(&response[..size]);
    text.starts_with("HTTP/1.1 200")
        && text.contains("\"ok\":true")
        && text.contains("\"app\":\"malaria-amplicon-nf\"")
}

fn backend_url_healthy(url: &str) -> bool {
    Url::parse(url)
        .ok()
        .filter(|parsed| matches!(parsed.host_str(), Some("127.0.0.1") | Some("localhost")))
        .and_then(|parsed| parsed.port())
        .is_some_and(backend_healthy)
}

fn wait_for_backend(child: &mut Child, port: u16, timeout: Duration) -> Result<(), String> {
    let start = Instant::now();
    while start.elapsed() < timeout {
        if let Ok(Some(status)) = child.try_wait() {
            return Err(format!(
                "malaria-amplicon-nf backend exited before it became ready ({status})."
            ));
        }
        if backend_healthy(port) {
            return Ok(());
        }
        thread::sleep(Duration::from_millis(250));
    }
    Err("malaria-amplicon-nf backend did not become ready within 45 seconds.".into())
}

fn launcher_log_path() -> PathBuf {
    let base = env::var_os("LOCALAPPDATA")
        .or_else(|| env::var_os("XDG_STATE_HOME"))
        .map(PathBuf::from)
        .unwrap_or_else(env::temp_dir);
    base.join("malaria-amplicon-nf")
        .join("launcher-backend.log")
}

fn open_launcher_log() -> Result<File, String> {
    let path = launcher_log_path();
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|err| format!("Could not create launcher log directory: {err}"))?;
    }
    OpenOptions::new()
        .create(true)
        .append(true)
        .open(&path)
        .map_err(|err| format!("Could not open launcher log {}: {err}", path.display()))
}

fn expand_window_for_app(window: &WebviewWindow) -> Result<(), String> {
    window
        .set_min_size(Some(LogicalSize::new(
            APP_WINDOW_MIN_WIDTH,
            APP_WINDOW_MIN_HEIGHT,
        )))
        .map_err(|err| format!("Could not set app window minimum size: {err}"))?;
    window
        .set_size(LogicalSize::new(APP_WINDOW_WIDTH, APP_WINDOW_HEIGHT))
        .map_err(|err| format!("Could not resize app window: {err}"))?;
    window
        .center()
        .map_err(|err| format!("Could not center app window: {err}"))?;
    Ok(())
}

fn detect_runtime_info() -> RuntimeInfo {
    if cfg!(target_os = "windows") {
        let registered = match registered_wsl_distributions() {
            Ok(registered) => registered,
            Err(_) if !wsl_platform_installed() => {
                return RuntimeInfo {
                    mode: "wsl_missing".into(),
                    message: "Windows Subsystem for Linux is not installed.".into(),
                    detail: "Install WSL once, restart Windows if requested, then open malaria-amplicon-nf again.".into(),
                };
            }
            Err(error) => {
                return RuntimeInfo {
                    mode: "wsl_distribution_missing".into(),
                    message: "The WSL Linux environment could not be checked.".into(),
                    detail: error.into(),
                };
            }
        };
        let Some(distro_name) = find_runtime_wsl_distro(&registered) else {
            return RuntimeInfo {
                mode: "wsl_distribution_missing".into(),
                message: "Ubuntu is not installed in WSL.".into(),
                detail: "Install the standard Ubuntu distribution once, then malaria-amplicon-nf will install only its workflow packages.".into(),
            };
        };
        if !wsl_distro_ready(&distro_name) {
            RuntimeInfo {
                mode: "wsl_distribution_missing".into(),
                message: "The WSL Linux environment could not start.".into(),
                detail: "Restart WSL or Windows, then check again.".into(),
            }
        } else if wsl_runtime_matches_expected_version_for(&distro_name) {
            RuntimeInfo {
                mode: "wsl".into(),
                message: "WSL malaria-amplicon-nf runtime detected.".into(),
                detail: format!("malaria-amplicon-nf runtime {RUNTIME_VERSION} is ready.").into(),
            }
        } else if wsl_runtime_installed_for(&distro_name) {
            RuntimeInfo {
                mode: "outdated".into(),
                message: "Installed malaria-amplicon-nf runtime needs an update.".into(),
                detail: format!(
                    "Install the managed malaria-amplicon-nf runtime {RUNTIME_VERSION} inside WSL."
                )
                .into(),
            }
        } else {
            RuntimeInfo {
                mode: "missing".into(),
                message: "malaria-amplicon-nf runtime is not installed.".into(),
                detail: "Install the managed malaria-amplicon-nf runtime inside WSL.".into(),
            }
        }
    } else if native_runtime_matches_expected_version() {
        RuntimeInfo {
            mode: "native".into(),
            message: "malaria-amplicon-nf runtime detected.".into(),
            detail: format!("malaria-amplicon-nf runtime {RUNTIME_VERSION} is ready.").into(),
        }
    } else if native_runtime_installed() {
        RuntimeInfo {
            mode: "outdated".into(),
            message: "Installed malaria-amplicon-nf runtime needs an update.".into(),
            detail: format!("Install the managed malaria-amplicon-nf runtime {RUNTIME_VERSION} for this computer.").into(),
        }
    } else {
        RuntimeInfo {
            mode: "missing".into(),
            message: "malaria-amplicon-nf runtime is not installed.".into(),
            detail: "Install the managed malaria-amplicon-nf runtime for this computer.".into(),
        }
    }
}

#[tauri::command]
async fn detect_runtime(state: State<'_, BackendState>) -> Result<RuntimeInfo, String> {
    let info = tauri::async_runtime::spawn_blocking(detect_runtime_info)
        .await
        .map_err(|err| format!("Could not check the local runtime: {err}"))?;
    if let Ok(mut verified) = state.runtime_verified.lock() {
        *verified = matches!(info.mode.as_str(), "wsl" | "native");
    }
    Ok(info)
}

fn strip_terminal_codes(value: &str) -> String {
    let mut output = String::with_capacity(value.len());
    let mut characters = value.chars().peekable();
    while let Some(character) = characters.next() {
        if character == '\u{1b}' && characters.peek() == Some(&'[') {
            characters.next();
            for code in characters.by_ref() {
                if ('@'..='~').contains(&code) {
                    break;
                }
            }
            continue;
        }
        if character == '\t' || !character.is_control() {
            output.push(character);
        }
    }
    output.trim().to_string()
}

fn installer_progress_from_log_line(line: &str) -> Option<RuntimeProgress> {
    let (phase, title, detail) = match line.trim() {
        "== Installing bundled app files ==" => (
            "install-files",
            "Installing workflow files",
            "Copying the bundled malaria-amplicon-nf workflow into the managed runtime.",
        ),
        "== Downloading release files ==" => (
            "download-runtime",
            "Downloading workflow files",
            "Downloading the versioned malaria-amplicon-nf workflow bundle.",
        ),
        "== Verifying checksum ==" => (
            "verify-download",
            "Verifying workflow files",
            "Checking the downloaded workflow bundle before installation.",
        ),
        "== Installing app files ==" => (
            "install-files",
            "Installing workflow files",
            "Installing the verified malaria-amplicon-nf workflow bundle.",
        ),
        "== Installing micromamba ==" => (
            "install-micromamba",
            "Preparing the package manager",
            "Installing or reusing the managed Micromamba package manager.",
        ),
        "== Creating managed runtime ==" => (
            "create-runtime",
            "Installing the analysis environment",
            "Installing the pinned Nextflow, Python, R, and bioinformatics packages.",
        ),
        "== Installing downstream R analysis packages ==" => (
            "install-r-packages",
            "Installing analysis modules",
            "Preparing DINEMITES, Dcifer, and their pinned R dependencies.",
        ),
        "== Creating launcher ==" => (
            "create-launcher",
            "Creating the local launcher",
            "Connecting the managed environment to the desktop application.",
        ),
        "== Checking PATH ==" => (
            "check-path",
            "Finalizing the launcher",
            "Checking the managed command path.",
        ),
        "== Verifying malaria-amplicon-nf ==" => (
            "verify-runtime",
            "Verifying the installation",
            "Running final workflow and dependency checks.",
        ),
        "== Setup complete ==" => (
            "packages-ready",
            "Workflow packages ready",
            "Starting malaria-amplicon-nf...",
        ),
        _ => return None,
    };

    Some(RuntimeProgress {
        phase: phase.into(),
        title: title.into(),
        detail: detail.into(),
    })
}

fn spawn_runtime_log_reader<R>(
    reader: R,
    app: AppHandle,
    stream: &'static str,
) -> thread::JoinHandle<Vec<String>>
where
    R: Read + Send + 'static,
{
    thread::spawn(move || {
        let mut retained = Vec::new();
        for line in BufReader::new(reader).lines().map_while(Result::ok) {
            let line = strip_terminal_codes(&line);
            if line.is_empty() {
                continue;
            }
            if stream == "stdout" {
                if let Some(progress) = installer_progress_from_log_line(&line) {
                    let _ = app.emit("runtime-progress", progress);
                }
            }
            let _ = app.emit(
                "runtime-log",
                RuntimeLog {
                    line: line.clone(),
                    stream: stream.into(),
                },
            );
            retained.push(line);
            if retained.len() > 160 {
                retained.remove(0);
            }
        }
        retained
    })
}

#[cfg(test)]
mod installer_progress_tests {
    use super::installer_progress_from_log_line;

    #[test]
    fn maps_real_installer_checkpoints_to_progress_phases() {
        let runtime = installer_progress_from_log_line("== Creating managed runtime ==").unwrap();
        assert_eq!(runtime.phase, "create-runtime");

        let complete = installer_progress_from_log_line("== Setup complete ==").unwrap();
        assert_eq!(complete.phase, "packages-ready");
    }

    #[test]
    fn ignores_regular_package_output() {
        assert!(installer_progress_from_log_line("Downloading package metadata").is_none());
    }
}

fn wait_with_runtime_log(
    mut child: Child,
    app: &AppHandle,
) -> Result<(std::process::ExitStatus, Vec<String>), String> {
    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| "Could not capture installer output.".to_string())?;
    let stderr = child
        .stderr
        .take()
        .ok_or_else(|| "Could not capture installer errors.".to_string())?;
    let stdout_reader = spawn_runtime_log_reader(stdout, app.clone(), "stdout");
    let stderr_reader = spawn_runtime_log_reader(stderr, app.clone(), "stderr");
    let status = child
        .wait()
        .map_err(|err| format!("Could not wait for the installer: {err}"))?;
    let mut retained = stdout_reader.join().unwrap_or_default();
    retained.extend(stderr_reader.join().unwrap_or_default());
    if retained.len() > 160 {
        retained.drain(..retained.len() - 160);
    }
    Ok((status, retained))
}

#[tauri::command]
async fn install_runtime(app: AppHandle) -> Result<(), String> {
    tauri::async_runtime::spawn_blocking(move || {
        let _ = app.emit(
            "runtime-progress",
            RuntimeProgress {
                phase: "packages-install".into(),
                title: "Installing workflow packages".into(),
                detail: "Installing the pinned Nextflow, Python, and R analysis environment. This is a one-time step.".into(),
            },
        );
        let _ = app.emit(
            "runtime-log",
            RuntimeLog {
                line: "Starting the bundled workflow installer...".into(),
                stream: "status".into(),
            },
        );
        let result = run_runtime_installer(&app);
        if result.is_ok() {
            let _ = app.emit(
                "runtime-progress",
                RuntimeProgress {
                    phase: "packages-ready".into(),
                    title: "Workflow packages ready".into(),
                    detail: "Starting malaria-amplicon-nf...".into(),
                },
            );
        }
        result
    })
        .await
        .map_err(|err| format!("Runtime installer task failed: {err}"))?
}

#[tauri::command]
async fn install_wsl(app: AppHandle) -> Result<(), String> {
    tauri::async_runtime::spawn_blocking(move || {
        let progress_app = app.clone();
        let progress = move |phase: &str, title: &str, detail: &str| {
            let _ = progress_app.emit(
                "runtime-progress",
                RuntimeProgress {
                    phase: phase.into(),
                    title: title.into(),
                    detail: detail.into(),
                },
            );
        };
        run_wsl_installer(&app, &progress)
    })
    .await
    .map_err(|err| format!("WSL installer task failed: {err}"))?
}

fn run_wsl_installer(app: &AppHandle, progress: &dyn Fn(&str, &str, &str)) -> Result<(), String> {
    if !cfg!(target_os = "windows") {
        return Err("WSL installation is only available on Windows.".into());
    }
    progress(
        "environment-detect",
        "Checking Ubuntu",
        "Looking for an existing WSL Linux distribution...",
    );
    if wsl_distribution_ready() {
        return Ok(());
    }
    let platform_installed = wsl_platform_installed();
    if registered_wsl_distributions()
        .ok()
        .and_then(|registered| find_runtime_wsl_distro(&registered))
        .is_some()
    {
        return Err("The installed Ubuntu WSL distribution could not start. Restart WSL or Windows, then retry setup.".into());
    }

    progress(
        "ubuntu-install",
        "Installing Ubuntu for WSL",
        "Adding the standard Ubuntu distribution required by the workflow. This happens only when WSL has no Linux distribution.",
    );
    let mut command = if platform_installed {
        let mut command = Command::new("wsl.exe");
        command.args([
            "--install",
            "--distribution",
            DEFAULT_WSL_DISTRO_NAME,
            "--no-launch",
            "--web-download",
        ]);
        command
    } else {
        progress(
            "wsl-platform-install",
            "Enabling WSL",
            "Windows will request administrator approval. A restart may be required before Ubuntu can start.",
        );
        let script = format!(
            "$process = Start-Process -FilePath 'wsl.exe' -ArgumentList @('--install','--distribution','{}','--no-launch','--web-download') -Verb RunAs -Wait -PassThru; exit $process.ExitCode",
            DEFAULT_WSL_DISTRO_NAME
        );
        let mut command = Command::new("powershell.exe");
        command.args(["-NoProfile", "-NonInteractive", "-Command", &script]);
        command
    };
    hide_command_window(&mut command);
    let status = command
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|err| format!("Could not install Ubuntu for WSL: {err}"))?;
    let (status, output) = wait_with_runtime_log(status, app)?;
    if wsl_distro_ready(DEFAULT_WSL_DISTRO_NAME) {
        progress(
            "ubuntu-ready",
            "Ubuntu is ready",
            "Installing the malaria-amplicon-nf workflow packages next...",
        );
        Ok(())
    } else if status.success() && !platform_installed {
        Err("Windows needs a restart to finish enabling WSL. Restart the computer, then open malaria-amplicon-nf again.".into())
    } else {
        let detail = output.into_iter().rev().take(8).collect::<Vec<_>>();
        Err(format!(
            "Ubuntu installation did not complete ({status}). {}",
            detail.into_iter().rev().collect::<Vec<_>>().join(" ")
        ))
    }
}

fn run_runtime_installer(app: &AppHandle) -> Result<(), String> {
    let bundled_runtime = bundled_runtime_directory(app)?;
    let mut command = if cfg!(target_os = "windows") {
        let bundled_runtime_wsl = windows_path_for_wsl(&bundled_runtime)?;
        let mut command = runtime_wsl_command();
        command
            .arg("env")
            .arg(format!(
                "SIMPLSEQ_BUNDLED_RUNTIME_DIR={bundled_runtime_wsl}"
            ))
            .args(["bash", "-s", "--"]);
        command
    } else {
        let mut command = Command::new("bash");
        command.env("SIMPLSEQ_BUNDLED_RUNTIME_DIR", &bundled_runtime);
        command.args(["-s", "--"]);
        command
    };

    hide_command_window(&mut command);
    let mut child = command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|err| format!("Could not run malaria-amplicon-nf installer: {err}"))?;
    child
        .stdin
        .take()
        .ok_or_else(|| "Could not open the runtime installer input.".to_string())?
        .write_all(INSTALL_SCRIPT.as_bytes())
        .map_err(|err| format!("Could not send the bundled runtime installer to bash: {err}"))?;
    let (status, output) = wait_with_runtime_log(child, app)?;

    if status.success() {
        Ok(())
    } else {
        let summary = output
            .into_iter()
            .rev()
            .take(8)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>()
            .join(" ");
        Err(if summary.trim().is_empty() {
            "malaria-amplicon-nf runtime installer failed. Check the managed runtime install log for details.".into()
        } else {
            format!("malaria-amplicon-nf runtime installer failed: {summary}")
        })
    }
}

#[tauri::command]
fn navigate_to_backend(window: &WebviewWindow, url: &str) -> Result<(), String> {
    let parsed = Url::parse(url).map_err(|err| format!("Invalid backend URL: {err}"))?;
    expand_window_for_app(window)?;
    window
        .navigate(parsed)
        .map_err(|err| format!("Could not open malaria-amplicon-nf in the desktop window: {err}"))
}

#[tauri::command]
fn start_backend(
    window: WebviewWindow,
    state: State<BackendState>,
) -> Result<BackendLaunch, String> {
    if let Some(url) = state.url.lock().map_err(|err| err.to_string())?.clone() {
        if backend_url_healthy(&url) {
            navigate_to_backend(&window, &url)?;
            return Ok(BackendLaunch { url });
        }
        stop_backend_process(&state)?;
    }

    let runtime_verified = *state
        .runtime_verified
        .lock()
        .map_err(|err| err.to_string())?;
    if !runtime_verified {
        let runtime = detect_runtime_info();
        let runtime_verified = matches!(runtime.mode.as_str(), "wsl" | "native");
        *state
            .runtime_verified
            .lock()
            .map_err(|err| err.to_string())? = runtime_verified;
        if !runtime_verified {
            return Err(runtime.detail.to_string());
        }
    }

    cleanup_stale_wsl_backends();
    let picker_bridge_dir = state
        .picker_bridge_dir
        .lock()
        .map_err(|err| err.to_string())?
        .clone()
        .ok_or_else(|| "Desktop picker bridge is not ready.".to_string())?;
    let port = find_port()?;
    let url = format!("http://127.0.0.1:{port}");
    let mut command = if cfg!(target_os = "windows") {
        let pid_file = format!("/tmp/malaria-amplicon-nf-backend-{port}.pid");
        let picker_bridge_wsl = shell_quote(&windows_path_for_wsl(&picker_bridge_dir)?);
        let script = format!(
            "printf '%s\\n' $$ > {pid_file}; export SIMPLSEQ_PICKER_BRIDGE_DIR={picker_bridge_wsl}; exec ~/.local/bin/simplseq run --host 127.0.0.1 --port {port} --no-browser"
        );
        let mut command = runtime_wsl_command();
        command.args(["bash", "-lc", &script]);
        command
    } else {
        let command_path = native_simplseq_path().unwrap_or_else(|| "simplseq".into());
        let mut command = Command::new(command_path);
        command.args([
            "run",
            "--host",
            "127.0.0.1",
            "--port",
            &port.to_string(),
            "--no-browser",
        ]);
        command.env("SIMPLSEQ_PICKER_BRIDGE_DIR", &picker_bridge_dir);
        command
    };

    hide_command_window(&mut command);
    let mut backend_log = open_launcher_log()?;
    let backend_stderr = backend_log
        .try_clone()
        .map_err(|err| format!("Could not prepare launcher error log: {err}"))?;
    let _ = writeln!(
        backend_log,
        "\n=== starting backend on 127.0.0.1:{port} ({RUNTIME_VERSION}) ==="
    );
    let mut child = command
        .stdin(Stdio::null())
        .stdout(Stdio::from(backend_log))
        .stderr(Stdio::from(backend_stderr))
        .spawn()
        .map_err(|err| format!("Could not start malaria-amplicon-nf backend: {err}"))?;

    if let Err(error) = wait_for_backend(&mut child, port, Duration::from_secs(45)) {
        let _ = child.kill();
        let _ = child.wait();
        return Err(format!(
            "{error} Launcher log: {}",
            launcher_log_path().display()
        ));
    }

    *state.child.lock().map_err(|err| err.to_string())? = Some(child);
    *state.url.lock().map_err(|err| err.to_string())? = Some(url.clone());
    navigate_to_backend(&window, &url)?;
    Ok(BackendLaunch { url })
}

fn stop_backend_process(state: &BackendState) -> Result<(), String> {
    let url = state.url.lock().map_err(|err| err.to_string())?.take();
    if cfg!(target_os = "windows") {
        if let Some(port) = url
            .as_deref()
            .and_then(|value| Url::parse(value).ok())
            .and_then(|value| value.port())
        {
            let script = format!(
                "pid_file=/tmp/malaria-amplicon-nf-backend-{port}.pid; if test -s \"$pid_file\"; then kill \"$(cat \"$pid_file\")\" 2>/dev/null || true; fi; rm -f \"$pid_file\""
            );
            let mut command = runtime_wsl_command();
            command.args(["bash", "-lc", &script]);
            hide_command_window(&mut command);
            let _ = command
                .stdin(Stdio::null())
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status();
        }
    }
    if let Some(mut child) = state.child.lock().map_err(|err| err.to_string())?.take() {
        let _ = child.kill();
        let _ = child.wait();
    }
    Ok(())
}

#[tauri::command]
fn stop_backend(state: State<BackendState>) -> Result<(), String> {
    stop_backend_process(&state)
}

pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .manage(BackendState::default())
        .invoke_handler(tauri::generate_handler![
            detect_runtime,
            install_wsl,
            install_runtime,
            start_backend,
            stop_backend
        ])
        .setup(|app| {
            let picker_bridge =
                start_picker_bridge(app.handle().clone()).map_err(std::io::Error::other)?;
            let state = app.state::<BackendState>();
            *state
                .picker_bridge_dir
                .lock()
                .map_err(|err| std::io::Error::other(err.to_string()))? = Some(picker_bridge);
            let window = app.get_webview_window("main").expect("main window");
            window.set_title("malaria-amplicon-nf")?;
            style_native_titlebar(&window).map_err(std::io::Error::other)?;
            Ok(())
        })
        .on_window_event(|window, event| {
            if matches!(event, WindowEvent::CloseRequested { .. }) {
                let state = window.state::<BackendState>();
                let _ = stop_backend_process(&state);
            }
        })
        .run(tauri::generate_context!())
        .expect("error while running malaria-amplicon-nf desktop app");
}
