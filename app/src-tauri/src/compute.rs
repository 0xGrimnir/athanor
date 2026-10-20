//! Optional GPU compute pack.
//!
//! GPU4PySCF runs only on Linux, so on Windows the pack is a private WSL2
//! distribution ("AthanorCompute") that this module installs and manages:
//!
//! 1. Make sure WSL2 is enabled (one administrator prompt if it is not).
//! 2. Download the pack's parts from the public release feed with Windows'
//!    own curl.exe, verifying each part's SHA-256 against the manifest.
//! 3. Import the root filesystem with `wsl --import`.
//! 4. Run the GPU self-test (GPU vs CPU at identical settings). The GPU tier
//!    is offered only if it passes.
//! 5. Run a second science core inside the distribution. It reads the same
//!    `athanor_core` package the app ships, through WSL's /mnt/<drive> view,
//!    so the GPU core always matches the app version.

use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};
use std::time::Duration;

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use tauri::{AppHandle, Emitter, Manager};

use crate::science::{CoreHandle, LaunchSpec};

#[cfg(windows)]
use std::os::windows::process::CommandExt;
#[cfg(windows)]
const CREATE_NO_WINDOW: u32 = 0x0800_0000;

/// Pack version this app build expects. Bump together with a new pack release.
pub const PACK_VERSION: &str = "1";
pub const DISTRO: &str = "AthanorCompute";
const FEED: &str = "https://github.com/0xGrimnir/athanor-releases/releases/download";

#[derive(Clone, Serialize, Default)]
pub struct ComputeStatus {
    /// unsupported | off | needs-restart | downloading | installing | testing | ready | failed
    pub stage: String,
    pub message: String,
    pub progress: Option<f64>,
    pub gpu: Option<String>,
    pub selftest: Option<Value>,
}

#[derive(Deserialize)]
struct Part {
    name: String,
    size: u64,
    sha256: String,
}

#[derive(Deserialize)]
struct Manifest {
    distro: String,
    python: String,
    parts: Vec<Part>,
    sha256: String,
    size: u64,
}

#[derive(Serialize, Deserialize)]
struct Installed {
    version: String,
    python: String,
    selftest: Value,
}

pub struct Compute {
    pub core: Arc<CoreHandle>,
    status: Mutex<ComputeStatus>,
    busy: Mutex<bool>,
}

impl Compute {
    pub fn new() -> Self {
        Compute { core: Arc::new(CoreHandle::new()), status: Mutex::new(ComputeStatus::default()), busy: Mutex::new(false) }
    }

    pub fn status(&self) -> ComputeStatus {
        self.status.lock().unwrap().clone()
    }

    fn set(&self, app: &AppHandle, st: ComputeStatus) {
        *self.status.lock().unwrap() = st.clone();
        let _ = app.emit("compute-status", &st);
    }

    fn update(&self, app: &AppHandle, stage: &str, message: &str, progress: Option<f64>) {
        let mut st = self.status();
        st.stage = stage.into();
        st.message = message.into();
        st.progress = progress;
        self.set(app, st);
    }

    /// Called at startup: report the current state and start the GPU core if
    /// the pack is installed and its self-test passed.
    pub fn init(&self, app: &AppHandle) {
        if !cfg!(windows) {
            self.set(app, ComputeStatus { stage: "unsupported".into(), message: "The GPU compute pack is available on Windows.".into(), ..Default::default() });
            return;
        }
        let gpu = nvidia_gpu();
        let base = ComputeStatus { gpu: gpu.clone(), ..Default::default() };
        let installed = read_installed(app);
        match installed {
            Some(inst) if inst.version == PACK_VERSION && distro_exists() && passed(&inst.selftest) => {
                self.set(app, ComputeStatus { stage: "testing".into(), message: "Starting GPU core".into(), selftest: Some(inst.selftest.clone()), ..base });
                match self.start_core(app, &inst.python) {
                    Ok(()) => self.update(app, "ready", "GPU compute ready", None),
                    Err(e) => self.update(app, "failed", &format!("GPU core did not start: {e}"), None),
                }
            }
            Some(inst) if inst.version != PACK_VERSION => self.set(app, ComputeStatus {
                stage: "off".into(), message: format!("A new compute pack ({PACK_VERSION}) is required by this version of Athanor."), ..base }),
            Some(inst) if !passed(&inst.selftest) => self.set(app, ComputeStatus {
                stage: "failed".into(), message: "The GPU self-test did not pass; the GPU tier is off.".into(), selftest: Some(inst.selftest), ..base }),
            _ => self.set(app, ComputeStatus {
                stage: "off".into(),
                message: match &gpu { Some(g) => format!("GPU compute is not installed. Detected: {g}"), None => "No NVIDIA GPU driver detected.".into() },
                ..base }),
        }
    }

    fn start_core(&self, app: &AppHandle, python: &str) -> Result<(), String> {
        let core_dir = app.path().resource_dir().map_err(|e| e.to_string())?.join("runtime").join("core");
        let linux_core = to_wsl_path(&core_dir)?;
        let spec = LaunchSpec {
            program: PathBuf::from("wsl.exe"),
            args: vec![
                "-d".into(), DISTRO.into(), "--cd".into(), "/root".into(), "--exec".into(), "/usr/bin/env".into(),
                format!("PYTHONPATH={linux_core}"), "PYTHONUNBUFFERED=1".into(), "PYTHONUTF8=1".into(),
                "LD_LIBRARY_PATH=/usr/lib/wsl/lib".into(),
                python.into(), "-m".into(), "athanor_core.server".into(),
            ],
            env: [("WSL_UTF8".to_string(), "1".to_string())].into_iter().collect(),
            cwd: std::env::temp_dir(),
            stop: Some((PathBuf::from("wsl.exe"), vec!["-d".into(), DISTRO.into(), "--exec".into(),
                        "pkill".into(), "-f".into(), "athanor_core.server".into()])),
        };
        self.core.start(spec)
    }

    /// Full installation. Runs on a background thread; progress goes out as
    /// "compute-status" events.
    pub fn enable(self: &Arc<Self>, app: AppHandle) -> Result<(), String> {
        {
            let mut b = self.busy.lock().unwrap();
            if *b {
                return Err("the compute pack is already being installed".into());
            }
            *b = true;
        }
        let me = self.clone();
        std::thread::spawn(move || {
            let r = me.enable_inner(&app);
            if let Err(e) = r {
                let stage = if e.starts_with("RESTART:") { "needs-restart" } else { "failed" };
                me.update(&app, stage, e.trim_start_matches("RESTART:"), None);
            }
            *me.busy.lock().unwrap() = false;
        });
        Ok(())
    }

    fn enable_inner(&self, app: &AppHandle) -> Result<(), String> {
        if !cfg!(windows) {
            return Err("The GPU compute pack is available on Windows.".into());
        }
        if nvidia_gpu().is_none() {
            return Err("No NVIDIA GPU driver detected. Install the current NVIDIA driver, then try again.".into());
        }
        // 1. WSL2
        self.update(app, "installing", "Checking WSL2", None);
        if !wsl_ready() {
            self.update(app, "installing", "Enabling WSL2: approve the Windows administrator prompt", None);
            let ok = elevate_wsl_install()?;
            if !ok || !wsl_ready() {
                return Err("RESTART:WSL2 was enabled. Restart Windows, then click Enable GPU compute again to continue.".into());
            }
        }
        // 2-3. Download and import, unless this pack version is already imported
        let dir = app.path().app_local_data_dir().map_err(|e| e.to_string())?.join("compute");
        let python = match existing_pack(&dir) {
            Some(py) => {
                self.update(app, "testing", "Compute pack already installed; re-running the self-test", None);
                py
            }
            None => self.install_pack(app, &dir)?,
        };
        // 4. Self-test on this machine's GPU
        self.update(app, "testing", "Starting GPU core", None);
        self.start_core(app, &python)?;
        self.update(app, "testing", "Running GPU self-test (GPU vs CPU on water and caffeine)", None);
        let report = self.core.call("gpu_selftest", json!({}))?;
        let inst = Installed { version: PACK_VERSION.into(), python: python.clone(), selftest: report.clone() };
        fs::write(dir.join("installed.json"), serde_json::to_vec_pretty(&inst).unwrap()).map_err(|e| e.to_string())?;
        let mut st = self.status();
        st.selftest = Some(report.clone());
        st.progress = None;
        if passed(&report) {
            st.stage = "ready".into();
            st.message = "GPU compute ready: self-test passed".into();
        } else {
            self.core.shutdown();
            st.stage = "failed".into();
            st.message = report.get("error").and_then(Value::as_str)
                .map(|e| format!("GPU self-test failed: {e}"))
                .unwrap_or_else(|| "GPU self-test failed: GPU and CPU results disagree; the GPU tier stays off".into());
        }
        self.set(app, st);
        Ok(())
    }

    /// Download, verify and import the pack; returns the interpreter path inside the distribution.
    fn install_pack(&self, app: &AppHandle, dir: &Path) -> Result<String, String> {
        // 2. Download
        let dl = dir.join("download");
        fs::create_dir_all(&dl).map_err(|e| e.to_string())?;
        let base = format!("{FEED}/compute-{PACK_VERSION}");
        self.update(app, "downloading", "Fetching compute pack manifest", Some(0.0));
        let man_path = dl.join("manifest.json");
        curl(&format!("{base}/manifest.json"), &man_path, false)?;
        let manifest: Manifest = serde_json::from_slice(&fs::read(&man_path).map_err(|e| e.to_string())?)
            .map_err(|e| format!("manifest is not valid: {e}"))?;
        if manifest.distro != DISTRO {
            return Err(format!("unexpected distribution name in manifest: {}", manifest.distro));
        }
        let mut done: u64 = 0;
        for p in &manifest.parts {
            let dest = dl.join(&p.name);
            let already = fs::metadata(&dest).map(|m| m.len()).unwrap_or(0);
            if already == p.size && sha256_file(&dest)? == p.sha256 {
                done += p.size;
                continue;
            }
            if already > p.size {
                let _ = fs::remove_file(&dest);
            }
            let total = manifest.size as f64;
            let msg = format!("Downloading compute pack ({:.1} GB)", manifest.size as f64 / 1e9);
            let mut child = curl_spawn(&format!("{base}/{}", p.name), &dest)?;
            loop {
                if let Some(status) = child.try_wait().map_err(|e| e.to_string())? {
                    if !status.success() {
                        return Err(format!("download of {} failed (curl exit {:?}); click Enable again to resume", p.name, status.code()));
                    }
                    break;
                }
                let now = fs::metadata(&dest).map(|m| m.len()).unwrap_or(0);
                self.update(app, "downloading", &msg, Some((done + now) as f64 / total));
                std::thread::sleep(Duration::from_millis(700));
            }
            if sha256_file(&dest)? != p.sha256 {
                let _ = fs::remove_file(&dest);
                return Err(format!("{} failed its SHA-256 check and was deleted; click Enable again", p.name));
            }
            done += p.size;
        }
        // 3. Assemble, verify, import
        self.update(app, "installing", "Verifying compute pack", Some(1.0));
        let tarball = dl.join("athanor-compute.tar.gz");
        {
            let mut out = File::create(&tarball).map_err(|e| e.to_string())?;
            let mut h = Sha256::new();
            let mut buf = vec![0u8; 1 << 22];
            for p in &manifest.parts {
                let mut f = File::open(dl.join(&p.name)).map_err(|e| e.to_string())?;
                loop {
                    let n = f.read(&mut buf).map_err(|e| e.to_string())?;
                    if n == 0 {
                        break;
                    }
                    h.update(&buf[..n]);
                    out.write_all(&buf[..n]).map_err(|e| e.to_string())?;
                }
            }
            if hex(&h.finalize()) != manifest.sha256 {
                return Err("assembled compute pack failed its SHA-256 check".into());
            }
        }
        for p in &manifest.parts {
            let _ = fs::remove_file(dl.join(&p.name));
        }
        self.core.shutdown();
        if distro_exists() {
            self.update(app, "installing", "Removing previous compute pack", None);
            let _ = wsl(&["--unregister", DISTRO]);
        }
        let target = dir.join("distro");
        let _ = fs::remove_dir_all(&target);
        fs::create_dir_all(&target).map_err(|e| e.to_string())?;
        self.update(app, "installing", "Installing into WSL2 (a few minutes)", None);
        let (ok, out) = wsl(&["--import", DISTRO, &target.to_string_lossy(), &tarball.to_string_lossy(), "--version", "2"])?;
        let _ = fs::remove_file(&tarball);
        if !ok {
            return Err(format!("wsl --import failed: {out}"));
        }
        let marker = Installed { version: PACK_VERSION.into(), python: manifest.python.clone(), selftest: Value::Null };
        fs::write(dir.join("pack.json"), serde_json::to_vec_pretty(&marker).unwrap()).map_err(|e| e.to_string())?;
        Ok(manifest.python)
    }

    /// Unregister the distribution and delete everything the pack installed.
    pub fn remove(&self, app: &AppHandle) -> Result<(), String> {
        self.core.shutdown();
        if distro_exists() {
            let (ok, out) = wsl(&["--unregister", DISTRO])?;
            if !ok {
                return Err(format!("wsl --unregister failed: {out}"));
            }
        }
        let dir = app.path().app_local_data_dir().map_err(|e| e.to_string())?.join("compute");
        let _ = fs::remove_dir_all(&dir);
        self.init(app);
        Ok(())
    }
}

/// The interpreter path if this pack version is already imported and usable.
/// Installs made by 0.2.0-0.2.2 have no pack.json marker; pack 1 is the only
/// version those could have imported, so a working interpreter identifies it.
fn existing_pack(dir: &Path) -> Option<String> {
    if !distro_exists() {
        return None;
    }
    let python = match fs::read(dir.join("pack.json")).ok().and_then(|b| serde_json::from_slice::<Installed>(&b).ok()) {
        Some(m) if m.version == PACK_VERSION => m.python,
        Some(_) => return None,
        None if PACK_VERSION == "1" => "/opt/athanor/venv/bin/python".to_string(),
        None => return None,
    };
    // find_spec checks installation without importing gpu4pyscf (whose import initialises CUDA).
    let probe = "import importlib.util as u, sys; sys.exit(0 if u.find_spec('pyscf') and u.find_spec('gpu4pyscf') else 3)";
    match wsl(&["-d", DISTRO, "--exec", python.as_str(), "-c", probe]) {
        Ok((true, _)) => Some(python),
        _ => None,
    }
}

fn passed(report: &Value) -> bool {
    report.get("passed").and_then(Value::as_bool).unwrap_or(false)
}

fn read_installed(app: &AppHandle) -> Option<Installed> {
    let p = app.path().app_local_data_dir().ok()?.join("compute").join("installed.json");
    serde_json::from_slice(&fs::read(p).ok()?).ok()
}

fn hidden(cmd: &mut Command) -> &mut Command {
    #[cfg(windows)]
    cmd.creation_flags(CREATE_NO_WINDOW);
    cmd
}

/// Run wsl.exe with UTF-8 output; returns (success, combined output).
fn wsl(args: &[&str]) -> Result<(bool, String), String> {
    let mut cmd = Command::new("wsl.exe");
    cmd.args(args).env("WSL_UTF8", "1").stdin(Stdio::null());
    let out = hidden(&mut cmd).output().map_err(|e| format!("could not run wsl.exe: {e}"))?;
    let text = format!("{}{}", String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr));
    Ok((out.status.success(), text.replace('\0', "").trim().to_string()))
}

/// The Store version of WSL (the one with `--version`) is what supports
/// `--install --no-distribution` and WSL2 imports; the legacy inbox stub is not enough.
fn wsl_ready() -> bool {
    matches!(wsl(&["--version"]), Ok((true, _)))
}

fn distro_exists() -> bool {
    match wsl(&["--list", "--quiet"]) {
        Ok((true, out)) => out.lines().any(|l| l.trim() == DISTRO),
        _ => false,
    }
}

/// `wsl --install --no-distribution` needs administrator rights: ask through UAC.
fn elevate_wsl_install() -> Result<bool, String> {
    let script = "try { $p = Start-Process -FilePath 'wsl.exe' -ArgumentList '--install','--no-distribution' \
                  -Verb RunAs -Wait -PassThru; exit $p.ExitCode } catch { exit 1223 }";
    let mut cmd = Command::new("powershell.exe");
    cmd.args(["-NoProfile", "-NonInteractive", "-Command", script]);
    let st = hidden(&mut cmd).status().map_err(|e| format!("could not start PowerShell: {e}"))?;
    match st.code() {
        Some(1223) => Err("WSL2 was not enabled: the administrator prompt was declined.".into()),
        Some(0) => Ok(true),
        _ => Ok(false),
    }
}

/// Name and driver of the first NVIDIA GPU, from the driver's own nvidia-smi.
fn nvidia_gpu() -> Option<String> {
    let mut cmd = Command::new("nvidia-smi");
    cmd.args(["--query-gpu=name,driver_version", "--format=csv,noheader"]);
    let out = hidden(&mut cmd).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8_lossy(&out.stdout);
    let line = s.lines().next()?.trim();
    let mut it = line.splitn(2, ',');
    let name = it.next()?.trim();
    let drv = it.next().map(str::trim).unwrap_or("");
    Some(format!("{name}, driver {drv}"))
}

fn curl_command(url: &str, dest: &Path, resume: bool) -> Command {
    // curl.exe ships with Windows 10 (1803+) and 11.
    let mut cmd = Command::new("curl.exe");
    cmd.args(["-L", "--fail", "-sS", "--retry", "3", "--retry-delay", "2"]);
    if resume {
        cmd.args(["-C", "-"]);
    }
    cmd.arg("-o").arg(dest).arg(url).stdin(Stdio::null()).stdout(Stdio::null()).stderr(Stdio::piped());
    hidden(&mut cmd);
    cmd
}

fn curl(url: &str, dest: &Path, resume: bool) -> Result<(), String> {
    let out = curl_command(url, dest, resume).output().map_err(|e| format!("could not run curl.exe: {e}"))?;
    if !out.status.success() {
        return Err(format!("download failed ({url}): {}", String::from_utf8_lossy(&out.stderr).trim()));
    }
    Ok(())
}

fn curl_spawn(url: &str, dest: &Path) -> Result<std::process::Child, String> {
    curl_command(url, dest, true).spawn().map_err(|e| format!("could not run curl.exe: {e}"))
}

fn sha256_file(p: &Path) -> Result<String, String> {
    let mut f = File::open(p).map_err(|e| e.to_string())?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 22];
    loop {
        let n = f.read(&mut buf).map_err(|e| e.to_string())?;
        if n == 0 {
            break;
        }
        h.update(&buf[..n]);
    }
    Ok(hex(&h.finalize()))
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}

/// C:\Users\me\x -> /mnt/c/Users/me/x (WSL's default automount).
pub fn to_wsl_path(p: &Path) -> Result<String, String> {
    let s = p.to_string_lossy().replace('/', "\\");
    let s = s.strip_prefix("\\\\?\\").unwrap_or(&s).to_string();
    let mut chars = s.chars();
    match (chars.next(), chars.next()) {
        (Some(d), Some(':')) if d.is_ascii_alphabetic() => {
            let rest = &s[2..];
            Ok(format!("/mnt/{}{}", d.to_ascii_lowercase(), rest.replace('\\', "/")))
        }
        _ => Err(format!("cannot map {s} into WSL (expected a drive-letter path)")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wsl_paths() {
        assert_eq!(to_wsl_path(Path::new(r"C:\Users\Me\AppData\Local\Athanor\runtime\core")).unwrap(),
                   "/mnt/c/Users/Me/AppData/Local/Athanor/runtime/core");
        assert_eq!(to_wsl_path(Path::new(r"\\?\D:\x y\z")).unwrap(), "/mnt/d/x y/z");
        assert!(to_wsl_path(Path::new(r"\\server\share")).is_err());
    }
}
