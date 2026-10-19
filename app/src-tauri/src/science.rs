//! Science-core process management.
//!
//! The app ships the Python science environment as a single archive
//! (`runtime/env.tar.gz`, built with conda-pack) plus the `athanor_core`
//! package (`runtime/core`). On first launch, or after an update that
//! changes the engines, the archive is extracted into the user's local app
//! data folder and relocated by `runtime/unpack_env.py`. The core then runs as a child
//! process speaking JSON-RPC over stdin/stdout: no network port is opened.

use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};
use std::sync::Mutex;

use serde::Serialize;
use serde_json::{json, Value};
use tauri::{AppHandle, Manager};

#[cfg(windows)]
use std::os::windows::process::CommandExt;
#[cfg(windows)]
const CREATE_NO_WINDOW: u32 = 0x0800_0000;

#[derive(Clone, Serialize)]
pub struct Status {
    pub stage: String,
    pub message: String,
}

impl Status {
    pub fn new(stage: &str, message: &str) -> Self {
        Status { stage: stage.into(), message: message.into() }
    }
}

#[derive(Clone)]
pub struct LaunchSpec {
    /// Program and arguments that start `athanor_core.server` speaking JSON-RPC on stdio.
    pub program: PathBuf,
    pub args: Vec<String>,
    pub env: HashMap<String, String>,
    pub cwd: PathBuf,
    /// Extra command run after killing the child (e.g. stop the process inside WSL).
    pub stop: Option<(PathBuf, Vec<String>)>,
}

impl LaunchSpec {
    pub fn python(python: PathBuf, env: HashMap<String, String>, cwd: PathBuf) -> Self {
        LaunchSpec { program: python, args: vec!["-m".into(), "athanor_core.server".into()], env, cwd, stop: None }
    }
}

struct Io {
    stdin: ChildStdin,
    stdout: BufReader<ChildStdout>,
    next_id: u64,
}

pub struct CoreHandle {
    io: Mutex<Option<Io>>,
    child: Mutex<Option<Child>>,
    spec: Mutex<Option<LaunchSpec>>,
}

impl CoreHandle {
    pub fn new() -> Self {
        CoreHandle { io: Mutex::new(None), child: Mutex::new(None), spec: Mutex::new(None) }
    }

    pub fn start(&self, spec: LaunchSpec) -> Result<(), String> {
        let mut cmd = Command::new(&spec.program);
        cmd.args(&spec.args)
            .current_dir(&spec.cwd)
            .envs(&spec.env)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null());
        #[cfg(windows)]
        cmd.creation_flags(CREATE_NO_WINDOW);
        let mut child = cmd.spawn().map_err(|e| format!("could not start the science core: {e}"))?;
        let stdin = child.stdin.take().ok_or("no stdin")?;
        let mut stdout = BufReader::new(child.stdout.take().ok_or("no stdout")?);
        // First line is the "ready" notification.
        let mut first = String::new();
        stdout
            .read_line(&mut first)
            .map_err(|e| format!("science core did not start: {e}"))?;
        if first.trim().is_empty() {
            return Err("science core exited during startup".into());
        }
        *self.child.lock().unwrap() = Some(child);
        *self.io.lock().unwrap() = Some(Io { stdin, stdout, next_id: 1 });
        *self.spec.lock().unwrap() = Some(spec);
        Ok(())
    }

    pub fn call(&self, method: &str, params: Value) -> Result<Value, String> {
        let mut guard = self.io.lock().unwrap();
        let io = guard.as_mut().ok_or("science core is not running yet")?;
        let id = io.next_id;
        io.next_id += 1;
        let req = json!({"jsonrpc": "2.0", "id": id, "method": method, "params": params});
        writeln!(io.stdin, "{req}").map_err(|e| format!("science core is not reachable: {e}"))?;
        io.stdin.flush().map_err(|e| e.to_string())?;
        loop {
            let mut line = String::new();
            let n = io.stdout.read_line(&mut line).map_err(|e| e.to_string())?;
            if n == 0 {
                return Err("the science core stopped".into());
            }
            let msg: Value = match serde_json::from_str(line.trim()) {
                Ok(v) => v,
                Err(_) => continue,
            };
            if msg.get("id").and_then(Value::as_u64) != Some(id) {
                continue; // notification or stale response
            }
            if let Some(err) = msg.get("error") {
                return Err(err.get("message").and_then(Value::as_str).unwrap_or("error").to_string());
            }
            return Ok(msg.get("result").cloned().unwrap_or(Value::Null));
        }
    }

    fn run_stop_hook(&self) {
        let hook = self.spec.lock().unwrap().as_ref().and_then(|s| s.stop.clone());
        if let Some((prog, args)) = hook {
            let mut cmd = Command::new(prog);
            cmd.args(args).stdin(Stdio::null()).stdout(Stdio::null()).stderr(Stdio::null());
            #[cfg(windows)]
            cmd.creation_flags(CREATE_NO_WINDOW);
            let _ = cmd.status();
        }
    }

    pub fn is_running(&self) -> bool {
        self.io.lock().map(|g| g.is_some()).unwrap_or(false)
    }

    pub fn restart(&self) -> Result<(), String> {
        if let Some(mut c) = self.child.lock().unwrap().take() {
            let _ = c.kill();
            let _ = c.wait();
        }
        self.run_stop_hook();
        // Any in-flight call now sees EOF and releases the io lock.
        let spec = self.spec.lock().unwrap().clone().ok_or("science core was never started")?;
        *self.io.lock().unwrap() = None;
        self.start(spec)
    }

    pub fn shutdown(&self) {
        let had_child = if let Some(mut c) = self.child.lock().unwrap().take() {
            let _ = c.kill();
            let _ = c.wait();
            true
        } else {
            false
        };
        if had_child {
            self.run_stop_hook();
        }
        // An in-flight call holds the lock until it sees EOF; it clears itself then.
        if let Ok(mut g) = self.io.try_lock() {
            *g = None;
        }
    }
}

/// Work out how to launch the core. Development: set ATHANOR_PYTHON (and
/// optionally ATHANOR_CORE_DIR). Installed app: unpack the bundled runtime.
pub fn launch_spec(app: &AppHandle, emit: &dyn Fn(Status)) -> Result<LaunchSpec, String> {
    if let Ok(py) = std::env::var("ATHANOR_PYTHON") {
        let core_dir = std::env::var("ATHANOR_CORE_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|_| PathBuf::from("../../core"));
        let mut env = HashMap::new();
        env.insert("PYTHONPATH".into(), core_dir.to_string_lossy().into_owned());
        env.insert("PYTHONUTF8".into(), "1".into());
        return Ok(LaunchSpec::python(py.into(), env, core_dir));
    }

    let res = app.path().resource_dir().map_err(|e| e.to_string())?.join("runtime");
    let env_id = fs::read_to_string(res.join("env_id.txt"))
        .map_err(|e| format!("installation is incomplete (env_id.txt missing): {e}"))?
        .trim()
        .to_string();
    let data = app.path().app_local_data_dir().map_err(|e| e.to_string())?;
    let envs = data.join("envs");
    let env_dir = envs.join(&env_id);
    let marker = env_dir.join(".athanor-ready");

    if !marker.exists() {
        emit(Status::new("preparing", "Preparing science engines (first launch after install or update; about a minute)"));
        let _ = fs::remove_dir_all(&env_dir);
        fs::create_dir_all(&env_dir).map_err(|e| format!("could not create {}: {e}", env_dir.display()))?;
        extract(&res.join("env.tar.gz"), &env_dir)?;
        emit(Status::new("preparing", "Relocating engine paths"));
        unpack(&env_dir, &res.join("unpack_env.py"))?;
        fs::write(&marker, env_id.as_bytes()).map_err(|e| e.to_string())?;
        // Remove environments left by previous versions.
        if let Ok(rd) = fs::read_dir(&envs) {
            for e in rd.flatten() {
                if e.file_name().to_string_lossy() != env_id {
                    let _ = fs::remove_dir_all(e.path());
                }
            }
        }
    }

    let core_dir = res.join("core");
    let scratch = data.join("scratch");
    let _ = fs::create_dir_all(&scratch);
    Ok(LaunchSpec::python(python_path(&env_dir), runtime_env(&env_dir, &core_dir, &scratch), scratch))
}

fn python_path(env: &Path) -> PathBuf {
    if cfg!(windows) { env.join("python.exe") } else { env.join("bin").join("python") }
}

/// The variables `conda activate` would set, so DLLs and data files resolve
/// without a shell activation step.
fn runtime_env(env: &Path, core: &Path, scratch: &Path) -> HashMap<String, String> {
    let mut m = HashMap::new();
    let sep = if cfg!(windows) { ";" } else { ":" };
    let mut parts: Vec<PathBuf> = if cfg!(windows) {
        vec![
            env.to_path_buf(),
            env.join("Library").join("mingw-w64").join("bin"),
            env.join("Library").join("usr").join("bin"),
            env.join("Library").join("bin"),
            env.join("Scripts"),
            env.join("bin"),
        ]
    } else {
        vec![env.join("bin")]
    };
    if let Ok(p) = std::env::var("PATH") {
        parts.extend(std::env::split_paths(&p));
    }
    m.insert("PATH".into(), parts.iter().map(|p| p.to_string_lossy().into_owned()).collect::<Vec<_>>().join(sep));
    m.insert("CONDA_PREFIX".into(), env.to_string_lossy().into_owned());
    m.insert("PYTHONPATH".into(), core.to_string_lossy().into_owned());
    m.insert("PYTHONNOUSERSITE".into(), "1".into());
    m.insert("PYTHONUTF8".into(), "1".into());
    m.insert("PSI_SCRATCH".into(), scratch.to_string_lossy().into_owned());
    m
}

fn extract(archive: &Path, dest: &Path) -> Result<(), String> {
    // Windows 10 1803+ and every Linux/macOS ship a `tar` that reads .tar.gz.
    let mut cmd = Command::new("tar");
    cmd.arg("-xzf").arg(archive).arg("-C").arg(dest);
    #[cfg(windows)]
    cmd.creation_flags(CREATE_NO_WINDOW);
    let out = cmd.output().map_err(|e| format!("could not run tar: {e}"))?;
    if !out.status.success() {
        return Err(format!("extracting the science engines failed: {}", String::from_utf8_lossy(&out.stderr)));
    }
    Ok(())
}

fn unpack(env: &Path, unpacker: &Path) -> Result<(), String> {
    // Rewrites the build machine's paths baked into the packed environment.
    // unpack_env.py replaces conda-pack's conda-unpack, which fails on Windows
    // by opening DLLs its own interpreter has loaded for writing.
    let mut cmd = Command::new(python_path(env));
    cmd.arg("-I").arg("-S").arg(unpacker).arg(env);
    cmd.envs(runtime_env(env, env, env));
    #[cfg(windows)]
    cmd.creation_flags(CREATE_NO_WINDOW);
    let out = cmd.output().map_err(|e| format!("could not run the relocation step: {e}"))?;
    if !out.status.success() {
        return Err(format!(
            "relocating the science engines failed: {}{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    Ok(())
}
