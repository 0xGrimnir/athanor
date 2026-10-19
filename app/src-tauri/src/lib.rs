mod compute;
mod science;

use std::sync::Arc;

use serde::Serialize;
use serde_json::Value;
use tauri::{Emitter, Manager, State};
use tauri_plugin_updater::UpdaterExt;

use crate::compute::{Compute, ComputeStatus};
use crate::science::{CoreHandle, Status};

/// Forward one JSON-RPC call to the science core and return its `result`.
/// Requests for the GPU tier go to the science core inside the compute pack.
#[tauri::command]
async fn core_call(
    core: State<'_, Arc<CoreHandle>>,
    compute: State<'_, Arc<Compute>>,
    method: String,
    params: Value,
) -> Result<Value, String> {
    let gpu = params.get("tier").and_then(Value::as_str) == Some("dft-gpu");
    let core = if gpu {
        let c = compute.inner().core.clone();
        if !c.is_running() {
            return Err("GPU compute is not ready; enable it from the GPU button".into());
        }
        c
    } else {
        core.inner().clone()
    };
    tauri::async_runtime::spawn_blocking(move || core.call(&method, params))
        .await
        .map_err(|e| format!("internal error: {e}"))?
}

/// Kill and relaunch the science cores (used by the Stop button).
#[tauri::command]
async fn core_restart(core: State<'_, Arc<CoreHandle>>, compute: State<'_, Arc<Compute>>) -> Result<(), String> {
    let core = core.inner().clone();
    let gpu = compute.inner().core.clone();
    tauri::async_runtime::spawn_blocking(move || {
        if gpu.is_running() {
            let _ = gpu.restart();
        }
        core.restart()
    })
    .await
    .map_err(|e| format!("internal error: {e}"))?
}

#[tauri::command]
fn compute_status(compute: State<'_, Arc<Compute>>) -> ComputeStatus {
    compute.inner().status()
}

#[tauri::command]
fn compute_enable(app: tauri::AppHandle, compute: State<'_, Arc<Compute>>) -> Result<(), String> {
    compute.inner().enable(app)
}

#[tauri::command]
async fn compute_remove(app: tauri::AppHandle, compute: State<'_, Arc<Compute>>) -> Result<(), String> {
    let c = compute.inner().clone();
    tauri::async_runtime::spawn_blocking(move || c.remove(&app))
        .await
        .map_err(|e| format!("internal error: {e}"))?
}

#[derive(Serialize)]
struct UpdateInfo {
    current: String,
    version: String,
    notes: Option<String>,
    date: Option<String>,
}

/// Ask the release feed whether a newer signed version exists.
#[tauri::command]
async fn update_check(app: tauri::AppHandle) -> Result<Option<UpdateInfo>, String> {
    let updater = app.updater().map_err(|e| e.to_string())?;
    let update = updater.check().await.map_err(|e| e.to_string())?;
    Ok(update.map(|u| UpdateInfo {
        current: u.current_version.clone(),
        version: u.version.clone(),
        notes: u.body.clone(),
        date: u.date.map(|d| d.to_string()),
    }))
}

/// Download the update, verify its signature against the key built into this
/// app, run the installer and restart. Progress is reported as events.
#[tauri::command]
async fn update_install(app: tauri::AppHandle, core: State<'_, Arc<CoreHandle>>) -> Result<(), String> {
    let updater = app.updater().map_err(|e| e.to_string())?;
    let Some(update) = updater.check().await.map_err(|e| e.to_string())? else {
        return Err("no update is available".into());
    };
    let progress = app.clone();
    let mut received: u64 = 0;
    let bytes = update
        .download(
            move |chunk, total| {
                received += chunk as u64;
                let _ = progress.emit("update-progress", (received, total));
            },
            || {},
        )
        .await
        .map_err(|e| e.to_string())?;
    // The installer replaces files the science cores may hold open.
    core.inner().shutdown();
    if let Some(c) = app.try_state::<Arc<Compute>>() {
        c.core.shutdown();
    }
    update.install(bytes).map_err(|e| e.to_string())?;
    app.restart();
}

pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_updater::Builder::new().build())
        .setup(|app| {
            let handle = app.handle().clone();
            let core = Arc::new(CoreHandle::new());
            app.manage(core.clone());
            let compute = Arc::new(Compute::new());
            app.manage(compute.clone());
            let chandle = app.handle().clone();
            std::thread::spawn(move || compute.init(&chandle));
            // Prepare the runtime (first launch unpacks the engines) and start the
            // science core off the UI thread; progress is reported as events.
            std::thread::spawn(move || {
                let emit = |s: Status| {
                    let _ = handle.emit("core-status", &s);
                };
                match science::launch_spec(&handle, &emit) {
                    Ok(spec) => {
                        emit(Status::new("starting", "Starting science core"));
                        match core.start(spec) {
                            Ok(()) => emit(Status::new("ready", "Science core ready")),
                            Err(e) => emit(Status::new("error", &e)),
                        }
                    }
                    Err(e) => emit(Status::new("error", &e)),
                }
            });
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![
            core_call, core_restart, update_check, update_install, compute_status, compute_enable, compute_remove
        ])
        .on_window_event(|window, event| {
            if let tauri::WindowEvent::Destroyed = event {
                if let Some(core) = window.try_state::<Arc<CoreHandle>>() {
                    core.shutdown();
                }
                if let Some(c) = window.try_state::<Arc<Compute>>() {
                    c.core.shutdown();
                }
            }
        })
        .run(tauri::generate_context!())
        .expect("error while running Athanor");
}
