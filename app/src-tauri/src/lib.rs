mod science;

use std::sync::Arc;

use serde_json::Value;
use tauri::{Emitter, Manager, State};

use crate::science::{CoreHandle, Status};

/// Forward one JSON-RPC call to the science core and return its `result`.
#[tauri::command]
async fn core_call(core: State<'_, Arc<CoreHandle>>, method: String, params: Value) -> Result<Value, String> {
    let core = core.inner().clone();
    tauri::async_runtime::spawn_blocking(move || core.call(&method, params))
        .await
        .map_err(|e| format!("internal error: {e}"))?
}

/// Kill and relaunch the science core (used by the Stop button).
#[tauri::command]
async fn core_restart(core: State<'_, Arc<CoreHandle>>) -> Result<(), String> {
    let core = core.inner().clone();
    tauri::async_runtime::spawn_blocking(move || core.restart())
        .await
        .map_err(|e| format!("internal error: {e}"))?
}

pub fn run() {
    tauri::Builder::default()
        .setup(|app| {
            let handle = app.handle().clone();
            let core = Arc::new(CoreHandle::new());
            app.manage(core.clone());
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
        .invoke_handler(tauri::generate_handler![core_call, core_restart])
        .on_window_event(|window, event| {
            if let tauri::WindowEvent::Destroyed = event {
                if let Some(core) = window.try_state::<Arc<CoreHandle>>() {
                    core.shutdown();
                }
            }
        })
        .run(tauri::generate_context!())
        .expect("error while running Athanor");
}
