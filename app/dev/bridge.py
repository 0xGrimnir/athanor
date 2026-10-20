"""Development bridge: serves app/ui and forwards /rpc to the science core.

Mirrors what the Tauri shell does (spawn `python -m athanor_core.server`, talk
JSON-RPC over pipes) so the interface can be exercised in a normal browser.
Usage:  python app/dev/bridge.py [port]
"""
from __future__ import annotations

import http.server
import itertools
import json
import os
import socketserver
import subprocess
import sys
import threading
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
UI = ROOT / "app" / "ui"


class Core:
    def __init__(self):
        self.lock = threading.Lock()
        self.ids = itertools.count(1)
        self.proc = None
        self.start()

    def start(self):
        env = dict(os.environ, PYTHONPATH=str(ROOT / "core"), PYTHONUNBUFFERED="1")
        self.proc = subprocess.Popen([sys.executable, "-m", "athanor_core.server"], stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE, stderr=sys.stderr, text=True, env=env, bufsize=1)
        self.proc.stdout.readline()  # "ready" notification

    def restart(self):
        try:
            self.proc.kill()
        except Exception:
            pass
        self.start()

    def call(self, method, params):
        with self.lock:
            rid = next(self.ids)
            self.proc.stdin.write(json.dumps({"jsonrpc": "2.0", "id": rid, "method": method, "params": params}) + "\n")
            self.proc.stdin.flush()
            while True:  # skip notifications (e.g. progress) until our response
                line = self.proc.stdout.readline()
                if not line:
                    raise RuntimeError("science core exited")
                msg = json.loads(line)
                if msg.get("id") == rid:
                    return msg


CORE = Core()


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=str(UI), **k)

    def log_message(self, *a):
        pass

    def do_POST(self):
        n = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(n) or b"{}"
        if self.path == "/rpc":
            req = json.loads(body)
            try:
                resp = CORE.call(req["method"], req.get("params") or {})
            except Exception as ex:
                resp = {"error": {"message": str(ex)}}
        elif self.path == "/restart":
            CORE.restart(); resp = {"result": True}
        else:
            self.send_error(404); return
        data = json.dumps(resp).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)


class Server(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True
    allow_reuse_address = True


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8790
    print(f"Athanor dev bridge on http://127.0.0.1:{port}", flush=True)
    Server(("127.0.0.1", port), Handler).serve_forever()
