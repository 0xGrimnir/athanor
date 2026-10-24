import subprocess, shutil
print(shutil.which("deamSim"))
r = subprocess.run(["deamSim"], capture_output=True, text=True); print((r.stdout + r.stderr)[:600])
r = subprocess.run(["micromamba", "list", "gargammel"], capture_output=True, text=True); print(r.stdout[-400:], r.stderr[-200:])
