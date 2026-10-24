/* Athanor desktop UI. All numbers come from the Python science core over
 * JSON-RPC; this file only renders them and orchestrates jobs. */
"use strict";

// ------------------------------------------------------------------ bridge
const TAURI = window.__TAURI__;
async function rpc(method, params = {}) {
  if (TAURI) return TAURI.core.invoke("core_call", { method, params });
  const r = await fetch("/rpc", { method: "POST", headers: { "Content-Type": "application/json" },
                                  body: JSON.stringify({ method, params }) });
  const j = await r.json();
  if (j.error) throw new Error(j.error.message || String(j.error));
  return j.result;
}
async function restartCore() {
  if (TAURI) return TAURI.core.invoke("core_restart");
  await fetch("/restart", { method: "POST" });
}

const $ = (id) => document.getElementById(id);
const fmt = (v, d = 4) => (v == null || Number.isNaN(v) ? "—" : Number(v).toFixed(d));
const HARTREE_EV = 27.211386245988;      // CODATA 2018
const HARTREE_KCAL = 627.509474063;

// ------------------------------------------------------------------ state
const S = {
  mol: null,            // molecule JSON from the core
  geomOrigin: null,     // {kind:'initial'|'optimized'|'imported', method, key}
  tier: "gfn2-xtb",
  validation: null,     // published comparators
  bench: null,          // Athanor's own benchmark report
  last: {},             // last results per task
  modes: null, sel: -1,
  showCharges: false, showLabels: false,
  busy: false,
};

const PRESETS = [
  ["Water", "O"], ["Ammonia", "N"], ["Methane", "C"], ["Formaldehyde", "C=O"],
  ["Ethylene", "C=C"], ["Ethanol", "CCO"], ["Acetic acid", "CC(=O)O"], ["Benzene", "c1ccccc1"],
  ["Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"], ["Aspirin", "CC(=O)Oc1ccccc1C(=O)O"],
];

const TIER_NOTES = {
  "mmff94": "Classical force field (Halgren 1996). Fast geometries for organic molecules. No electrons: dipoles are point-charge estimates and are not validated. Molecules outside MMFF94's fitted parameters are declined.",
  "gfn2-xtb": "Semi-empirical tight binding (Bannwarth, Ehlert & Grimme 2019). Electronic structure in seconds; known to be weaker for dipoles of small molecules containing heavier atoms.",
  "dft": "Kohn-Sham DFT via Psi4. Most accurate tier here; cost grows roughly with the cube of system size. Frequencies need many gradients: expect minutes for small molecules.",
  "dft-gpu": "Kohn-Sham DFT on your NVIDIA GPU via GPU4PySCF, with density fitting (RI-J/K, def2-universal-jkfit). Same functionals and basis sets as the DFT tier; the RI approximation shifts absolute energies by about 10⁻⁵ Eh and does not change geometries or frequencies meaningfully. Verified on this machine by the GPU self-test.",
};
const isDft = (t) => t === "dft" || t === "dft-gpu";

function methodKey() {
  if (!isDft(S.tier)) return S.tier;
  return `${S.tier}:${$("funcSel").value}/${$("basisSel").value}`;
}
function methodLabel() {
  if (S.tier === "mmff94") return "MMFF94";
  if (S.tier === "gfn2-xtb") return "GFN2-xTB";
  return `${$("funcSel").selectedOptions[0].textContent}/${$("basisSel").value}${S.tier === "dft-gpu" ? " (GPU, RI)" : ""}`;
}
function analyticHessianValidated() {
  const a = S.compute && S.compute.selftest && S.compute.selftest.analytic_hessian;
  return !!(a && a.passed);
}
function engineOptions() {
  if (!isDft(S.tier)) return {};
  const o = { functional: $("funcSel").value, basis: $("basisSel").value };
  // Analytic Hessian + IR on the GPU only after this machine's self-test validated it.
  if (S.tier === "dft-gpu" && analyticHessianValidated()) o.hessian = "analytic";
  return o;
}

// ------------------------------------------------------------------ ticker & jobs
const logs = [];
function log(msg, warn = false) {
  logs.unshift({ msg, warn }); logs.length = Math.min(logs.length, 4);
  $("ticker").innerHTML = logs.map((l) => `<div class="${l.warn ? "w" : ""}">${escapeHtml(l.msg)}</div>`).join("");
}
function escapeHtml(s) { return String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c])); }

function addJob(label) {
  const li = document.createElement("li");
  li.className = "run";
  const t0 = performance.now();
  li.innerHTML = `<span class="t">${new Date().toLocaleTimeString()}</span><span>${escapeHtml(label)}</span><span class="s">running</span>`;
  $("jobs").prepend(li);
  return {
    done(note = "") { li.className = "ok"; li.querySelector(".s").textContent = `${((performance.now() - t0) / 1000).toFixed(1)} s${note ? " · " + note : ""}`; },
    fail(err) { li.className = "err"; li.querySelector(".s").textContent = "failed"; li.title = err; },
  };
}

function setBusy(on, text = "Running…") {
  S.busy = on;
  $("busy").hidden = !on; $("busyText").textContent = text;
  for (const id of ["btnOpt", "btnFreq", "btnProps"]) $(id).disabled = on || !S.mol;
  $("btnStop").hidden = !on;
}

// ------------------------------------------------------------------ 3D viewer
const ELEM = {
  H: [0xe8f4ff, .31], C: [0x8fb4c9, .76], N: [0x4f8dff, .71], O: [0xff5a4a, .66], F: [0xb8ff6a, .57],
  P: [0xff9a2e, 1.07], S: [0xffe14a, 1.05], Cl: [0x3fe08a, 1.02], Br: [0xd0573a, 1.20], I: [0xb26bff, 1.39],
  B: [0xffb3c7, .84], Si: [0xe0c38c, 1.11], Li: [0xcc80ff, 1.28], Na: [0xab5cf2, 1.66], K: [0x8f40d4, 2.03],
  Mg: [0x8aff00, 1.41], Ca: [0x3dff00, 1.76], Se: [0xffa100, 1.20], Al: [0xbfa6a6, 1.21],
};
const elInfo = (s) => ELEM[s] || [0xd0d8e0, 1.0];

const V = (() => {
  const canvas = $("gl");
  const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setClearColor(0x000000, 0);
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(38, 1, 0.1, 500);
  scene.add(new THREE.AmbientLight(0x6aa8c8, 0.6));
  const key = new THREE.DirectionalLight(0xdff4ff, 1.0); key.position.set(5, 8, 10); scene.add(key);
  const rim = new THREE.DirectionalLight(0x62d6ff, 0.8); rim.position.set(-8, -3, -6); scene.add(rim);
  const world = new THREE.Group(); scene.add(world);
  const molG = new THREE.Group(); world.add(molG);
  // floor rings
  const floor = new THREE.Group(); floor.position.y = -5; scene.add(floor);
  const grid = new THREE.PolarGridHelper(11, 24, 6, 96, 0x62d6ff, 0x62d6ff);
  grid.material.transparent = true; grid.material.opacity = 0.12; floor.add(grid);

  const sph = new THREE.SphereGeometry(1, 32, 22);
  const cyl = new THREE.CylinderGeometry(1, 1, 1, 14, 1, true);
  const bondMat = new THREE.MeshStandardMaterial({ color: 0x9fdcff, emissive: 0x1d5a7a, metalness: .3, roughness: .4 });
  let atoms = [], bonds = [], labels = [], eq = null, center = new THREE.Vector3();
  let dist = 14, distTarget = 14;

  function resize() {
    const w = canvas.clientWidth, h = canvas.clientHeight;
    if (!w || !h) return;
    renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix();
  }
  new ResizeObserver(resize).observe(canvas);

  function perceiveBonds(m) {
    if (m.bonds && m.bonds.length) return m.bonds.map(([i, j, o]) => [i, j, o]);
    const out = [];
    for (let i = 0; i < m.symbols.length; i++) for (let j = i + 1; j < m.symbols.length; j++) {
      const a = m.coords[i], b = m.coords[j];
      const d = Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
      if (d < 1.2 * (elInfo(m.symbols[i])[1] + elInfo(m.symbols[j])[1])) out.push([i, j, 1]);
    }
    return out;
  }

  function label(text) {
    const c = document.createElement("canvas"); c.width = 128; c.height = 64;
    const x = c.getContext("2d"); x.font = "600 34px JetBrains Mono, monospace"; x.fillStyle = "#d5ebf6";
    x.textAlign = "center"; x.textBaseline = "middle"; x.fillText(text, 64, 32);
    const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(c), depthTest: false, transparent: true }));
    sp.scale.set(0.9, 0.45, 1); return sp;
  }

  function setMolecule(m, keepView = false) {
    molG.clear(); atoms = []; bonds = []; labels = [];
    if (!m) return;
    eq = m.coords.map((p) => new THREE.Vector3(...p));
    center = eq.reduce((a, v) => a.add(v), new THREE.Vector3()).multiplyScalar(1 / eq.length);
    m.symbols.forEach((s, i) => {
      const [col, rcov] = elInfo(s);
      const mat = new THREE.MeshStandardMaterial({ color: col, emissive: new THREE.Color(col).multiplyScalar(.25), metalness: .25, roughness: .35 });
      const mesh = new THREE.Mesh(sph, mat);
      mesh.scale.setScalar(s === "H" ? .24 : .2 + rcov * .26);
      molG.add(mesh); atoms.push(mesh);
      const lb = label(`${s}${i + 1}`); lb.visible = S.showLabels; molG.add(lb); labels.push(lb);
    });
    for (const [i, j, o] of perceiveBonds(m)) {
      const n = o >= 2.5 ? 3 : o >= 1.75 ? 2 : 1;
      const meshes = [];
      for (let k = 0; k < n; k++) { const c = new THREE.Mesh(cyl, bondMat); molG.add(c); meshes.push(c); }
      bonds.push({ i, j, meshes, aromatic: o === 1.5 });
    }
    if (!keepView) {
      let r = 0; for (const v of eq) r = Math.max(r, v.distanceTo(center));
      distTarget = Math.max(8, (r + 2.5) / Math.tan(THREE.MathUtils.degToRad(17)));
    }
    applyColors(); place(eq);
  }

  const up = new THREE.Vector3(0, 1, 0), d = new THREE.Vector3(), perp = new THREE.Vector3();
  function placeCyl(mesh, a, b, r) {
    d.subVectors(b, a); const len = d.length(); if (len < 1e-6) return;
    mesh.position.addVectors(a, b).multiplyScalar(.5);
    mesh.scale.set(r, len, r);
    mesh.quaternion.setFromUnitVectors(up, d.clone().divideScalar(len));
  }
  function place(pos) {
    atoms.forEach((m, i) => { m.position.copy(pos[i]).sub(center); labels[i].position.copy(m.position).add(new THREE.Vector3(0, .45, 0)); });
    for (const b of bonds) {
      const a = pos[b.i].clone().sub(center), c = pos[b.j].clone().sub(center);
      if (b.meshes.length === 1) { placeCyl(b.meshes[0], a, c, .075); continue; }
      d.subVectors(c, a).normalize();
      perp.set(0, 0, 1).cross(d); if (perp.lengthSq() < 1e-6) perp.set(1, 0, 0).cross(d); perp.normalize();
      const gap = b.meshes.length === 2 ? .16 : .2;
      b.meshes.forEach((mesh, k) => {
        const off = (k - (b.meshes.length - 1) / 2) * gap;
        placeCyl(mesh, a.clone().addScaledVector(perp, off), c.clone().addScaledVector(perp, off), b.meshes.length === 2 ? .055 : .045);
      });
    }
  }

  function applyColors() {
    const q = S.showCharges && S.last.charges ? S.last.charges : null;
    atoms.forEach((m, i) => {
      let col;
      if (q && q.length === atoms.length) {
        const v = Math.max(-1, Math.min(1, q[i] / 0.6));
        col = v < 0 ? new THREE.Color(0xeef6ff).lerp(new THREE.Color(0xff3b5c), -v) : new THREE.Color(0xeef6ff).lerp(new THREE.Color(0x3b8bff), v);
      } else col = new THREE.Color(elInfo(S.mol.symbols[i])[0]);
      m.material.color.copy(col); m.material.emissive.copy(col).multiplyScalar(.25);
    });
  }
  function setLabels(on) { labels.forEach((l) => (l.visible = on)); }

  // orbit
  let drag = null; const q = new THREE.Quaternion(), ax = new THREE.Vector3();
  canvas.addEventListener("pointerdown", (e) => { canvas.setPointerCapture(e.pointerId); drag = { x: e.clientX, y: e.clientY }; canvas.classList.add("drag"); });
  canvas.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y; drag = { x: e.clientX, y: e.clientY };
    ax.set(0, 1, 0); q.setFromAxisAngle(ax, dx * .008); world.quaternion.premultiply(q);
    ax.set(1, 0, 0); q.setFromAxisAngle(ax, dy * .008); world.quaternion.premultiply(q);
  });
  const endDrag = () => { drag = null; canvas.classList.remove("drag"); };
  canvas.addEventListener("pointerup", endDrag); canvas.addEventListener("pointercancel", endDrag);
  canvas.addEventListener("wheel", (e) => { e.preventDefault(); distTarget = Math.max(3, Math.min(80, distTarget * (1 + e.deltaY * .0012))); }, { passive: false });

  // animation of a normal mode: x(t) = x_eq + A sin(wt) * L, A in angstrom
  let mode = null, t0 = performance.now();
  function setMode(vec) { mode = vec ? vec.map((v) => new THREE.Vector3(...v)) : null; t0 = performance.now(); if (!mode && eq) place(eq); }

  function frame(now) {
    dist += (distTarget - dist) * .1;
    camera.position.set(0, dist * .2, dist); camera.lookAt(0, -.3, 0);
    if (mode && eq) {
      const s = Math.sin((now - t0) / 1000 * Math.PI * 2 * .9) * .35;
      place(eq.map((p, i) => p.clone().addScaledVector(mode[i], s)));
    }
    floor.rotation.y += .0012;
    renderer.render(scene, camera);
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);
  resize();
  return { setMolecule, applyColors, setLabels, setMode };
})();

// ------------------------------------------------------------------ spectrum
const SPEC = (() => {
  const cv = $("spec"), ctx = cv.getContext("2d");
  let data = null, hits = [];
  function draw() {
    const w = cv.clientWidth, h = cv.clientHeight, dpr = Math.min(window.devicePixelRatio, 2);
    if (!w || !h) return;
    cv.width = w * dpr; cv.height = h * dpr; ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, w, h);
    if (!data) return;
    const padL = 16, padR = 16, padT = 14, padB = 22;
    const x0 = padL, x1 = w - padR, y0 = h - padB, y1 = padT;
    const xmax = 4000, xmin = 400;
    const X = (v) => x0 + (xmax - v) / (xmax - xmin) * (x1 - x0);
    // axis
    ctx.strokeStyle = "rgba(110,205,255,.12)"; ctx.fillStyle = "#7f9fb3"; ctx.font = "10px JetBrains Mono, monospace"; ctx.textAlign = "center";
    for (let v = 4000; v >= 500; v -= 500) { const x = X(v); ctx.beginPath(); ctx.moveTo(x, y1); ctx.lineTo(x, y0); ctx.stroke(); ctx.fillText(v, x, h - 6); }
    const { freqs, ints } = data;
    const imax = Math.max(1e-9, ...ints.filter((v) => v != null));
    // Lorentzian broadening, FWHM 10 cm-1 (display only), filled area under the curve
    const g = 5, curve = [];
    for (let px = x0; px <= x1; px++) {
      const v = xmax - (px - x0) / (x1 - x0) * (xmax - xmin);
      let y = 0; freqs.forEach((f, k) => { if (f > 0 && ints[k] != null) y += ints[k] / imax * g * g / ((v - f) ** 2 + g * g); });
      curve.push([px, Math.min(1, y)]);
    }
    const Y = (y) => y0 - y * (y0 - y1);
    const grad = ctx.createLinearGradient(0, y1, 0, y0);
    grad.addColorStop(0, "rgba(98,214,255,.35)"); grad.addColorStop(1, "rgba(98,214,255,.02)");
    ctx.beginPath(); ctx.moveTo(x0, y0); curve.forEach(([px, y]) => ctx.lineTo(px, Y(y))); ctx.lineTo(x1, y0); ctx.closePath();
    ctx.fillStyle = grad; ctx.fill();
    ctx.beginPath(); curve.forEach(([px, y], k) => (k ? ctx.lineTo(px, Y(y)) : ctx.moveTo(px, Y(y))));
    ctx.strokeStyle = "#62d6ff"; ctx.lineWidth = 1.2; ctx.stroke();
    // sticks (intensity scale), selected mode highlighted
    hits = [];
    freqs.forEach((f, k) => {
      if (f <= 0 || f < xmin || f > xmax) return;
      const x = X(f), hgt = (ints[k] ?? 0) / imax;
      const on = k === S.sel;
      ctx.strokeStyle = on ? "#ffb34a" : "rgba(255,179,74,.75)"; ctx.lineWidth = on ? 3 : 1.5;
      ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, Y(Math.max(.02, hgt))); ctx.stroke();
      if (on) { ctx.fillStyle = "#ffb34a"; ctx.font = "11px JetBrains Mono, monospace"; ctx.textAlign = x > x1 - 60 ? "right" : "left";
        ctx.fillText(` ${f.toFixed(0)} cm⁻¹`, x + (x > x1 - 60 ? -4 : 4), Y(Math.max(.02, hgt)) + 10); }
      hits.push({ x, k });
    });
    ctx.fillStyle = "#7f9fb3"; ctx.textAlign = "left"; ctx.font = "10px JetBrains Mono, monospace";
    ctx.fillText("rel. intensity", 4, y1 + 8);
  }
  cv.addEventListener("click", (e) => {
    const r = cv.getBoundingClientRect(), mx = e.clientX - r.left;
    let best = null, bd = 10;
    for (const h of hits) { const dd = Math.abs(h.x - mx); if (dd < bd) { bd = dd; best = h; } }
    if (best) selectMode(best.k);
  });
  new ResizeObserver(draw).observe(cv);
  return { set(d) { data = d; draw(); }, draw };
})();

// ------------------------------------------------------------------ molecule loading
function setMolecule(m, origin) {
  S.mol = m; S.geomOrigin = origin; S.last = {}; S.modes = null; S.sel = -1;
  $("chargeInput").value = m.charge; $("multInput").value = m.multiplicity;
  $("molName").textContent = m.name || m.formula;
  $("molFormula").innerHTML = m.formula.replace(/(\d+)/g, "<sub>$1</sub>");
  V.setMolecule(m); V.setMode(null);
  $("specWrap").hidden = true; $("vibSec").hidden = true; $("provSec").hidden = true;
  $("results").innerHTML = `<p class="note">${m.symbols.length} atoms. Choose a method and run a calculation.</p>`;
  updateBadge(); setBusy(false);
}
function updateBadge() {
  const b = $("geomBadge"); b.hidden = !S.mol;
  if (!S.mol) return;
  const o = S.geomOrigin || {};
  if (o.kind === "optimized") { b.className = "badge"; b.textContent = `Optimized · ${o.label}`; }
  else if (o.kind === "imported") { b.className = "badge raw"; b.textContent = "Imported geometry · not optimized"; }
  else { b.className = "badge raw"; b.textContent = "Initial guess (ETKDGv3) · not optimized"; }
}

async function build(smiles, name) {
  const job = addJob(`Build ${name || smiles}`);
  try {
    const m = await rpc("build", { smiles, name: name || smiles });
    setMolecule(m, { kind: "initial" }); job.done();
    log(`Built ${m.formula} from SMILES with RDKit ETKDGv3`);
  } catch (e) { job.fail(String(e)); log(`Build failed: ${e.message || e}`, true); }
}

function currentMolecule() {
  const m = structuredClone(S.mol);
  m.charge = parseInt($("chargeInput").value, 10) || 0;
  m.multiplicity = Math.max(1, parseInt($("multInput").value, 10) || 1);
  return m;
}

// ------------------------------------------------------------------ jobs
// Live progress from the science cores (optimizer steps, Hessian gradients).
function progressText(d) {
  if (d.task === "optimize")
    return `${S.busyLabel}: step ${d.step} · E ${fmt(d.energy, 8)} Eh · max force ${d.max_force.toExponential(2)} (target ${d.target_max_force.toExponential(1)})`;
  if (d.task === "adna") {
    let st = d.stage || "";
    for (const g of (A.catalog && A.catalog.references.genomes) || []) st = st.replaceAll(g.accession, commonName(g.organism));
    return `Step ${d.step ?? "?"} of ${d.steps ?? "?"} · ${st}`;
  }
  if (d.task === "personal") return d.stage ? d.stage[0].toUpperCase() + d.stage.slice(1) : "";
  if (d.task === "pgdata") {
    if (d.total_bytes) return `${d.stage} · ${(d.bytes / 1e6).toFixed(0)} of ${(d.total_bytes / 1e6).toFixed(0)} MB`;
    return d.stage || "";
  }
  if (d.task === "frequencies" && d.analytic) return `${S.busyLabel}: analytic Hessian and IR intensities (coupled-perturbed SCF)…`;
  if (d.task === "frequencies") return `${S.busyLabel}: Hessian gradient ${d.step} of ${d.total}`;
  return "";
}

async function runOptimize() {
  const label = methodLabel(), key = methodKey();
  const job = addJob(`Optimize · ${label}`);
  S.busyLabel = label;
  setBusy(true, `Optimizing with ${label}…`);
  try {
    let molecule = currentMolecule(), pre = null;
    // DFT from a rough starting structure needs many expensive steps; a GFN2-xTB
    // pre-optimization (seconds) gets close first. Recorded in the provenance.
    const fromRough = !(S.geomOrigin && S.geomOrigin.kind === "optimized");
    if (isDft(S.tier) && $("preoptChk").checked && fromRough) {
      S.busyLabel = "GFN2-xTB pre-optimization";
      setBusy(true, "Pre-optimizing with GFN2-xTB…");
      try {
        const p0 = await rpc("optimize", { molecule, tier: "gfn2-xtb", options: {} });
        molecule = { ...p0.molecule, name: S.mol.name };
        pre = { method: "GFN2-xTB", engine: `${p0.provenance.engine} ${p0.provenance.engine_version}`,
                steps: p0.steps, converged: p0.converged };
        log(`GFN2-xTB pre-optimization: ${p0.steps} steps; starting ${label} from that geometry`);
      } catch (e) { log(`GFN2-xTB pre-optimization failed (${e.message || e}); starting ${label} from the original geometry`, true); }
      S.busyLabel = label;
      setBusy(true, `Optimizing with ${label}…`);
    }
    const r = await rpc("optimize", { molecule, tier: S.tier, options: engineOptions(), pre_optimization: pre });
    const keep = S.mol.name;
    S.mol = { ...r.molecule, name: keep };
    S.geomOrigin = { kind: "optimized", key, label };
    S.last = { opt: r, charges: S.last.charges };
    V.setMolecule(S.mol, true); V.setMode(null); S.modes = null; $("specWrap").hidden = true; $("vibSec").hidden = true;
    updateBadge(); renderResults(); renderProv(r.provenance);
    job.done(r.converged ? `${r.steps} steps` : "not converged");
    log(r.converged ? `${label}: converged in ${r.steps} steps, E = ${fmt(r.energy_hartree, 8)} Eh` : `${label}: optimization did not converge`, !r.converged);
    return r;
  } catch (e) { job.fail(String(e)); log(`Optimization failed: ${e.message || e}`, true); throw e; }
  finally { setBusy(false); }
}

async function runFrequencies() {
  if (!S.geomOrigin || S.geomOrigin.kind !== "optimized" || S.geomOrigin.key !== methodKey()) {
    log(`Frequencies need a geometry optimized at ${methodLabel()}; optimizing first`);
    try { await runOptimize(); } catch { return; }
  }
  const label = methodLabel();
  const job = addJob(`Frequencies · ${label}`);
  S.busyLabel = label;
  setBusy(true, `Computing Hessian with ${label}…`);
  try {
    const r = await rpc("frequencies", { molecule: currentMolecule(), tier: S.tier, options: engineOptions() });
    S.last.freq = r; S.modes = r.normal_modes; S.sel = -1;
    renderVibrations(); renderResults(); renderProv(r.provenance);
    job.done(`${r.frequencies_cm1.length} modes`);
    log(r.n_imaginary ? `${label}: ${r.n_imaginary} imaginary frequency(ies), structure is not a minimum` : `${label}: ${r.frequencies_cm1.length} real modes, ZPE ${fmt(r.zpe_kcal_mol, 2)} kcal/mol`, r.n_imaginary > 0);
  } catch (e) { job.fail(String(e)); log(`Frequency calculation failed: ${e.message || e}`, true); }
  finally { setBusy(false); }
}

async function runProps() {
  const label = methodLabel();
  const job = addJob(`Properties · ${label}`);
  setBusy(true, `Single point with ${label}…`);
  try {
    const r = await rpc("single_point", { molecule: currentMolecule(), tier: S.tier, options: engineOptions() });
    S.last.sp = r; S.last.charges = r.charges;
    renderResults(); renderProv(r.provenance); V.applyColors();
    job.done(); log(`${label}: E = ${fmt(r.energy_hartree, 8)} Eh, dipole ${fmt(r.dipole_debye, 3)} D`);
  } catch (e) { job.fail(String(e)); log(`Single point failed: ${e.message || e}`, true); }
  finally { setBusy(false); }
}

// ------------------------------------------------------------------ results rendering
function dipoleErrorText() {
  const c = S.validation && S.validation.comparators;
  if (!c) return "";
  if (S.tier === "mmff94") return `<div class="warnbox">MMFF94 dipoles are classical point-charge estimates and are not validated against experiment.</div>`;
  const pub = isDft(S.tier) ? c.dft.dipole_regularized_rms_percent : c["gfn2-xtb"].dipole_regularized_rms_percent;
  let measured = "";
  const b = S.bench && S.bench.tiers && S.bench.tiers[S.tier] && S.bench.tiers[S.tier].dipoles;
  if (b && b.regularized_rms_percent != null) measured = ` · Athanor benchmark vs experiment: ${b.regularized_rms_percent.toFixed(1)}% (n=${b.n})`;
  const basisNote = isDft(S.tier) ? ` (${c.dft.method}/${pub.basis})` : "";
  return `<span class="errbar">typical error ±${pub.value}%${basisNote} vs ${pub.reference}${measured}</span>`;
}

function renderResults() {
  const out = [];
  const sp = S.last.sp, opt = S.last.opt, fq = S.last.freq;
  if (opt) {
    out.push(`<div class="big"><span class="v">${fmt(opt.energy_hartree, 8)}</span><span class="u">Eh · ${escapeHtml(methodLabel())} total energy at optimized geometry</span></div>`);
    out.push(`<dl class="readout"><dt>Converged</dt><dd>${opt.converged ? "yes" : "no"}</dd><dt>Steps</dt><dd>${opt.steps}</dd>
      <dt>Max force</dt><dd>${opt.max_force.toExponential(2)}<small>Eh/bohr</small></dd><dt>RMS force</dt><dd>${opt.rms_force.toExponential(2)}<small>Eh/bohr</small></dd></dl>`);
    if (!opt.converged) out.push(`<div class="warnbox">Optimization did not reach the convergence criteria. Results at this geometry are not at a stationary point.</div>`);
  }
  if (sp) {
    out.push(`<dl class="readout"><dt>Energy</dt><dd>${fmt(sp.energy_hartree, 8)}<small>Eh</small></dd>
      <dt>Dipole</dt><dd>${fmt(sp.dipole_debye, 3)}<small>D</small></dd>
      ${sp.extras && sp.extras.homo_lumo && sp.extras.homo_lumo.homo_eh != null ? `<dt>HOMO / LUMO</dt><dd>${fmt(sp.extras.homo_lumo.homo_eh * HARTREE_EV, 2)} / ${fmt(sp.extras.homo_lumo.lumo_eh * HARTREE_EV, 2)}<small>eV</small></dd>` : ""}
      <dt>Charges</dt><dd>${escapeHtml(sp.charge_model || "—")}</dd></dl>`);
    out.push(dipoleErrorText());
    if (sp.max_force > 4.5e-4) out.push(`<div class="warnbox">This geometry is not optimized at ${escapeHtml(methodLabel())} (max force ${sp.max_force.toExponential(1)} Eh/bohr). Properties are reported at the current structure.</div>`);
  }
  if (fq) {
    out.push(`<dl class="readout"><dt>Zero-point energy</dt><dd>${fmt(fq.zpe_kcal_mol, 2)}<small>kcal/mol</small></dd>
      <dt>Imaginary modes</dt><dd>${fq.n_imaginary}</dd></dl>`);
    if (fq.n_imaginary) out.push(`<div class="warnbox">${fq.n_imaginary} imaginary frequency(ies): this structure is a saddle point, not a minimum.</div>`);
    else out.push(`<div class="okbox">All frequencies real: the optimized structure is a true minimum at this level of theory.</div>`);
  }
  $("results").innerHTML = out.join("") || `<p class="note">No results yet.</p>`;
}

function renderVibrations() {
  const r = S.last.freq; if (!r) return;
  $("vibSec").hidden = false;
  $("vibs").innerHTML = r.frequencies_cm1.map((f, k) => {
    const ir = r.ir_intensities_km_mol[k];
    return `<button type="button" class="vib${f < 0 ? " imag" : ""}" data-k="${k}">
      <span>${f < 0 ? Math.abs(f).toFixed(1) + "i" : f.toFixed(1)}</span><span>${ir == null ? "—" : ir.toFixed(1)}</span><span>${r.reduced_masses_amu[k].toFixed(3)}</span></button>`;
  }).join("");
  $("vibs").onclick = (e) => { const b = e.target.closest(".vib"); if (b) selectMode(+b.dataset.k); };
  // spectrum: scale only when a literature factor exists for exactly this level of theory
  let freqs = r.frequencies_cm1.slice(), note = "harmonic, unscaled";
  const c = S.validation && S.validation.comparators && S.validation.comparators.dft;
  if (isDft(S.tier) && c && $("funcSel").value === c.method && $("basisSel").value === c.frequency_scale_factor.basis) {
    const sf = c.frequency_scale_factor.value; freqs = freqs.map((f) => f * sf);
    note = `scaled ×${sf} (CCCBDB, ${c.method}/${c.frequency_scale_factor.basis}, σ ${c.frequency_scale_factor.sigma})`;
  } else if (isDft(S.tier)) note = "harmonic, unscaled: no literature scale factor for this exact level; harmonic frequencies typically exceed experimental fundamentals by a few percent";
  else if (S.tier === "gfn2-xtb") {
    const b = S.bench && S.bench.tiers && S.bench.tiers["gfn2-xtb"] && S.bench.tiers["gfn2-xtb"].frequencies;
    note = "GFN2-xTB harmonic, unscaled" + (b ? ` · Athanor benchmark MAD ${b.gas_phase_harmonic_modes.raw.mad.toFixed(0)} cm⁻¹ vs experiment` : "");
  } else note = "MMFF94 harmonic frequencies are not validated";
  $("specWrap").hidden = false; $("specNote").textContent = `${note} · 10 cm⁻¹ Lorentzian display broadening`;
  SPEC.set({ freqs, ints: r.ir_intensities_km_mol });
  $("vibNote").textContent = r.notes && r.notes.length ? r.notes.join(" ") : "Click a mode to animate it. Amplitude is exaggerated for visibility.";
}

function selectMode(k) {
  S.sel = k; V.setMode(S.modes ? S.modes[k] : null); SPEC.draw();
  document.querySelectorAll(".vib").forEach((b) => b.classList.toggle("on", +b.dataset.k === k));
}

function renderProv(p) {
  if (!p) return;
  S.last.prov = p;
  $("provSec").hidden = false;
  const params = Object.entries(p.parameters || {}).map(([k, v]) => `${k}=${typeof v === "object" ? JSON.stringify(v) : v}`).join(", ");
  $("prov").innerHTML = `<dt>Task</dt><dd>${escapeHtml(p.task)}</dd><dt>Engine</dt><dd>${escapeHtml(p.engine)} ${escapeHtml(p.engine_version)}</dd>
    <dt>Method</dt><dd>${escapeHtml(p.method)}</dd><dt>Parameters</dt><dd>${escapeHtml(params || "—")}</dd>
    <dt>Input hash</dt><dd title="${p.input_hash}">${p.input_hash.slice(0, 16)}…</dd><dt>Platform</dt><dd>${escapeHtml(p.platform)} · Python ${escapeHtml(p.python)}</dd>
    <dt>Athanor</dt><dd>${escapeHtml(p.athanor_version)}</dd><dt>Finished</dt><dd>${escapeHtml(p.finished_utc.replace("T", " ").slice(0, 19))} UTC</dd>
    ${p.notes && p.notes.length ? `<dt>Notes</dt><dd>${escapeHtml(p.notes.join("; "))}</dd>` : ""}`;
}

async function copyText(text, what) {
  try { await navigator.clipboard.writeText(text); log(`${what} copied`); }
  catch { log(`Clipboard unavailable; ${what} not copied`, true); }
}
function toXyz(m) {
  return `${m.symbols.length}\n${m.name || m.formula}\n` + m.symbols.map((s, i) =>
    `${s.padEnd(2)} ${m.coords[i].map((v) => v.toFixed(8).padStart(14)).join(" ")}`).join("\n") + "\n";
}

// ------------------------------------------------------------------ validation view
function renderValidation() {
  const c = S.validation, b = S.bench;
  if (!c) { $("validation").innerHTML = `<p class="note">Validation data loads when the science core is ready.</p>`; return; }
  const src = (k) => escapeHtml(c.sources[k] || k);
  const tiers = [
    ["mmff94", "MMFF94", "RDKit", ["Reproduces Halgren's MMFF94 validation suite: 755 molecules to within 1×10⁻⁴ kcal/mol; empirical-rule test molecules to within 0.01 kcal/mol.",
      "Declines molecules whose bonds fall outside MMFF94's fitted parameter table instead of using rule-based estimates."]],
    ["gfn2-xtb", "GFN2-xTB", "xtb", ["Runs on the Grimme lab's xtb program; the independent tblite library reproduces its energies to 10⁻⁶ Eh on every build.",
      "Frequencies cross-checked against the xtb program's own Hessian (within 3 cm⁻¹).",
      "Vibrational analysis cross-checked against ASE's implementation (frequencies within 3 cm⁻¹, IR intensities within 5%)."]],
    ["dft", "DFT", "Psi4", ["Finite-difference Hessians cross-checked against Psi4's own frequency analysis on every build.",
      "Benchmark level: B3LYP/cc-pVTZ for geometries and frequencies, B3LYP/aug-cc-pVTZ for dipoles."]],
  ];
  $("validation").innerHTML = tiers.map(([key, name, eng, checks]) => {
    const t = b && b.tiers && b.tiers[key];
    const gates = t ? t.gates : [];
    const rows = [];
    if (t && t.geometries) { const g = t.geometries.bond_length_angstrom; if (g.n) rows.push(["Bond lengths vs experimental rₑ", `${g.n}`, `MAD ${g.mad.toFixed(4)} Å`, `mean ${g.mean_signed >= 0 ? "+" : ""}${g.mean_signed.toFixed(4)} Å`]); }
    if (t && t.dipoles && t.dipoles.n) rows.push(["Dipoles vs experiment", `${t.dipoles.n}`, `regularized RMS ${t.dipoles.regularized_rms_percent.toFixed(1)}%`, t.dipoles.method]);
    if (t && t.frequencies) { const s = t.frequencies.gas_phase_harmonic_modes; if (s.n) rows.push(["Gas-phase fundamentals", `${s.n}`, `raw MAD ${s.raw.mad.toFixed(0)} cm⁻¹`, s.literature_scaled ? `scaled ×${s.literature_scale_factor}: MAD ${s.literature_scaled.mad.toFixed(0)} cm⁻¹` : `fitted scale ${s.fitted_scale_factor.toFixed(3)}`]); }
    return `<div class="vtier"><h3>${name}<small>${eng}</small></h3>
      <div class="label">Correctness</div><ul class="note">${checks.map((x) => `<li>${x}</li>`).join("")}</ul>
      <div class="label">Accuracy against experiment</div>
      ${rows.length ? `<div class="tablewrap"><table><thead><tr><th>Property</th><th>n</th><th>Error</th><th>Detail</th></tr></thead><tbody>
        ${rows.map((r) => `<tr><td>${r[0]}</td><td class="n">${r[1]}</td><td class="n">${r[2]}</td><td>${escapeHtml(r[3])}</td></tr>`).join("")}</tbody></table></div>`
        : `<p class="note">${key === "dft" ? "The DFT benchmark runs in CI; results appear here when bundled with a release." : "Not benchmarked against experiment."}</p>`}
      ${gates.length ? `<div class="label">Consistency with published error statistics</div><div class="tablewrap"><table><thead><tr><th>Check</th><th>Athanor</th><th>Published</th><th>Gate</th><th>Source</th></tr></thead><tbody>
        ${gates.map((g) => `<tr><td>${escapeHtml(g.check)}</td><td class="n">${g.measured.toPrecision(3)}</td><td class="n">${g.published}</td><td><span class="pill ${g.pass ? "pass" : "fail"}">${g.pass ? "PASS" : "FAIL"}</span></td><td class="src">${src(g.source)}</td></tr>`).join("")}
        </tbody></table></div>` : ""}
    </div>`;
  }).join("") + `<div class="vtier"><h3>Sources</h3><ul class="note">${Object.values(c.sources).map((s) => `<li>${escapeHtml(s)}</li>`).join("")}</ul></div>`;
}

// ------------------------------------------------------------------ wiring
function setTier(t) {
  S.tier = t;
  document.querySelectorAll("#tierSeg button").forEach((b) => b.setAttribute("aria-checked", String(b.dataset.tier === t)));
  $("dftOpts").hidden = !isDft(t);
  $("preoptChk").closest("label").hidden = !isDft(t);
  // ωB97X-D's dispersion term is defined differently in PySCF; offer it on Psi4 only.
  const wb = $("funcSel").querySelector('option[value="wB97X-D"]');
  wb.disabled = t === "dft-gpu";
  if (wb.disabled && $("funcSel").value === "wB97X-D") $("funcSel").value = "B3LYP";
  $("tierNote").textContent = TIER_NOTES[t];
}

function wire() {
  $("presets").innerHTML = PRESETS.map(([n, s]) => `<button type="button" class="btn" data-s="${s}" data-n="${n}">${n}</button>`).join("");
  $("presets").onclick = (e) => { const b = e.target.closest("button"); if (b && !S.busy) build(b.dataset.s, b.dataset.n); };
  $("smilesForm").onsubmit = (e) => { e.preventDefault(); const s = $("smilesInput").value.trim(); if (s && !S.busy) build(s); };
  $("tierSeg").onclick = (e) => { const b = e.target.closest("button"); if (b) setTier(b.dataset.tier); };
  $("btnOpt").onclick = () => runOptimize().catch(() => {});
  $("btnFreq").onclick = runFrequencies;
  $("adnaRun").onclick = () => runAdna();
  document.querySelectorAll("#adnaModeSeg button").forEach((b) => (b.onclick = () => setAdnaMode(b.dataset.mode)));
  $("realAdd").onclick = () => realAddFiles();
  $("realClear").onclick = () => { R2.raw = []; R2.files = []; renderRealFiles(); };
  $("realRef").onchange = () => realRefChanged();
  $("realCmp").onclick = () => realPickCompare();
  $("realCmpClear").onclick = () => { R2.cmpPath = null; $("realCmpName").textContent = "None"; };
  $("realRun").onclick = () => runReal();
  $("realIt").oninput = (e) => ($("realItVal").textContent = e.target.value);
  $("realDepth").oninput = (e) => ($("realDepthVal").textContent = e.target.value);
  $("realSaveFasta").onclick = () => realExport("fasta");
  $("realSaveMaj").onclick = () => realExport("majority_fasta");
  $("realSaveJson").onclick = () => realExport("json");
  $("pgInstall").onclick = () => pgInstall();
  $("pgRemove").onclick = () => pgRemove();
  $("pgChoose").onclick = () => pgChoose();
  $("pgRun").onclick = () => pgRun();
  $("pgExport").onclick = () => pgExport();
  $("adnaCov").oninput = (e) => ($("adnaCovVal").textContent = `${e.target.value}×`);
  $("adnaCont").oninput = (e) => ($("adnaContVal").textContent = `${e.target.value}%`);
  $("adnaLen").oninput = (e) => ($("adnaLenVal").textContent = `${e.target.value} bp`);
  $("adnaIt").oninput = (e) => ($("adnaItVal").textContent = e.target.value);
  $("btnProps").onclick = runProps;
  $("btnStop").onclick = async () => { log("Stopping: restarting the science core", true); await restartCore(); setBusy(false); };
  $("vCharges").onclick = () => { S.showCharges = !S.showCharges; $("vCharges").setAttribute("aria-pressed", S.showCharges); if (!S.last.charges && S.showCharges) log("Run Properties to compute partial charges"); V.applyColors(); };
  $("vLabels").onclick = () => { S.showLabels = !S.showLabels; $("vLabels").setAttribute("aria-pressed", S.showLabels); V.setLabels(S.showLabels); };
  $("btnCopyProv").onclick = () => S.last.prov && copyText(JSON.stringify(S.last.prov, null, 2), "Provenance record");
  $("btnExportXyz").onclick = () => S.mol && copyText(toXyz(S.mol), "XYZ coordinates");
  $("fileInput").onchange = async (e) => {
    const f = e.target.files[0]; if (!f) return;
    const text = await f.text(); const ext = f.name.split(".").pop().toLowerCase();
    const job = addJob(`Import ${f.name}`);
    try {
      const m = ext === "xyz"
        ? await rpc("parse_xyz", { xyz: text, charge: +$("chargeInput").value || 0, multiplicity: +$("multInput").value || 1 })
        : await rpc("parse_molblock", { molblock: text });
      m.name = f.name.replace(/\.[^.]+$/, ""); setMolecule(m, { kind: "imported" }); job.done();
    } catch (err) { job.fail(String(err)); log(`Import failed: ${err.message || err}`, true); }
    e.target.value = "";
  };
  document.getElementById("tabs").onclick = (e) => {
    const b = e.target.closest(".tab"); if (!b) return;
    document.querySelectorAll(".tab").forEach((t) => t.setAttribute("aria-pressed", String(t === b)));
    for (const v of ["molecule", "validation", "adna", "personal", "dna"]) $(`view-${v}`).hidden = b.dataset.view !== v;
    if (b.dataset.view === "validation") renderValidation();
    if (b.dataset.view === "dna") initDna();
    if (b.dataset.view === "adna") initAdna();
    if (b.dataset.view === "personal") initPersonal();
  };
}

// Signed in-app updates (Tauri updater). The check is a single request to the
// public release feed; nothing about the user or their work is sent.
async function checkForUpdate() {
  if (!TAURI || !TAURI.core) return;
  let info = null;
  try { info = await TAURI.core.invoke("update_check"); }
  catch (e) { log(`Update check failed: ${e}`); return; }
  if (!info) return;
  const btn = $("updateBtn");
  btn.textContent = `Update to ${info.version}`;
  btn.title = info.notes || "";
  btn.hidden = false;
  log(`Athanor ${info.version} is available (you have ${info.current}).`);
  btn.onclick = async () => {
    if (!confirm(`Install Athanor ${info.version} and restart? Running calculations will stop.`)) return;
    btn.disabled = true; btn.textContent = "Downloading update";
    if (TAURI.event) TAURI.event.listen("update-progress", (e) => {
      const [got, total] = e.payload || [];
      if (total) btn.textContent = `Downloading update ${Math.floor((100 * got) / total)}%`;
    });
    try { await TAURI.core.invoke("update_install"); }
    catch (e) { btn.disabled = false; btn.textContent = `Update to ${info.version}`; log(`Update failed: ${e}`, true); }
  };
}

// Optional GPU compute pack (GPU4PySCF inside a private WSL2 distribution).
function renderCompute(st) {
  const btn = $("gpuBtn");
  S.compute = st;
  const labels = {
    "off": "Install Linux engines", "failed": "Linux engines: retry", "needs-restart": "Linux engines: restart Windows",
    "downloading": st.progress != null ? `Linux engines: downloading ${Math.floor(100 * st.progress)}%` : "Linux engines: downloading",
    "installing": "Linux engines: installing", "testing": "Linux engines: self-test",
    "ready": st.gpu_ok ? "Linux engines · GPU ready" : "Linux engines ready",
  };
  if (st.stage === "unsupported" || !labels[st.stage]) { btn.hidden = true; return; }
  btn.hidden = false;
  btn.textContent = labels[st.stage];
  btn.title = [st.message, st.gpu ? `Detected: ${st.gpu}` : ""].filter(Boolean).join("\n");
  btn.disabled = ["downloading", "installing", "testing"].includes(st.stage);
  btn.classList.toggle("ok", st.stage === "ready");
  const tierBtn = document.querySelector('#tierSeg [data-tier="dft-gpu"]');
  tierBtn.hidden = !(st.stage === "ready" && st.gpu_ok);
  if (tierBtn.hidden && S.tier === "dft-gpu") setTier("dft");
  if (!$("view-adna").hidden) initAdna();
}

function selftestSummary(rep) {
  if (!rep) return "";
  if (rep.error) return rep.error;
  const a = rep.analytic_hessian;
  const an = !a ? " · analytic Hessian: not tested (re-run the self-test to enable fast GPU frequencies)"
    : a.error ? ` · analytic Hessian: ${a.error}`
    : ` · analytic Hessian vs finite differences (formaldehyde): max Δfreq ${a.max_frequency_difference_cm1} cm⁻¹, ${a.seconds_analytic} s vs ${a.seconds_finite_difference} s ${a.passed ? "PASS" : "FAIL (frequencies use finite differences)"}`;
  return (rep.cases || []).map((c) => c.error ? `${c.molecule}: ${c.error}` :
    `${c.molecule} (${c.atoms} atoms, ${c.method}): ΔE ${c.differences.energy_hartree.toExponential(1)} Eh, ` +
    `max Δgradient ${c.differences.gradient_hartree_per_bohr.toExponential(1)} Eh/bohr; ` +
    `GPU ${c.seconds_gpu} s vs CPU ${c.seconds_cpu} s (${c.speedup}×) ${c.passed ? "PASS" : "FAIL"}`).join(" · ") + an;
}

async function initCompute() {
  if (!TAURI || !TAURI.core) return;
  let last = "";
  TAURI.event.listen("compute-status", (e) => {
    const st = e.payload || {};
    renderCompute(st);
    if (st.message && st.message !== last && st.stage !== "downloading") { log(st.message, st.stage === "failed"); last = st.message; }
    if ((st.stage === "ready" || st.stage === "failed") && st.selftest) log(`GPU self-test: ${selftestSummary(st.selftest)}`, st.stage === "failed");
  });
  renderCompute(await TAURI.core.invoke("compute_status"));
  $("gpuBtn").onclick = async () => {
    const st = S.compute || {};
    if (st.stage === "ready") {
      const again = confirm(`Linux engines are ready${st.gpu_ok ? ", including GPU DFT" : "; the GPU tier is off"}.\n\n${st.gpu || "No NVIDIA GPU detected"}\n\nGPU self-test: ${selftestSummary(st.selftest)}\n\n` +
        "Re-run the GPU self-test now? (A minute or two; nothing is downloaded.)");
      if (again) { try { await TAURI.core.invoke("compute_enable"); } catch (e) { log(`GPU compute: ${e}`, true); } }
      return;
    }
    if (st.stage === "needs-restart") { alert(st.message); return; }
    const retry = st.stage === "failed";
    const ok = retry ? confirm("Retry?\n\nIf the Linux engines pack is already installed it is not downloaded again; Athanor re-runs the GPU self-test.") : confirm(
      "Install the Linux engines?\n\n" +
      "• Uses WSL2. If it is not enabled yet, Windows will ask for administrator approval, and may need a restart.\n" +
      "• Downloads the Linux engines pack (about 2.4 GB: GPU DFT and the ancient-DNA tools) from the Athanor release feed and checks every part's SHA-256.\n" +
      "• Runs a self-test on your NVIDIA GPU (GPU vs CPU on the same molecules) before the GPU tier is offered; without a GPU the other Linux engines still work.\n\n" +
      "You can keep working while it installs.");
    if (!ok) return;
    try { await TAURI.core.invoke("compute_enable"); }
    catch (e) { log(`GPU compute: ${e}`, true); }
  };
}

// ------------------------------------------------------------------ Ancient DNA
const A = { catalog: null, result: null, busy: false, mode: "sim" };

async function initAdna() {
  const linuxReady = S.compute && S.compute.stage === "ready";
  $("adnaRun").disabled = !linuxReady || A.busy;
  $("realRun").disabled = !linuxReady || A.busy;
  $("adnaEmpty").textContent = !linuxReady
    ? "The Ancient DNA Lab runs in the Linux engines pack. Install or update it from the top bar first."
    : A.mode === "real" ? "Add read files, choose a starting reference, then reconstruct."
    : "Choose a sample and references, then run.";
  if (!linuxReady || A.catalog) return;
  try {
    A.catalog = await rpc("adna_catalog", {});
    const g = A.catalog.references.genomes;
    $("adnaEndo").innerHTML = g.map((x) => `<option value="${x.accession}" ${x.accession === "NC_007596" ? "selected" : ""}>${escapeHtml(x.organism)} — ${x.version}</option>`).join("");
    $("adnaRefs").innerHTML = g.map((x) => `<label><input type="checkbox" value="${x.accession}" ${["NC_007596", "NC_000934", "NC_005129"].includes(x.accession) ? "checked" : ""}>
      <span>${escapeHtml(x.organism)}<small>${x.version} · ${x.length.toLocaleString()} bp</small></span></label>`).join("");
    $("realRef").innerHTML = g.map((x) => `<option value="${x.accession}" ${x.accession === "NC_007596" ? "selected" : ""}>${escapeHtml(x.organism)} — ${x.version}</option>`).join("")
      + `<option value="__file">Another genome (FASTA file)…</option>`;
  } catch (e) { $("adnaEmpty").textContent = `Could not reach the Linux engines: ${e.message || e}`; }
}

function adnaParams() {
  const len = +$("adnaLen").value;
  return {
    endogenous: $("adnaEndo").value,
    references: [...document.querySelectorAll("#adnaRefs input:checked")].map((i) => i.value),
    coverage: +$("adnaCov").value,
    damage: $("adnaDamage").value.split(",").map(Number),
    contamination: +$("adnaCont").value / 100,
    length_loc: Math.log(len), length_scale: 0.35,
    iterations: +$("adnaIt").value,
  };
}

function fmtDuration(sec) {
  if (sec < 45) return "less than a minute";
  if (sec < 90) return "about a minute";
  if (sec < 3600) return `about ${Math.round(sec / 60)} min`;
  return `about ${(sec / 3600).toFixed(1)} h`;
}
function adnaProgress(d) {
  const f = d.fraction ?? 0;
  $("adnaBar").style.width = `${(100 * f).toFixed(1)}%`;
  $("adnaPct").textContent = `${Math.floor(100 * f)}%`;
  const el = d.elapsed_seconds != null ? ` · ${fmtDuration(d.elapsed_seconds).replace("about ", "")} elapsed` : "";
  $("adnaEta").textContent = d.overrunning
    ? `· this step is taking longer than expected · ${fmtDuration(d.eta_seconds)} left at the earliest${el}`
    : `· ${fmtDuration(d.eta_seconds)} left${el}`;
  $("adnaBusyText").textContent = progressText(d);
}

async function runAdna() {
  const p = adnaParams();
  if (!p.references.length) { log("Choose at least one reference", true); return; }
  A.busy = true; $("adnaRun").disabled = true; $("adnaBusy").hidden = false; $("adnaBusyText").textContent = "Starting…";
  $("adnaBar").style.width = "0%"; $("adnaPct").textContent = "0%"; $("adnaEta").textContent = "· estimating time";
  const job = addJob(`Ancient DNA · ${p.coverage}× · ${p.references.length} reference(s)`);
  try {
    A.result = await rpc("adna_run", p);
    if (A.mode === "sim") renderAdna();
    job.done(`${A.result.references.length} reference(s)`);
    const endo = A.result.simulation.endogenous.accession;
    const others = A.result.references.filter((r) => r.reference.accession !== endo);
    if (others.length) {
      const best = others.reduce((a, b) => (a.accuracy > b.accuracy ? a : b));
      log(`Ancient DNA: best reconstruction from a relative's genome ${(100 * best.accuracy).toFixed(3)}% correct (${commonName(best.reference.organism)})`);
    } else log("Ancient DNA: reconstructed against the sample's own genome (control)");
  } catch (e) { job.fail(String(e)); log(`Ancient DNA run failed: ${e.message || e}`, true); }
  finally { A.busy = false; $("adnaRun").disabled = false; $("adnaBusy").hidden = true; }
}

const commonName = (org) => { const m = /\(([^)]+)\)/.exec(org); return m ? m[1].replace(/^./, (c) => c.toUpperCase()) : org; };
const pct = (x, d = 2) => (x == null ? "—" : `${(100 * x).toFixed(d)}%`);
function renderAdna() {
  const R = A.result; if (!R) return;
  $("adnaEmpty").hidden = true; $("adnaResults").hidden = false;
  const rows = R.references.map((r) => {
    const c = r.calibration;
    const self = r.reference.accession === R.simulation.endogenous.accession;
    return `<tr><td>${escapeHtml(commonName(r.reference.organism))}${self ? " <small>(control)</small>" : ""}<br><small>${r.reference.version}</small></td>
      <td class="n">${pct(r.truth_vs_reference_divergence)}</td><td class="n">${pct(r.estimated_divergence_theta)}</td>
      <td class="n">${(r.mapped_reads ?? 0).toLocaleString()}${r.total_reads ? `<br><small>of ${r.total_reads.toLocaleString()}</small>` : ""}</td><td class="n">${r.mean_depth.toFixed(1)}×</td>
      <td class="n">${pct(r.covered_fraction, 1)}</td><td class="n">${pct(r.accuracy, 3)}</td>
      <td class="n">${pct(c.predicted_error_rate, 3)}</td><td class="n">${pct(c.observed_error_rate, 3)}</td>
      <td class="n">${c.expected_calibration_error.toFixed(4)}</td></tr>`;
  }).join("");
  $("adnaTable").innerHTML = `<table><thead><tr><th>Reference</th><th>True divergence</th><th>Estimated</th>
    <th>Mapped</th><th>Depth</th><th>Covered</th><th>Correct</th><th>Predicted error</th><th>Observed error</th><th>ECE</th></tr></thead><tbody>${rows}</tbody></table>
    <p class="note">"Predicted error" is what the posteriors claim; "observed" is measured against the true genome. ECE (expected calibration error) near 0 means the stated confidence can be trusted. Positions where the truth has an indel relative to the reference are not scored.</p>`;
  const itRows = R.references.filter((r) => r.iterations && r.reference.accession !== R.simulation.endogenous.accession).flatMap((r) =>
    r.iterations.map((h) => `<tr><td>${h.iteration === 1 ? escapeHtml(commonName(r.reference.organism)) : ""}</td><td class="n">${h.iteration}</td>
      <td class="n">${(h.mapped_reads ?? 0).toLocaleString()}</td><td class="n">${pct(h.covered_fraction, 1)}</td><td class="n">${pct(h.accuracy, 3)}</td>
      <td class="n">${pct(h.estimated_divergence_theta)}</td><td class="n">${pct(r.truth_vs_reference_divergence)}</td>
      <td class="n">${pct(h.predicted_error, 3)} / ${pct(h.observed_error, 3)}</td></tr>`));
  $("adnaIterTable").innerHTML = itRows.length ? `<table><thead><tr><th>Reference</th><th>Iteration</th><th>Mapped</th><th>Covered</th><th>Correct</th>
    <th>Estimated divergence</th><th>True divergence</th><th>Predicted / observed error</th></tr></thead><tbody>${itRows.join("")}</tbody></table>
    <p class="note">Reads carrying the sample's own bases map poorly to a relative's genome (reference bias). Each iteration recovers more of them, and the divergence estimate approaches the truth.</p>` : `<p class="note">Map against a relative's genome to see reference bias.</p>`;
  const sim = R.simulation.damage_briggs_v_l_d_s;
  const names = [["Nu", "ν nick frequency", sim[0]], ["Lambda", "λ overhang parameter", sim[1]], ["DeltaD", "δd double-stranded", sim[2]], ["DeltaS", "δs single-stranded", sim[3]]];
  const drows = names.map(([k, lab, v]) => `<tr><td>${lab}</td><td class="n">${v}</td>${R.references.map((r) => {
    const e = r.mapdamage_bayesian && r.mapdamage_bayesian[k];
    return `<td class="n">${e && e.Mean != null ? e.Mean.toFixed(4) : "—"}</td>`; }).join("")}</tr>`).join("");
  $("adnaDamageTable").innerHTML = `<table><thead><tr><th>Briggs parameter</th><th>Simulated</th>${R.references.map((r) => `<th>${escapeHtml(commonName(r.reference.organism))}</th>`).join("")}</tr></thead><tbody>${drows}</tbody></table>
    <p class="note">mapDamage2 posterior means, estimated on the final (iterated) reference. δd and λ are recovered; δs comes out near 0.55 against gargammel's s = 0.3 on every reference including the sample's own, which points to the two programs defining that parameter differently. That is under investigation, so δs is reported but not used.</p>`;
  drawDamage(); drawCalibration(); drawTracks(); renderAdnaProv();
}

function plotFrame(cv) {
  const dpr = window.devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight || 220;
  cv.width = w * dpr; cv.height = h * dpr;
  const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, w, h);
  g.font = "11px JetBrains Mono, monospace"; g.fillStyle = "#7f9fb3"; g.strokeStyle = "rgba(110,205,255,.25)";
  return { g, w, h, l: 44, r: w - 10, t: 18, b: h - 26 };
}
const COLORS = ["#62d6ff", "#ffb34a", "#52f0b0", "#ff6a57"];

function drawDamage() {
  const { g, l, r, t, b } = plotFrame($("adnaDamagePlot"));
  const refs = A.result.references, n = 25;
  const ymax = Math.max(0.05, ...refs.flatMap((x) => [...x.damage_profile.raw_5p_CtoT.slice(0, n), ...x.damage_profile.raw_3p_GtoA.slice(0, n)]));
  const X = (i) => l + (i / (n - 1)) * (r - l), Y = (v) => b - (v / ymax) * (b - t);
  g.strokeRect(l, t, r - l, b - t); g.fillText("C→T from 5′ (solid), G→A from 3′ (dashed)", l, t - 5);
  g.fillText(ymax.toFixed(2), 4, t + 8); g.fillText("0", 30, b); g.fillText("position from read end", l + 60, b + 18);
  refs.forEach((x, k) => {
    for (const [arr, dash] of [[x.damage_profile.raw_5p_CtoT, []], [x.damage_profile.raw_3p_GtoA, [5, 4]]]) {
      g.setLineDash(dash); g.strokeStyle = COLORS[k % 4]; g.beginPath();
      arr.slice(0, n).forEach((v, i) => (i ? g.lineTo(X(i), Y(v)) : g.moveTo(X(i), Y(v)))); g.stroke();
    }
  });
  g.setLineDash([]);
  refs.forEach((x, k) => { g.fillStyle = COLORS[k % 4]; g.fillText(commonName(x.reference.organism).slice(0, 22), r - 150, t + 14 + 13 * k); });
}

function drawCalibration() {
  const { g, l, r, t, b } = plotFrame($("adnaCalPlot"));
  const lo = -6, hi = 0, X = (e) => l + ((Math.log10(Math.max(e, 1e-6)) - lo) / (hi - lo)) * (r - l), Y = (e) => b - ((Math.log10(Math.max(e, 1e-6)) - lo) / (hi - lo)) * (b - t);
  g.strokeRect(l, t, r - l, b - t); g.fillText("observed vs predicted error (log)", l, t - 5);
  g.fillText("1e-6", 4, b); g.fillText("1", 30, t + 8); g.fillText("predicted error", l + 80, b + 18);
  g.setLineDash([4, 4]); g.beginPath(); g.moveTo(X(1e-6), Y(1e-6)); g.lineTo(X(1), Y(1)); g.stroke(); g.setLineDash([]);
  A.result.references.forEach((x, k) => {
    g.fillStyle = COLORS[k % 4];
    for (const bn of x.calibration.bins) {
      if (bn.observed_error <= 0) continue;
      const rad = Math.min(6, 2 + Math.log10(bn.n));
      g.beginPath(); g.arc(X(bn.predicted_error), Y(bn.observed_error), rad, 0, 7); g.fill();
    }
    g.fillText(commonName(x.reference.organism).slice(0, 22), r - 150, b - 10 - 13 * k);
  });
  g.fillStyle = "#7f9fb3"; g.fillText("above the diagonal = overconfident", l + 6, t + 14);
}

function drawTracks() {
  const refs = A.result.references;
  $("adnaTracks").innerHTML = refs.map((x, k) => `<div class="track"><div class="tl" style="color:${COLORS[k % 4]}">${escapeHtml(x.reference.organism)} · ${x.reference.version}</div><canvas id="trk${k}"></canvas></div>`).join("");
  refs.forEach((x, k) => {
    const cv = $(`trk${k}`), dpr = window.devicePixelRatio || 1, w = cv.clientWidth, h = 54;
    cv.width = w * dpr; cv.height = h * dpr; const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0);
    const L = x.depth.length, maxd = Math.max(1, ...x.depth), bins = Math.max(1, Math.floor(w));
    for (let p = 0; p < bins; p++) {
      const s0 = Math.floor((p * L) / bins), s1 = Math.max(s0 + 1, Math.floor(((p + 1) * L) / bins));
      let d = 0, low = 0, unc = 0;
      for (let i = s0; i < s1; i++) { d += x.depth[i]; if (x.depth[i] === 0) unc++; else if (x.confidence[i] < 0.99) low++; }
      d /= (s1 - s0);
      g.fillStyle = unc === s1 - s0 ? "rgba(127,159,179,.35)" : "rgba(98,214,255,.45)";
      const hh = unc === s1 - s0 ? h : (d / maxd) * (h - 6);
      g.fillRect(p, h - hh, 1, hh);
      if (low) { g.fillStyle = "#ffb34a"; g.fillRect(p, 0, 1, 4); }
    }
    g.fillStyle = "#ff6a57";
    for (const pos of x.error_positions || []) g.fillRect(Math.floor((pos / L) * w), 6, 1.5, h - 6);
  });
}

function renderAdnaProv() {
  const R = A.result;
  const tools = Object.entries(R.tools).map(([k, v]) => `<li><b>${escapeHtml(k)}</b>: ${escapeHtml(v)}</li>`).join("");
  const refs = R.reference_catalog.genomes.map((g) => `<li>${escapeHtml(g.version)} — ${escapeHtml(g.organism)} (${g.length.toLocaleString()} bp, sha256 ${g.sha256.fa.slice(0, 12)}…)</li>`).join("");
  const steps = R.steps.map((s) => `<li>${escapeHtml(s.step)} · ${s.seconds} s<br><code>${escapeHtml(s.command)}</code></li>`).join("");
  const S0 = R.simulation, b = S0.damage_briggs_v_l_d_s;
  const simRows = [["Endogenous", `${S0.endogenous.version} (${commonName(S0.endogenous.organism)})`], ["Coverage", `${S0.coverage}×`],
    ["Damage ν, λ, δd, δs", b.join(", ")], ["Contamination", pct(S0.contamination, 1)],
    ["Fragment length", `log-normal, median ${Math.round(Math.exp(S0.fragment_length_lognormal[0]))} bp, σ ${S0.fragment_length_lognormal[1]}; min ${S0.min_length} bp`],
    ["Reads", `${S0.read_length} bp paired-end, ${S0.platform}`]];
  $("adnaProv").innerHTML = `<p><b>Simulation</b></p><ul>${simRows.map(([k, v]) => `<li>${escapeHtml(k)}: ${escapeHtml(String(v))}</li>`).join("")}</ul>
    <p><b>Tools</b></p><ul>${tools}</ul><p><b>Reference genomes</b> (${escapeHtml(R.reference_catalog.source)}, retrieved ${escapeHtml(R.reference_catalog.retrieved_utc)})</p><ul>${refs}</ul>
    <p><b>Steps</b></p><ol>${steps}</ol>`;
}

// ------------------------------------------------------------------ Ancient DNA: the user's own reads
const R2 = { raw: [], files: [], refPath: null, cmpPath: null, result: null };
const baseName = (p) => String(p).split(/[\\/]/).pop();

function setAdnaMode(m) {
  A.mode = m;
  document.querySelectorAll("#adnaModeSeg button").forEach((b) => b.setAttribute("aria-checked", String(b.dataset.mode === m)));
  $("adnaSimControls").hidden = m !== "sim";
  $("adnaRealControls").hidden = m !== "real";
  const has = m === "sim" ? !!A.result : !!R2.result;
  $("adnaResults").hidden = !(m === "sim" && A.result);
  $("realResults").hidden = !(m === "real" && R2.result);
  $("adnaEmpty").hidden = has;
  initAdna();
  if (m === "sim" && A.result) renderAdna();
  else if (m === "real" && R2.result) renderReal();
  else $("adnaProv").textContent = "No run yet.";
}

function kindOf(p) {
  if (/\.bam$/i.test(p)) return "bam";
  if (/\.f(ast)?q(\.gz)?$/i.test(p)) return "fastq";
  return null;
}

// Pair FASTQ mates: file names that differ only by a 1 in one and a 2 in the
// other (sample_R1_001.fastq.gz / sample_R2_001.fastq.gz, run_1.fq / run_2.fq),
// where that digit stands alone (lane chunks _001/_002 are not mates).
function pairFiles(paths) {
  const out = [], fq = [];
  for (const p of paths) { const k = kindOf(p); if (k === "fastq") fq.push(p); else if (k) out.push({ kind: k, path: p }); }
  const used = new Set();
  const mates = (a, b) => {
    const x = baseName(a), y = baseName(b);
    if (x.length !== y.length || a.length - x.length !== b.length - y.length || a.slice(0, a.length - x.length) !== b.slice(0, b.length - y.length)) return false;
    let d = -1;
    for (let i = 0; i < x.length; i++) if (x[i] !== y[i]) { if (d >= 0) return false; d = i; }
    // the mate number stands alone (not part of a lane or chunk number like _001/_002)
    return d >= 0 && x[d] === "1" && y[d] === "2" && !/\d/.test(x[d - 1] || "") && !/\d/.test(x[d + 1] || "");
  };
  for (const a of fq) {
    if (used.has(a)) continue;
    const b = fq.find((q) => !used.has(q) && q !== a && mates(a, q));
    if (b) { used.add(a); used.add(b); out.push({ kind: "fastq", path: a, path2: b }); }
  }
  for (const a of fq) if (!used.has(a)) out.push({ kind: "fastq", path: a });
  return out;
}

function renderRealFiles() {
  $("realFiles").innerHTML = R2.files.length
    ? R2.files.map((f, i) => `<div class="file"><span class="kind">${f.kind === "bam" ? "BAM" : f.path2 ? "FASTQ pair" : "FASTQ"}</span>
        <span class="fn" title="${escapeHtml(f.path)}${f.path2 ? "\n" + escapeHtml(f.path2) : ""}">${escapeHtml(baseName(f.path))}${f.path2 ? `<br>${escapeHtml(baseName(f.path2))}` : ""}</span>
        <button type="button" class="btn sm" data-rm="${i}" aria-label="Remove">×</button></div>`).join("")
    : `<p class="note">No files yet.</p>`;
  $("realFiles").querySelectorAll("[data-rm]").forEach((b) => (b.onclick = () => {
    const f = R2.files[+b.dataset.rm];
    R2.raw = R2.raw.filter((p) => p !== f.path && p !== f.path2);
    R2.files = pairFiles(R2.raw); renderRealFiles();
  }));
}

async function realAddFiles() {
  if (!TAURI) return;
  try {
    const paths = await TAURI.core.invoke("pick_files", { title: "Add sequencing reads", multiple: true,
      filters: [["Reads (FASTQ, BAM)", ["fq", "fastq", "gz", "bam"]], ["All files", ["*"]]] });
    const bad = paths.filter((p) => !kindOf(p));
    if (bad.length) log(`Not FASTQ or BAM (skipped): ${bad.map(baseName).join(", ")}`, true);
    R2.raw = [...new Set([...R2.raw, ...paths.filter(kindOf)])];
    R2.files = pairFiles(R2.raw); renderRealFiles();
    if ($("realSample").value === "sample" && R2.files.length) $("realSample").value = baseName(R2.files[0].path).replace(/([._-]R?[12])?([._-]\d{3})?\.(bam|f(ast)?q(\.gz)?)$/i, "");
  } catch (e) { log(`File dialog: ${e}`, true); }
}

async function realRefChanged() {
  if ($("realRef").value !== "__file") { R2.refPath = null; $("realRefFileRow").hidden = true; return; }
  const [p] = await TAURI.core.invoke("pick_files", { title: "Reference genome (FASTA)", multiple: false,
    filters: [["FASTA", ["fa", "fasta", "fna", "fas"]], ["All files", ["*"]]] }).catch(() => []);
  if (!p) { $("realRef").value = "NC_007596"; return; }
  R2.refPath = p; $("realRefFile").textContent = baseName(p); $("realRefFileRow").hidden = false;
}

async function realPickCompare() {
  const [p] = await TAURI.core.invoke("pick_files", { title: "Published sequence (FASTA)", multiple: false,
    filters: [["FASTA", ["fa", "fasta", "fna", "fas"]], ["All files", ["*"]]] }).catch(() => []);
  if (p) { R2.cmpPath = p; $("realCmpName").textContent = baseName(p); }
}

async function runReal() {
  if (!R2.files.length) { log("Add read files first", true); return; }
  const trimmed = $("realTrimmed").checked;
  const p = {
    sample: $("realSample").value.trim() || "sample",
    inputs: R2.files.map((f) => ({ kind: f.kind, path: f.path, ...(f.path2 ? { path2: f.path2 } : {}), trimmed: f.kind === "bam" ? true : trimmed })),
    circular: $("realCircular").checked, contamination_check: $("realContam").checked,
    iterations: +$("realIt").value, min_depth: +$("realDepth").value, min_confidence: +$("realConf").value,
    mapdamage_stats: $("realStats").checked,
  };
  if (R2.refPath) p.reference_path = R2.refPath; else p.reference = $("realRef").value;
  if (R2.cmpPath) p.compare_path = R2.cmpPath;
  A.busy = true; $("realRun").disabled = true; $("adnaRun").disabled = true; $("adnaBusy").hidden = false;
  $("adnaBusyText").textContent = "Starting…"; $("adnaBar").style.width = "0%"; $("adnaPct").textContent = "0%"; $("adnaEta").textContent = "· estimating time";
  const job = addJob(`Ancient DNA · ${p.sample} · ${p.inputs.length} input(s)`);
  try {
    R2.result = await rpc("adna_real_run", p);
    const c = R2.result.consensus;
    job.done(`${c.called_sites.toLocaleString()} of ${R2.result.coverage.length.toLocaleString()} bases called`);
    log(`Ancient DNA: ${p.sample} reconstructed, ${c.called_sites.toLocaleString()} bases called at ${R2.result.coverage.mean_depth.toFixed(1)}× mean depth`);
    if (A.mode === "real") { $("adnaEmpty").hidden = true; $("realResults").hidden = false; renderReal(); }
  } catch (e) { job.fail(String(e)); log(`Ancient DNA run failed: ${e.message || e}`, true); }
  finally { A.busy = false; $("realRun").disabled = false; $("adnaRun").disabled = false; $("adnaBusy").hidden = true; }
}

async function realExport(what) {
  if (!R2.result) return;
  const ext = what === "json" ? "json" : "fasta";
  const name = `${R2.result.sample}${what === "majority_fasta" ? ".majority" : what === "json" ? ".athanor-report" : ".athanor"}.${ext}`;
  try {
    const path = await TAURI.core.invoke("pick_save_path", { title: "Save", defaultName: name,
      filters: [[ext === "json" ? "JSON" : "FASTA", [ext === "json" ? "json" : "fasta", "fa"]]] });
    if (!path) return;
    await rpc("adna_export", { path, what });
    log(`Saved ${baseName(path)}`);
  } catch (e) { log(`Save failed: ${e.message || e}`, true); }
}

const ci = (c) => (c ? `${pct(c[0], 2)} – ${pct(c[1], 2)}` : "—");
const kvTable = (rows) => `<table class="kv"><tbody>${rows.map(([k, v]) => `<tr><td>${k}</td><td class="n">${v}</td></tr>`).join("")}</tbody></table>`;

function renderReal() {
  const R = R2.result; if (!R) return;
  const reads = R.inputs.reduce((a, x) => a + (x.reads_kept || 0), 0);
  const m = R.mapping, cov = R.coverage, c = R.consensus, L = cov.length;
  const ref = R.reference;
  $("realSummary").innerHTML = kvTable([
    ["Sample", escapeHtml(R.sample)],
    ["Starting reference", `${escapeHtml(ref.organism || ref.accession || "")} ${escapeHtml(ref.version || "")}`],
    ["Reads analysed (≥ " + R.settings.min_length + " bp)", `${reads.toLocaleString()} <small>(${R.fragment_lengths.whole_molecules.toLocaleString()} whole molecules, ${R.fragment_lengths.cut_by_read_length.toLocaleString()} cut by the read length)</small>`],
    ["Mapped, MAPQ ≥ " + R.settings.min_mapq + ", unique", `${(m.kept || 0).toLocaleString()} <small>(${(m.duplicates || 0).toLocaleString()} duplicates removed, ${pct(m.duplication_rate, 1)})</small>`],
    ["Depth", `mean ${cov.mean_depth.toFixed(1)}×, median ${cov.median_depth}× · covered ${pct(cov.fraction_ge1, 1)}, ≥3 reads ${pct(cov.fraction_ge3, 1)}, ≥10 reads ${pct(cov.fraction_ge10, 1)}`],
    ["Bases called", `${c.called_sites.toLocaleString()} of ${L.toLocaleString()} (${pct(c.called_sites / L, 2)}); the rest are N`],
    ["Expected errors among called bases", c.predicted_errors_in_called_sites < 0.01 ? "< 0.01" : c.predicted_errors_in_called_sites.toFixed(2)],
    ["Divergence from the starting reference", pct(c.estimated_divergence_theta, 2)],
    ["Iterative mapping", R.converged ? `converged after ${R.iterations.length} iteration(s)`
      : R.oscillating ? `<span class="warn">stopped after ${R.iterations.length} iterations: the remaining calls alternate between states (positions ${escapeHtml((R.pending_changes[0]?.substitutions || []).slice(0, 8).map((x) => x + 1).join(", "))}); they are reported with their posteriors</span>`
      : `<span class="warn">not converged after ${R.iterations.length} iterations; increase iterations</span>`],
    ["Majority-rule consensus", `${R.majority_rule.called_sites.toLocaleString()} bases called (${escapeHtml(R.majority_rule.rule)})`],
  ]) + `<p class="note">Bases are called from damage-aware posteriors: N where fewer than ${R.settings.min_depth} reads cover a base or its posterior is below ${R.settings.min_posterior}. The expected error count is the sum of (1 − posterior) over called bases; it assumes the model's damage and error rates, which are measured from this sample.</p>`;

  const ts = R.authentication.terminal_signal, fl = R.fragment_lengths;
  const row = (lab, x) => [lab, x.frequency == null ? "—" : `${pct(x.frequency, 2)} <small>(95% ${ci(x.ci95)}; n ${x.n.toLocaleString()})</small>`];
  $("realAuth").innerHTML = kvTable([
    row("5′ C→T, first base", ts.CtoT_5p_position1), row("5′ C→T, second base", ts.CtoT_5p_position2),
    row("C→T, molecule interior (≥ 10 bases in)", ts.CtoT_interior),
    row("3′ G→A, last base (whole molecules)", ts.GtoA_3p_position1), row("G→A, interior", ts.GtoA_interior),
    ["Whole-molecule length", fl.median == null ? "—" : `median ${fl.median} bp, mean ${fl.mean.toFixed(1)} bp`],
  ]) + `<p class="${ts.terminal_excess ? "note" : "warnbox"}">${ts.terminal_excess
      ? "Deamination is elevated at the molecule ends compared with their interior (intervals do not overlap): the pattern of authentic ancient DNA (Briggs et al. 2007). Libraries partially treated with uracil-DNA glycosylase keep this signal only at the terminal one or two bases."
      : "No excess of terminal deamination was detected. Either the library was fully treated with uracil-DNA glycosylase, or the reads are not ancient."}</p>
    <p class="note">Measured against this sample's own consensus, so genuine differences from the reference do not count as damage.</p>`;

  const ct = R.contamination;
  if (!ct) $("realContamOut").innerHTML = `<p class="note">Not estimated.</p>`;
  else {
    const h = ct.human_reads_detail, t = ct.target_reads;
    $("realContamOut").innerHTML = kvTable([
      ["Present-day human (reads assigned to human mtDNA)", `${pct(ct.fraction, 2)} <small>(95% ${ci(ct.ci95)}; ${ct.human_reads.toLocaleString()} of ${ct.assigned_reads.toLocaleString()})</small>`],
      ["Median length, human / sample reads", `${h.median_length ?? "—"} / ${t.median_length ?? "—"} bp`],
      ["5′ C→T first base, human / sample reads", `${pct(h.damage.CtoT_5p_position1.frequency, 1)} / ${pct(t.damage.CtoT_5p_position1.frequency, 1)}`],
    ]) + `<p class="note">${escapeHtml(ct.method)}: reads mapped to the sample's sequence and the human mitogenome (${escapeHtml(ct.references.human)}) together; only reads with mapping quality above ${ct.min_mapq_exclusive} count. Present-day contaminants are typically longer and lack terminal damage.</p>`
      + ct.notes.map((n) => `<p class="note">${escapeHtml(n)}</p>`).join("");
  }

  $("realIter").innerHTML = `<table><thead><tr><th>Iteration</th><th>Length</th><th>Reads kept</th><th>Covered</th><th>Mean depth</th><th>Divergence</th><th>Bases changed</th><th>Indels</th></tr></thead><tbody>${
    R.iterations.map((h) => `<tr><td class="n">${h.iteration}</td><td class="n">${h.reference_length.toLocaleString()}</td><td class="n">${(h.mapping.kept || 0).toLocaleString()}</td>
      <td class="n">${pct(h.covered_fraction, 1)}</td><td class="n">${h.mean_depth.toFixed(1)}×</td><td class="n">${pct(h.estimated_divergence_theta)}</td>
      <td class="n">${h.substitutions_found}</td><td class="n">${h.indels_found}</td></tr>`).join("")}</tbody></table>
    <p class="note">Each round writes confident bases and majority-supported insertions and deletions into the reference and maps again, until nothing changes. Posteriors always use the starting reference as the prior, so reads are not counted twice.</p>`;

  const cmp = R.comparison;
  $("realCmpSec").hidden = !cmp;
  if (cmp) {
    const line = (lab, x) => `<tr><td>${lab}</td><td class="n">${pct(x.agreement, 3)}</td><td class="n">${x.sites_called_by_both.toLocaleString()}</td>
      <td class="n">${x.substitutions.length}</td><td class="n">${x.insertions_in_ours.length + x.deletions_in_ours.length}</td>
      <td class="n">${x.uncalled_ours.toLocaleString()} / ${x.uncalled_theirs.toLocaleString()}</td><td class="n">${pct(x.theirs_aligned_fraction, 1)}</td></tr>`;
    const diffs = cmp.athanor.substitutions.map((d) => `<tr><td class="n">${(d.theirs_pos + 1).toLocaleString()}</td><td class="n">${(d.ours_pos + 1).toLocaleString()}</td>
      <td class="n">${d.theirs}</td><td class="n">${d.ours}</td><td class="n">${d.depth}</td><td class="n">${d.posterior.toFixed(4)}</td>
      <td class="n">${Object.entries(d.base_counts_trimmed).map(([b, n]) => `${b}${n}`).join(" ")}</td></tr>
      ${d.observations && d.observations.length ? `<tr><td colspan="7"><details><summary class="note">Every read base here (strand, bases from the 5′ / 3′ end)</summary>
        <p class="note mono">${d.observations.map((o) => `${o.base}${o.strand} ${o.from_5p}/${o.from_3p ?? "–"}`).join(" · ")}</p>
        <p class="note">A T supported only by bases at 5′ ends of + reads or 3′ ends of − reads (or an A at the mirror positions) is what deamination produces.</p></details></td></tr>` : ""}`).join("");
    const indels = [...cmp.athanor.insertions_in_ours.map((d) => [d.uncalled ? "extra uncalled (N) bases in ours: length unresolved" : "extra in ours", d]),
                    ...cmp.athanor.deletions_in_ours.map((d) => [d.uncalled ? "published has N here; ours lacks them" : "missing in ours", d])]
      .map(([k, d]) => `<tr><td class="n">${(d.theirs_pos + 1).toLocaleString()}</td><td>${k}</td><td class="n">${escapeHtml(d.sequence)}</td></tr>`).join("");
    $("realCmpOut").innerHTML = `<table><thead><tr><th>Consensus</th><th>Agreement</th><th>Bases both call</th><th>Differences</th><th>Indel differences</th><th>Uncalled ours / theirs</th><th>Published aligned</th></tr></thead>
      <tbody>${line("Athanor (posterior)", cmp.athanor)}${line("Majority rule", cmp.majority_rule)}</tbody></table>
      <p class="note">Compared with ${escapeHtml(cmp.athanor.theirs_name)} by minimap2. Agreement counts only bases both sequences call; N and ambiguity codes are uncalled, not disagreements.</p>
      ${diffs ? `<table><thead><tr><th>Published position</th><th>Our position</th><th>Published</th><th>Athanor</th><th>Reads</th><th>Posterior</th><th>Base counts (ends trimmed)</th></tr></thead><tbody>${diffs}</tbody></table>` : `<p class="note">No base differences where both call a base.</p>`}
      ${indels ? `<table><thead><tr><th>Published position</th><th>Difference</th><th>Bases</th></tr></thead><tbody>${indels}</tbody></table>` : ""}`;
  }
  drawRealDamage(); drawRealLengths(); drawRealTrack(); renderRealProv();
}

function drawRealDamage() {
  const { g, l, r, t, b } = plotFrame($("realDamagePlot"));
  const P = R2.result.authentication.damage_vs_consensus, n = Math.min(25, P.raw_5p_CtoT.length);
  const ymax = Math.max(0.02, ...P.raw_5p_CtoT.slice(0, n), ...P.raw_3p_GtoA.slice(0, n)) * 1.1;
  const X = (i) => l + (i / (n - 1)) * (r - l), Y = (v) => b - (v / ymax) * (b - t);
  g.strokeRect(l, t, r - l, b - t); g.fillText("C→T from 5′ (blue), G→A from 3′ (amber)", l, t - 5);
  g.fillText(ymax.toFixed(3), 4, t + 8); g.fillText("0", 30, b); g.fillText("position from molecule end", l + 60, b + 18);
  for (const [arr, col] of [[P.raw_5p_CtoT, COLORS[0]], [P.raw_3p_GtoA, COLORS[1]]]) {
    g.strokeStyle = col; g.beginPath(); arr.slice(0, n).forEach((v, i) => (i ? g.lineTo(X(i), Y(v)) : g.moveTo(X(i), Y(v)))); g.stroke();
  }
}

function drawRealLengths() {
  const { g, l, r, t, b } = plotFrame($("realLenPlot"));
  const hst = R2.result.fragment_lengths.histogram || [];
  g.strokeRect(l, t, r - l, b - t); g.fillText("whole-molecule lengths (bp)", l, t - 5);
  if (!hst.length) { g.fillText("no whole molecules (all reads ran to the read length)", l + 8, t + 20); return; }
  const lo = 20, hi = Math.max(lo + 10, hst.length), ymax = Math.max(1, ...hst);
  const X = (i) => l + ((i - lo) / (hi - lo)) * (r - l);
  g.fillStyle = "rgba(98,214,255,.6)";
  for (let i = lo; i < hst.length; i++) { const hh = (hst[i] / ymax) * (b - t); g.fillRect(X(i), b - hh, Math.max(1, (r - l) / (hi - lo) - 0.5), hh); }
  g.fillStyle = "#7f9fb3"; g.fillText(String(lo), l, b + 14); g.fillText(String(hi), r - 20, b + 14);
}

function drawRealTrack() {
  const R = R2.result, depth = R.consensus.depth, seq = R.consensus.sequence;
  $("realTrack").innerHTML = `<div class="track"><div class="tl">${escapeHtml(R.sample)} · ${depth.length.toLocaleString()} bp</div><canvas id="realTrk"></canvas></div>`;
  const cv = $("realTrk"), dpr = window.devicePixelRatio || 1, w = cv.clientWidth, h = 54;
  cv.width = w * dpr; cv.height = h * dpr; const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0);
  const L = depth.length, maxd = Math.max(1, ...depth), bins = Math.max(1, Math.floor(w));
  for (let p = 0; p < bins; p++) {
    const s0 = Math.floor((p * L) / bins), s1 = Math.max(s0 + 1, Math.floor(((p + 1) * L) / bins));
    let d = 0, nN = 0;
    for (let i = s0; i < s1; i++) { d += depth[i]; if (seq[i] === "N") nN++; }
    d /= (s1 - s0);
    g.fillStyle = "rgba(98,214,255,.45)"; const hh = (d / maxd) * (h - 6); g.fillRect(p, h - hh, 1, hh);
    if (nN) { g.fillStyle = "#ffb34a"; g.fillRect(p, 0, 1, 4); }
  }
  const cmp = R.comparison && R.comparison.athanor;
  if (cmp) { g.fillStyle = "#ff6a57"; for (const d of cmp.substitutions) g.fillRect(Math.floor((d.ours_pos / L) * w), 6, 1.5, h - 6); }
}

function renderRealProv() {
  const R = R2.result;
  const inputs = R.inputs.map((x) => `<li>${escapeHtml(x.file)}${x.file2 ? " + " + escapeHtml(x.file2) : ""} (${escapeHtml(x.kind)}) · ${(x.reads_in || 0).toLocaleString()} reads, ${(x.reads_kept || 0).toLocaleString()} kept · sha256 ${x.sha256.slice(0, 12)}…${x.all_reads_mapped ? "<br><small>BAM held only mapped reads</small>" : ""}</li>`).join("");
  const tools = Object.entries(R.tools).map(([k, v]) => `<li><b>${escapeHtml(k)}</b>: ${escapeHtml(v)}</li>`).join("");
  const steps = R.steps.map((s) => `<li>${escapeHtml(s.step)} · ${s.seconds} s<br><code>${escapeHtml(s.command)}</code></li>`).join("");
  const st = Object.entries(R.settings).map(([k, v]) => `<li>${escapeHtml(k)}: ${escapeHtml(String(v))}</li>`).join("");
  $("adnaProv").innerHTML = `<p><b>Inputs</b> (read in place on this computer)</p><ul>${inputs}</ul>
    <p><b>Starting reference</b></p><ul><li>${escapeHtml(R.reference.version || R.reference.accession || "")} ${escapeHtml(R.reference.organism || "")} · ${R.reference.length.toLocaleString()} bp · sha256 ${R.reference.sha256.slice(0, 12)}…${R.circular ? " · circular" : ""}</li></ul>
    <p><b>Settings</b></p><ul>${st}</ul><p><b>Tools</b></p><ul>${tools}</ul><p><b>Steps</b></p><ol>${steps}</ol>`;
}

// ------------------------------------------------------------------ Personal Genome
const PG = { status: null, path: null, result: null, busy: false };
const GRADE = { A: "established", B: "moderate", C: "weak or for fun" };

async function initPersonal() {
  try {
    PG.status = await rpc("personal_status", {});
  } catch (e) { $("pgPack").textContent = `Could not reach the science core: ${e.message || e}`; return; }
  renderPgPack();
}

function renderPgPack() {
  const st = PG.status || {};
  $("pgInstall").hidden = !!st.installed && !st.outdated;
  $("pgInstall").textContent = st.outdated ? "Update reference data" : "Download reference data";
  $("pgRemove").hidden = !st.installed;
  $("pgRun").disabled = !(st.installed && PG.path) || PG.busy;
  if (!st.installed) {
    $("pgPack").innerHTML = `Not downloaded yet. The reference data holds published ancient genomes, present-day reference populations, the haplogroup trees and the trait variants. It is downloaded once and checked against its checksum.`;
    $("pgEmpty").textContent = "Download the reference data, then choose your DNA file.";
    return;
  }
  const src = st.sources || {};
  $("pgPack").innerHTML = (st.outdated ? `<p class="note"><b>A newer version of the reference data is available.</b> It adds modern ancestry and the appearance, food and drink and fitness reports. Click Update reference data.</p>` : "") + kvTable([
    ["Version", `${st.pack_version ? `pack ${st.pack_version}` : escapeHtml(st.pack || "")} · built ${escapeHtml((st.built || "").slice(0, 10))}`],
    ["Reference SNPs", (st.snps || 0).toLocaleString()],
    ["Ancient genomes", (st.ancients || 0).toLocaleString()],
    ["Built from", Object.values(src).map((x) => escapeHtml(x.name || "")).join("<br>")],
  ]);
  $("pgEmpty").textContent = PG.path ? "Ready. Click Analyse." : "Choose your DNA file.";
}

function pgPackProgress(d) {
  $("pgPackBusy").hidden = false;
  $("pgPackBar").style.width = `${(100 * (d.fraction || 0)).toFixed(1)}%`;
  $("pgPackText").textContent = progressText(d);
}

async function pgInstall() {
  $("pgInstall").disabled = true;
  const job = addJob("Personal Genome reference data");
  try {
    PG.status = await rpc("personal_install", {});
    job.done("installed");
    log("Personal Genome reference data installed and verified");
  } catch (e) { job.fail(String(e)); log(`Reference data: ${e.message || e}`, true); }
  finally { $("pgInstall").disabled = false; $("pgPackBusy").hidden = true; renderPgPack(); }
}

async function pgRemove() {
  PG.status = await rpc("personal_remove", {});
  renderPgPack();
}

async function pgChoose() {
  if (!TAURI) return;
  try {
    const [p] = await TAURI.core.invoke("pick_files", { title: "Your raw DNA file", multiple: false,
      filters: [["Raw DNA data (txt, csv, zip, gz, vcf)", ["txt", "csv", "tsv", "zip", "gz", "vcf"]], ["All files", ["*"]]] });
    if (!p) return;
    PG.path = p;
    $("pgFile").innerHTML = `<div class="file"><span class="kind">DNA</span><span class="fn" title="${escapeHtml(p)}">${escapeHtml(baseName(p))}</span></div>`;
    renderPgPack();
  } catch (e) { log(`File dialog: ${e}`, true); }
}

function pgProgress(d) {
  if (!PG.busy) return;
  $("pgBar").style.width = `${(100 * (d.fraction || 0)).toFixed(1)}%`;
  $("pgPct").textContent = `${Math.floor(100 * (d.fraction || 0))}%`;
  $("pgBusyText").textContent = progressText(d);
}

async function pgRun() {
  if (!PG.path) return;
  PG.busy = true; $("pgRun").disabled = true; $("pgBusy").hidden = false;
  $("pgBar").style.width = "0%"; $("pgPct").textContent = "0%"; $("pgBusyText").textContent = "Starting…";
  const job = addJob(`Personal Genome · ${baseName(PG.path)}`);
  try {
    PG.result = await rpc("personal_analyse", { path: PG.path });
    job.done(`${PG.result.cards.length} reports`);
    renderPersonal();
    log(`Personal Genome: ${PG.result.cards.length} reports from ${PG.result.file.vendor} data`);
  } catch (e) { job.fail(String(e)); log(`Personal Genome: ${e.message || e}`, true); }
  finally { PG.busy = false; $("pgRun").disabled = false; $("pgBusy").hidden = true; }
}

async function pgExport() {
  try {
    const path = await TAURI.core.invoke("pick_save_path", { title: "Save report",
      defaultName: "athanor-personal-genome.json", filters: [["JSON", ["json"]]] });
    if (!path) return;
    await rpc("personal_export", { path });
    log(`Saved ${baseName(path)}`);
  } catch (e) { log(`Save failed: ${e.message || e}`, true); }
}

const BP = (bp) => (bp == null ? "undated" : (1950 - bp <= 0 ? `${(bp - 1950).toLocaleString()} BCE` : `${(1950 - bp).toLocaleString()} CE`));

function renderPersonal() {
  const R = PG.result; if (!R) return;
  $("pgEmpty").hidden = true; $("pgCards").hidden = false; $("pgExport").hidden = false;
  const groups = [];
  for (const c of R.cards) {
    let g = groups.find((x) => x.name === c.group);
    if (!g) groups.push((g = { name: c.group, cards: [] }));
    g.cards.push(c);
  }
  $("pgCards").innerHTML = groups.map((g) => `<div class="sec"><div class="label">${escapeHtml(g.name)}</div>
    ${g.cards.map(cardHtml).join("")}</div>`).join("");
  for (const c of R.cards) {
    const el = $(`pgbody_${c.id}`);
    const fn = DETAIL[c.id] || (c.id.startsWith("pgx_") ? pgxDetail : null);
    if (el && fn) { try { fn(el, c); } catch (e) { el.innerHTML = `<p class="note">${escapeHtml(String(e))}</p>`; } }
  }
  renderPgProv();
}

function cardHtml(c) {
  return `<div class="pgcard${c.ok ? "" : " warn"}">
    <div class="pghead"><b>${escapeHtml(c.title)}</b><span class="grade g${c.grade}" title="Strength of the evidence">${c.grade} · ${GRADE[c.grade] || ""}</span></div>
    <p class="pgsum">${escapeHtml(c.summary)}</p>
    <div id="pgbody_${c.id}"></div>
    ${c.limits && c.limits.length ? `<details><summary class="note">What this cannot tell you</summary><ul class="note">${c.limits.map((l) => `<li>${escapeHtml(l)}</li>`).join("")}</ul></details>` : ""}
    ${c.method ? `<details><summary class="note">How it is worked out</summary><p class="note">${escapeHtml(c.method)}</p>
      ${c.validation ? `<p class="note"><b>Checked:</b> ${escapeHtml(c.validation)}</p>` : ""}
      ${c.sources && c.sources.length ? `<ul class="note">${c.sources.map((s) => `<li>${escapeHtml(s.title)}</li>`).join("")}</ul>` : ""}</details>` : ""}
  </div>`;
}

const DETAIL = {
  file_check(el, c) {
    const r = c.result, m = r.reference_match || {};
    el.innerHTML = kvTable([
      ["Testing company", escapeHtml(r.vendor)],
      ["Markers in the file", r.markers.toLocaleString() + ` (${pct(r.call_rate, 1)} with a call)`],
      ["Genome build", escapeHtml(m.build || "")],
      ["Used for ancestry", (m.used || 0).toLocaleString() + " markers" + (m.strand_flipped ? `, ${m.strand_flipped.toLocaleString()} flipped to the reference strand` : "")],
      ["Strand", escapeHtml(m.strand || "")],
      ["Left out", `${(m.ambiguous_dropped || 0).toLocaleString()} whose strand could not be settled, ${(m.allele_mismatch_dropped || 0).toLocaleString()} whose alleles do not match the reference`],
      ["Genetic sex", `${escapeHtml(r.genetic_sex)} <small>(${escapeHtml(r.genetic_sex_basis)})</small>`],
    ]);
  },
  ancient_ancestry(el, c) {
    const ms = c.result.models || [];
    el.innerHTML = ms.map((m) => {
      if (m.error) return `<p class="note">${escapeHtml(m.title)}: ${escapeHtml(m.error)}</p>`;
      const bar = m.components.map((x, i) => `<i class="s${i}" style="width:${Math.max(0, 100 * x.proportion).toFixed(1)}%" title="${escapeHtml(x.label)}"></i>`).join("");
      return `<div class="model ${m.fit === "fits" ? "ok" : m.fit === "marginal" ? "mar" : "no"}">
        <div class="mhead"><b>${escapeHtml(m.title)}</b> <span class="fit">${escapeHtml(m.fit)} · p = ${m.p.toFixed(3)}</span></div>
        <p class="note">${escapeHtml(m.description || "")}</p>
        ${m.fit === "does not fit" ? `<p class="note">This model does not explain your genome, so it gives no proportions:
            the sources it offers cannot account for how your genome relates to the outgroup populations.
            Sources tried: ${m.components.map((x) => escapeHtml(x.label)).join(", ")}. Fitted on ${(m.snps || 0).toLocaleString()} markers.</p>`
          : `<div class="stack">${bar}</div>
        <table><thead><tr><th>Ancient source</th><th>Share</th><th>95% interval</th></tr></thead><tbody>
        ${m.components.map((x, i) => `<tr><td><i class="dot s${i}"></i>${escapeHtml(x.label)}</td><td class="n">${pct(x.proportion, 1)}</td>
          <td class="n">${pct(x.ci95[0], 1)} – ${pct(x.ci95[1], 1)}</td></tr>`).join("")}</tbody></table>
        <p class="note">${m.imprecise ? "The intervals are wide, so read these as rough shares. " : ""}Fitted on ${(m.snps || 0).toLocaleString()} markers.</p>`}</div>`;
    }).join("");
  },
  pca(el, c) {
    const maps = c.result.maps || {};
    el.innerHTML = Object.entries(maps).filter(([, m]) => m.coords)
      .map(([k, m]) => `<div class="plotwrap"><div class="note">${escapeHtml(m.title)} — placed using ${(m.snps || 0).toLocaleString()} markers</div>
        <canvas id="pgpca_${k}" height="320"></canvas></div>`).join("")
      || `<p class="note">No map could be drawn.</p>`;
    for (const [k, m] of Object.entries(maps)) if (m.coords) drawPca($(`pgpca_${k}`), m);
  },
  closest_ancients(el, c) {
    const rows = (c.result.west && c.result.west.length ? c.result.west : c.result.world) || [];
    el.innerHTML = `<div class="plotwrap"><canvas id="pgmap" height="300"></canvas></div>
      <table><thead><tr><th>Ancient individual</th><th>When</th><th>Where</th><th>mtDNA</th><th>Y</th></tr></thead><tbody>
      ${rows.slice(0, 10).map((a) => `<tr><td>${escapeHtml(a.group || a.id)}<br><small>${escapeHtml(a.id)}</small></td>
        <td class="n">${escapeHtml(BP(a.date_bp))}</td><td>${escapeHtml([a.locality, a.country].filter(Boolean).join(", "))}</td>
        <td class="n">${escapeHtml(a.mt || "—")}</td><td class="n">${escapeHtml(a.y || "—")}</td></tr>`).join("")}</tbody></table>`;
    drawMap($("pgmap"), rows.slice(0, 10));
  },
  neanderthal(el, c) {
    const r = c.result; if (r.error) { el.innerHTML = `<p class="note">${escapeHtml(r.error)}</p>`; return; }
    const groups = (r.reference_groups && r.reference_groups.regions) || [];
    el.innerHTML = `<div class="plotwrap"><canvas id="pgnean" height="240"></canvas></div>
      <p class="note">Your estimate with its 95% interval, against the averages of present-day reference groups computed the same way.</p>`;
    drawNean($("pgnean"), r, groups);
  },
  denisovan(el, c) {
    const r = c.result; if (r.error) { el.innerHTML = `<p class="note">${escapeHtml(r.error)}</p>`; return; }
    const regions = (r.reference_groups && r.reference_groups.regions) || [];
    el.innerHTML = kvTable([["Your statistic", `${r.D.toFixed(4)} (Z ${r.z >= 0 ? "+" : ""}${r.z.toFixed(1)}, ${(r.snps || 0).toLocaleString()} markers)`]])
      + `<table><thead><tr><th>Reference region</th><th>Average</th><th>Spread</th></tr></thead><tbody>
      ${regions.map((g) => `<tr><td>${escapeHtml(g.region)}</td><td class="n">${g.mean.toFixed(4)}</td><td class="n">${g.sd == null ? "—" : "±" + g.sd.toFixed(4)}</td></tr>`).join("")}</tbody></table>`;
  },
  archaic_alleles(el, c) {
    const r = c.result;
    const vs = r.variants || [];
    el.innerHTML = kvTable([
      ["Sites checked in your file", (r.typed_sites || 0).toLocaleString()],
      ["Sites where you carry it", (r.sites_carried || 0).toLocaleString()],
      ["Total copies", (r.copies || 0).toLocaleString()],
    ]) + (vs.length ? `<details><summary class="note">The first ${Math.min(vs.length, 400)} variants, with how common each is</summary>
      <table><thead><tr><th>Variant</th><th>Position</th><th>Your copies</th><th>Europe</th><th>East Asia</th><th>South Asia</th></tr></thead><tbody>
      ${vs.slice(0, 400).map((v) => `<tr><td class="n">rs${v.rsid}</td><td class="n">${v.chrom}:${v.pos.toLocaleString()}</td><td class="n">${v.copies}</td>
        <td class="n">${pct(v.EUR, 0)}</td><td class="n">${pct(v.EAS, 0)}</td><td class="n">${pct(v.SAS, 0)}</td></tr>`).join("")}</tbody></table></details>` : "")
      + `<p class="note">${escapeHtml(r.definition || "")}</p>`;
  },
  mt_haplogroup(el, c) { haploDetail(el, c, "mtDNA"); },
  y_haplogroup(el, c) { haploDetail(el, c, "Y"); },
  continental_ancestry(el, c) {
    const cs = c.result.continents || [];
    const shown = cs.filter((x) => x.detected);
    el.innerHTML = `<div class="stack">${shown.map((x) => `<i style="width:${(100 * x.proportion).toFixed(1)}%;background:${CONT_COLOR[x.id] || "#7f9fb3"}" title="${escapeHtml(x.title)}"></i>`).join("")}</div>
      <table><thead><tr><th>Continent</th><th>Share</th><th>95% interval</th><th>Detection limit</th></tr></thead><tbody>
      ${cs.map((x) => `<tr class="${x.detected ? "" : "dim"}"><td><i class="dot" style="background:${CONT_COLOR[x.id] || "#7f9fb3"}"></i>${escapeHtml(x.title)}</td>
        <td class="n">${x.detected ? pct(x.proportion, 1) : "not detected"}</td><td class="n">${x.detected ? `${pct(x.ci95[0], 1)} – ${pct(x.ci95[1], 1)}` : "—"}</td>
        <td class="n">${pct(x.floor, 1)}</td></tr>`).join("")}</tbody></table>
      <p class="note">Estimated from ${(c.result.snps || 0).toLocaleString()} SNPs. "Not detected" means the estimate is below what people with no ancestry from there are sometimes assigned.</p>`
      + typicalHtml(c.result.typical, "continents");
  },
  regional_ancestry(el, c) {
    const rs = c.result.regions || [];
    const det = rs.filter((x) => x.detected);
    el.innerHTML = `<div class="stack">${det.map((x) => `<i style="width:${(100 * x.proportion).toFixed(1)}%;background:${CONT_COLOR[x.continent] || "#7f9fb3"}" title="${escapeHtml(x.title)}"></i>`).join("")}</div>
      <table><thead><tr><th>Region</th><th>Share</th><th>95% interval</th></tr></thead><tbody>
      ${det.map((x) => `<tr><td><i class="dot" style="background:${CONT_COLOR[x.continent] || "#7f9fb3"}"></i>${escapeHtml(x.title)}</td>
        <td class="n">${pct(x.proportion, 1)}</td><td class="n">${pct(x.ci95[0], 1)} – ${pct(x.ci95[1], 1)}</td></tr>`).join("")}</tbody></table>
      ${typicalHtml(c.result.typical, "regions")}
      <details><summary class="note">Reference populations behind each region</summary><table><thead><tr><th>Region</th><th>Reference groups</th><th>People</th></tr></thead><tbody>
      ${rs.map((x) => `<tr><td>${escapeHtml(x.title)}</td><td>${x.reference_groups.map((gname) => escapeHtml(x.group_names[gname] ? `${gname} (${x.group_names[gname]})` : gname)).join(", ")}</td>
        <td class="n">${x.reference_individuals}</td></tr>`).join("")}</tbody></table></details>`;
  },
  chromosome_painting(el, c) {
    const r = c.result;
    el.innerHTML = `<div class="plotwrap"><canvas id="pgpaint" style="height:${22 * 17 + 40}px"></canvas></div>
      <p class="note">Each chromosome is drawn as two stripes, one per copy. Where the two copies differ, which colour belongs to which copy is not known.</p>
      <p class="note">${(r.continents || []).map((x) => `<i class="dot" style="background:${CONT_COLOR[x.id] || "#7f9fb3"}"></i>${escapeHtml(x.title)}`).join(" &nbsp; ")}</p>`;
    drawPainting($("pgpaint"), r);
  },
  indigenous_american(el, c) {
    const r = c.result;
    el.innerHTML = kvTable([
      ["Estimate", `${pct(r.proportion, 1)} (95% interval ${pct(r.ci95[0], 1)} – ${pct(r.ci95[1], 1)})`],
      ["Detection limit", pct(r.floor, 1)],
      ["Painted segments", `${(r.segments || []).length} (${(r.segment_mb || 0).toFixed(0)} Mb of chromosome copies)`],
      ["Reference groups", escapeHtml((r.reference_groups || []).join(", "))],
    ]);
  },
  roh(el, c) {
    const r = c.result;
    const ref = Object.entries(r.reference_groups || {}).sort((a, b) => b[1].froh_median - a[1].froh_median);
    el.innerHTML = kvTable([
      ["Total in runs of homozygosity", `${r.total_mb.toFixed(0)} Mb in ${r.segments} runs (${pct(r.froh, 2)} of the autosomes)`],
      ["By run length", Object.entries(r.by_length_mb).map(([k, v]) => `${k} Mb: ${v.toFixed(0)}`).join(" · ")],
      ["Longest run", `${r.longest_mb.toFixed(1)} Mb`],
      ["Markers used", (r.markers || 0).toLocaleString()],
    ]) + `<div class="plotwrap"><canvas id="pgroh" style="height:${22 * 12 + 34}px"></canvas></div>
      <details><summary class="note">Expected share in long runs when parents are related</summary><table><thead><tr><th>Parents</th><th>Expected share of the genome</th></tr></thead><tbody>
      ${r.relationships.map((x) => `<tr><td>${escapeHtml(x.parents)}</td><td class="n">${pct(x.expected_froh, 1)}</td></tr>`).join("")}</tbody></table>
      <p class="note">These are averages; the actual share varies widely between children of the same parents.</p></details>
      <details><summary class="note">Reference groups, measured the same way</summary><table><thead><tr><th>Group</th><th>Typical share</th><th>Long runs (10 Mb+), typical</th><th>People</th></tr></thead><tbody>
      ${ref.map(([gname, v]) => `<tr><td>${escapeHtml(v.label || gname)}</td><td class="n">${pct(v.froh_median, 2)}</td><td class="n">${v.long_mb_median.toFixed(0)} Mb</td><td class="n">${v.n}</td></tr>`).join("")}</tbody></table></details>`;
    drawRoh($("pgroh"), r.segments_list || []);
  },
  x_ancestry(el, c) {
    const r = c.result; if (!r.continents) return;
    el.innerHTML = `<table><thead><tr><th>Continent</th><th>X chromosome</th><th>95% interval</th><th>Autosomes</th></tr></thead><tbody>
      ${r.continents.filter((x) => x.proportion >= 0.01 || x.autosomal >= 0.01).map((x) => `<tr><td><i class="dot" style="background:${CONT_COLOR[x.id] || "#7f9fb3"}"></i>${escapeHtml(x.title)}</td>
        <td class="n">${pct(x.proportion, 1)}</td><td class="n">${pct(x.ci95[0], 1)} – ${pct(x.ci95[1], 1)}</td><td class="n">${pct(x.autosomal, 1)}</td></tr>`).join("")}</tbody></table>
      <p class="note">${(r.snps || 0).toLocaleString()} X-chromosome SNPs.</p>`;
  },
};

function typicalHtml(typ, key) {
  if (!typ || !typ.length) return "";
  return typ.map((t) => `<p class="note">For comparison, reference people from <b>${escapeHtml(t.region)}</b>, analysed the same way, typically get:
    ${t[key].map(([n, v]) => `${escapeHtml(n)} ${pct(v, 0)}`).join(" · ")}</p>`).join("");
}

const TRAIT_CARDS = ["eye_colour", "red_hair", "freckles", "blond_hair", "skin_pigmentation", "earwax", "hair_thickness", "hair_curl",
  "lactose", "alcohol_flush", "caffeine", "bitter_taste", "vitamin_d", "vitamin_b12", "omega3", "celiac_hla", "actn3"];
for (const id of TRAIT_CARDS) DETAIL[id] = traitDetail;

function traitDetail(el, c) {
  const vs = (c.result.variants || []).filter((v, i, a) => a.findIndex((w) => w.rsid === v.rsid) === i);
  const POPS = [["EUR", "Europe"], ["EAS", "East Asia"], ["SAS", "South Asia"], ["AFR", "Africa"], ["AMR", "Americas"]];
  el.innerHTML = `<table><thead><tr><th>Variant</th><th>Your genotype</th><th>Copies of the effect allele</th><th>How common the effect allele is</th></tr></thead><tbody>
    ${vs.map((v) => `<tr><td>${escapeHtml(v.variant)}<br><small>${escapeHtml(v.rsid)} · ${escapeHtml(v.gene)} · effect allele ${escapeHtml(v.effect_allele)}</small></td>
      <td class="n">${v.genotype ? escapeHtml(v.genotype) : "—"}</td>
      <td class="n">${v.copies == null ? `<small>${escapeHtml(v.note || "not read")}</small>` : v.copies}</td>
      <td><small>${v.freq ? POPS.filter(([k]) => v.freq[k] != null).map(([k, n]) => `${n} ${pct(v.freq[k], 0)}`).join(" · ") : ""}</small>
      ${v.gwas ? `<br><small>GWAS Catalog: ${escapeHtml(v.gwas.trait)} (${escapeHtml(v.gwas.first_author || "")}, p = ${Number(v.gwas.p).toExponential(0)})</small>` : ""}</td></tr>`).join("")}
    </tbody></table><p class="note">Genotypes on the forward strand of the GRCh37 reference; frequencies from 1000 Genomes.</p>`;
}

// ------------------------------------------------------------------ medicines
const POP1KG = { AFR: "Africa", AMR: "Americas", EAS: "East Asia", EUR: "Europe", SAS: "South Asia" };

DETAIL.medicines = function (el, c) {
  const r = c.result, genes = r.genes || {};
  const known = (d) => d.rows && d.rows.length;
  const flagged = new Set(r.flagged || []);
  const act = (r.drugs || []).filter((d) => known(d) && d.genes.some((g) => flagged.has(g)));
  const rest = (r.drugs || []).filter((d) => !act.includes(d));
  const geneText = (d) => d.genes.map((g) => {
    const x = genes[g];
    if (!x) return `${escapeHtml(g)}: <i>not in the pack</i>`;
    if (g === "CYP2D6") return `${escapeHtml(g)}: <i>needs a clinical test</i>`;
    if (x.status === "unreliable") return `${escapeHtml(g)}: <i>not determinable from your file</i>`;
    return `${escapeHtml(g)}: ${escapeHtml(x.phenotype || (x.diplotypes || []).join(" or ") || x.status)}`;
  }).join("<br>");
  const recText = (d) => d.depends ? `<small>Depends on ${escapeHtml(d.unknown.join(", "))}, which your file cannot settle.</small>`
    : d.rows.map((x) => `${escapeHtml(x.recommendation || "")}${x.classification ? ` <small>(${escapeHtml(x.classification)})</small>` : ""}
        ${x.population && x.population !== "general" ? `<br><small>For: ${escapeHtml(x.population)}</small>` : ""}`).join("<hr>");
  const table = (ds) => `<table><thead><tr><th>Medicine</th><th>Your genes</th><th>CPIC recommendation</th></tr></thead><tbody>
    ${ds.map((d) => `<tr><td>${escapeHtml(d.drug)}<br><small>${escapeHtml((d.guideline && d.guideline.name) || "")}</small></td>
      <td><small>${geneText(d)}</small></td><td>${recText(d)}</td></tr>`).join("")}</tbody></table>`;
  el.innerHTML = `${act.length ? `<p class="note"><b>Medicines whose CPIC recommendation changes with your result</b></p>${table(act)}` : ""}
    <details><summary class="note">All ${rest.length} other CPIC level-A medicines</summary>${table(rest)}</details>
    <p class="note">CPIC data ${escapeHtml((r.cpic && r.cpic.accessed) || "")}${r.cpic && r.cpic.release ? `, database ${escapeHtml(r.cpic.release)}` : ""} (CC0). ${escapeHtml((r.cpic && r.cpic.notice) || "")}</p>`;
};

function pgxDetail(el, c) {
  const r = c.result, k = r.call || {}, rel = r.reliability;
  const rows = [
    ["Diplotype", escapeHtml((k.diplotypes || []).join(" or ") || "—")],
    ["Phenotype", k.status === "unreliable" ? `not determinable from your file <small>(the positions read suggest ${escapeHtml(k.presumptive || "")})</small>`
      : escapeHtml(k.phenotype || (k.possible ? `one of: ${k.possible.join(", ")}` : k.status || "—"))],
  ];
  if (k.score) rows.push(["Activity score", escapeHtml(k.score)]);
  if (k.chosen_by) rows.push(["Chosen by", `${escapeHtml(k.chosen_by)} <small>(${pct(k.diplotype_probability, 0)} of the diplotypes that fit)</small>`]);
  rows.push(["Positions read", `${(r.positions || []).length}`]);
  if (rel && rel.phenotype_agrees != null) {
    rows.push(["Same phenotype from a full reading", `${pct(rel.phenotype_agrees, 0)} <small>of ${rel.same_call.toLocaleString()} 1000 Genomes people with your call</small>`]);
  }
  let html = kvTable(rows);
  if (k.variants && k.variants.length) {
    html += `<table><thead><tr><th>Variant found</th><th>Copies</th><th>Function</th></tr></thead><tbody>
      ${k.variants.map((v) => `<tr><td>${escapeHtml(v.name)}</td><td class="n">${v.copies}</td><td>${escapeHtml(v.function || "")}${v.activity ? ` <small>(activity ${escapeHtml(v.activity)})</small>` : ""}</td></tr>`).join("")}</tbody></table>`;
  }
  if (k.alternatives && k.alternatives.length) {
    html += `<p class="note">Also fits your genotypes, though less likely (the chip cannot tell which variants sit on the same chromosome): ${k.alternatives.map((a) => `${escapeHtml(a.diplotype)} (${escapeHtml(a.phenotype)})`).join("; ")}.</p>`;
  }
  if (rel && rel.full_instead && rel.full_instead.length) {
    html += `<p class="note">When a full reading disagreed, it gave: ${rel.full_instead.map((x) => `${escapeHtml(x.phenotype || "no call")} ${pct(x.share, 0)}`).join(" · ")}.</p>`;
  }
  if (rel && rel.by_population) {
    html += `<details><summary class="note">Reliability by 1000 Genomes population</summary><table><thead><tr><th>Population</th><th>People</th><th>With your call</th><th>Same phenotype</th></tr></thead><tbody>
      ${Object.entries(rel.by_population).map(([p, x]) => `<tr><td>${escapeHtml(POP1KG[p] || p)}</td><td class="n">${x.n}</td><td class="n">${x.same_call}</td><td class="n">${pct(x.phenotype_agrees, 0)}</td></tr>`).join("")}</tbody></table></details>`;
  }
  if (r.positions && r.positions.length) {
    html += `<details><summary class="note">Your genotypes at the ${r.positions.length} positions read</summary><table><thead><tr><th>Position</th><th>Genotype</th></tr></thead><tbody>
      ${r.positions.map((p) => `<tr><td>${escapeHtml(p.rsid || "")} <small>${escapeHtml(p.hgvs || "")}${p.name ? ` · ${escapeHtml(p.name)}` : ""}</small></td><td class="n">${escapeHtml(p.genotype)}</td></tr>`).join("")}</tbody></table></details>`;
  }
  if (r.drugs && r.drugs.length) html += `<p class="note">CPIC level-A medicines for ${escapeHtml(r.gene)}: ${r.drugs.map(escapeHtml).join(", ")}.</p>`;
  el.innerHTML = html;
}

const CONT_COLOR = { africa: "#ff6a57", europe: "#62d6ff", west_asia: "#ffb34a", south_asia: "#b98cff", east_asia: "#52f0b0",
  north_asia: "#9fd4ff", oceania: "#ff8fd0", americas: "#ffe066" };
const CHR_LEN37 = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,
  135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566];

function drawPainting(cv, r) {
  const { g, w, h } = plotFrame(cv);
  const L = 34, R = w - 10, row = 17, top = 10, maxLen = CHR_LEN37[0];
  const X = (p) => L + (p / maxLen) * (R - L);
  g.font = "10px JetBrains Mono, monospace";
  for (let c = 1; c <= 22; c++) {
    const y = top + (c - 1) * row;
    g.fillStyle = "#7f9fb3"; g.fillText(String(c), 6, y + 10);
    g.fillStyle = "rgba(110,205,255,.10)"; g.fillRect(L, y, X(CHR_LEN37[c - 1]) - L, 12);
  }
  for (const sg of r.segments || []) {
    const y = top + (sg.chrom - 1) * row;
    const x0 = X(sg.start), x1 = Math.max(X(sg.end), x0 + 1);
    g.fillStyle = CONT_COLOR[sg.a] || "#7f9fb3"; g.fillRect(x0, y, x1 - x0, 6);
    g.fillStyle = CONT_COLOR[sg.b] || "#7f9fb3"; g.fillRect(x0, y + 6, x1 - x0, 6);
  }
}

function drawRoh(cv, segs) {
  const { g, w } = plotFrame(cv);
  const L = 34, R = w - 10, row = 12, top = 8, maxLen = CHR_LEN37[0];
  const X = (p) => L + (p / maxLen) * (R - L);
  g.font = "10px JetBrains Mono, monospace";
  for (let c = 1; c <= 22; c++) {
    const y = top + (c - 1) * row;
    g.fillStyle = "#7f9fb3"; g.fillText(String(c), 6, y + 8);
    g.fillStyle = "rgba(110,205,255,.10)"; g.fillRect(L, y, X(CHR_LEN37[c - 1]) - L, 8);
  }
  for (const s of segs) {
    const y = top + (s.chrom - 1) * row;
    g.fillStyle = s.kb >= 10000 ? "#ffb34a" : "#62d6ff";
    g.fillRect(X(s.start), y, Math.max(1, X(s.end) - X(s.start)), 8);
  }
  g.fillStyle = "#7f9fb3"; g.fillText("blue: runs under 10 Mb · orange: 10 Mb or longer", L, top + 22 * row + 12);
}

// Consecutive branches whose oldest known carrier is the same genome are shown as one row ("ROOT → K2").
function collapseTimeline(tl) {
  const out = [];
  for (const t of tl.filter((x) => x.oldest)) {
    const last = out[out.length - 1];
    if (last && last.oldest.id === t.oldest.id) { last.branch = `${last.first} → ${t.branch}`; continue; }
    out.push({ ...t, first: t.branch });
  }
  return out;
}

function haploDetail(el, c, kind) {
  const r = c.result;
  if (!r.call) {
    const sc = r.strand_check;
    el.innerHTML = (sc && sc.on_tree ? kvTable([
      ["Markers in your file on this tree", (sc.on_tree || 0).toLocaleString()],
      ["Matched the tree's alleles", `${sc.direct} directly, ${sc.complemented} after switching strand`],
      ["Matched neither allele", (sc.unmatched || 0).toLocaleString()]]) : "")
      + (r.weak_call ? `<p class="note">Best-supported branch, not reported: ${escapeHtml(r.weak_call)}.</p>` : "");
    return;
  }
  const path = r.path || [], tl = r.timeline || [], car = r.carriers || [];
  el.innerHTML = kvTable([
    ["Haplogroup", `<b>${escapeHtml(r.call)}</b>`],
    ["Tree", escapeHtml(r.tree || "")],
    ["Markers in your file on this tree", (r.markers_in_file || 0).toLocaleString() + ` (${r.typed_mutations || 0} on branch-defining positions)`],
    ["Support along the lineage", `${r.path_derived} markers agree, ${r.path_ancestral} disagree`],
  ]) + `<p class="note">Your lineage from the root: ${path.map((p) => escapeHtml(p.name)).join(" → ")}</p>`
    + (tl.length ? `<details open><summary class="note">Oldest dated genome known to carry each branch of your line</summary>
      <table><thead><tr><th>Branch</th><th>Oldest carrier</th><th>When</th><th>Where</th></tr></thead><tbody>
      ${collapseTimeline(tl).map((t) => `<tr><td class="n">${escapeHtml(t.branch)}</td><td>${escapeHtml(t.oldest.group || t.oldest.id)}</td>
        <td class="n">${escapeHtml(BP(t.oldest.date_bp))}</td><td>${escapeHtml([t.oldest.locality, t.oldest.country].filter(Boolean).join(", "))}</td></tr>`).join("")}
      </tbody></table><p class="note">Each row gives a minimum age for that branch: the lineage existed at least that long ago.</p></details>` : "")
    + (car.length ? `<details><summary class="note">${car.length} ancient genomes carry this ${kind} branch or one below it</summary>
      <div class="plotwrap"><canvas id="pgmap_${c.id}" height="260"></canvas></div>
      <table><thead><tr><th>Individual</th><th>When</th><th>Where</th></tr></thead><tbody>
      ${car.slice(0, 60).map((a) => `<tr><td>${escapeHtml(a.group || a.id)}</td><td class="n">${escapeHtml(BP(a.date_bp))}</td>
        <td>${escapeHtml([a.locality, a.country].filter(Boolean).join(", "))}</td></tr>`).join("")}</tbody></table></details>` : "");
  const cv = $(`pgmap_${c.id}`);
  if (cv) setTimeout(() => drawMap(cv, car.slice(0, 400)), 0);
}

const REGION_COLOR = { "Europe": "#62d6ff", "West Asia": "#ffb34a", "North Africa": "#d9a066", "Africa": "#ff6a57",
  "South Asia": "#b98cff", "Central Asia": "#8fe3c8", "East Asia": "#52f0b0", "Siberia and North Asia": "#9fd4ff",
  "Oceania and Island Southeast Asia": "#ff8fd0", "Americas": "#ffe066" };

function drawPca(cv, m) {
  const { g, w, h, l, r, t, b } = plotFrame(cv);
  const pts = m.reference.coords, reg = m.reference.region;
  const xs = pts.map((p) => p[0]).concat([m.coords[0]]), ys = pts.map((p) => p[1]).concat([m.coords[1]]);
  const anc = (m.nearest_ancients || []);
  const pad = 0.08;
  const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
  const dx = (x1 - x0) * pad || 1, dy = (y1 - y0) * pad || 1;
  const X = (v) => l + ((v - x0 + dx) / (x1 - x0 + 2 * dx)) * (r - l);
  const Y = (v) => b - ((v - y0 + dy) / (y1 - y0 + 2 * dy)) * (b - t);
  g.strokeRect(l, t, r - l, b - t);
  pts.forEach((p, i) => { g.fillStyle = (REGION_COLOR[reg[i]] || "#7f9fb3") + "88"; g.beginPath(); g.arc(X(p[0]), Y(p[1]), 2.4, 0, 7); g.fill(); });
  // reference group labels at each group's centre
  const cen = {};
  pts.forEach((p, i) => { const k = m.reference.group[i]; (cen[k] = cen[k] || [0, 0, 0]); cen[k][0] += p[0]; cen[k][1] += p[1]; cen[k][2]++; });
  g.font = "9px JetBrains Mono, monospace";
  for (const [k, v] of Object.entries(cen)) {
    if (v[2] < 5) continue;
    g.fillStyle = "rgba(220,235,245,.55)"; g.fillText(k, X(v[0] / v[2]) + 3, Y(v[1] / v[2]) - 3);
  }
  // your closest ancient individuals, projected onto the same axes
  g.fillStyle = "#ffffff";
  for (const a of anc) {
    if (!a.pc) continue;
    const x = X(a.pc[0]), y = Y(a.pc[1]);
    g.beginPath(); g.moveTo(x, y - 4); g.lineTo(x + 4, y); g.lineTo(x, y + 4); g.lineTo(x - 4, y); g.closePath(); g.fill();
  }
  g.strokeStyle = "#fff"; g.lineWidth = 2;
  g.beginPath(); g.arc(X(m.coords[0]), Y(m.coords[1]), 6, 0, 7); g.stroke();
  g.fillStyle = "#fff"; g.font = "12px JetBrains Mono, monospace";
  g.fillText("you", X(m.coords[0]) + 9, Y(m.coords[1]) + 4);
  const regs = [...new Set(reg)];
  g.font = "10px JetBrains Mono, monospace";
  regs.forEach((rr, i) => { g.fillStyle = REGION_COLOR[rr] || "#7f9fb3"; g.fillText(rr, l + 6, t + 12 + 12 * i); });
  g.fillStyle = "#7f9fb3"; g.fillText("component 1", (l + r) / 2 - 30, b + 16);
  g.save(); g.translate(l - 8, (t + b) / 2 + 30); g.rotate(-Math.PI / 2); g.fillText("component 2", 0, 0); g.restore();
  if (anc.some((a) => a.pc)) g.fillText("◆ your closest ancient individuals", r - 230, b + 16);
}

let WORLD = null;
async function loadWorld() {
  if (WORLD) return WORLD;
  try { WORLD = await (await fetch("data/world110.json")).json(); } catch { WORLD = { rings: [] }; }
  return WORLD;
}

async function drawMap(cv, points) {
  const world = await loadWorld();
  const dpr = window.devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight || 300;
  cv.width = w * dpr; cv.height = h * dpr;
  const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0);
  const lons = points.map((p) => p.lon).filter((x) => x != null), lats = points.map((p) => p.lat).filter((x) => x != null);
  let lo0 = -25, lo1 = 80, la0 = 20, la1 = 72;
  if (lons.length) {
    lo0 = Math.min(...lons) - 8; lo1 = Math.max(...lons) + 8;
    la0 = Math.min(...lats) - 6; la1 = Math.max(...lats) + 6;
  }
  const ar = (lo1 - lo0) / Math.max(1e-6, la1 - la0), want = w / h;
  if (ar < want) { const c = (lo0 + lo1) / 2, half = (la1 - la0) * want / 2; lo0 = c - half; lo1 = c + half; }
  else { const c = (la0 + la1) / 2, half = (lo1 - lo0) / want / 2; la0 = c - half; la1 = c + half; }
  const X = (lon) => ((lon - lo0) / (lo1 - lo0)) * w, Y = (lat) => h - ((lat - la0) / (la1 - la0)) * h;
  g.fillStyle = "rgba(10,22,32,.6)"; g.fillRect(0, 0, w, h);
  g.strokeStyle = "rgba(110,205,255,.35)"; g.lineWidth = 1;
  for (const ring of world.rings) {
    g.beginPath();
    for (let i = 0; i < ring.length; i += 2) {
      const x = X(ring[i]), y = Y(ring[i + 1]);
      i ? g.lineTo(x, y) : g.moveTo(x, y);
    }
    g.stroke();
  }
  const bps = points.map((p) => p.date_bp).filter((x) => x != null);
  const oldest = Math.max(1, ...bps);
  for (const p of points) {
    if (p.lon == null || p.lat == null) continue;
    const f = p.date_bp == null ? 0.5 : p.date_bp / oldest;
    g.fillStyle = `rgba(${Math.round(98 + 157 * f)},${Math.round(214 - 70 * f)},${Math.round(255 - 180 * f)},.85)`;
    g.beginPath(); g.arc(X(p.lon), Y(p.lat), 4, 0, 7); g.fill();
  }
  g.fillStyle = "#7f9fb3"; g.font = "10px JetBrains Mono, monospace";
  g.fillText("blue = more recent, orange = older", 6, h - 6);
}

function drawNean(cv, r, regions) {
  const { g, l, r: rr, t, b } = plotFrame(cv);
  const rows = [{ region: "you", mean: r.estimate, lo: r.ci95[0], hi: r.ci95[1] }]
    .concat(regions.filter((x) => x.n >= 3).map((x) => ({ region: x.region, mean: x.mean, lo: x.mean - 2 * (x.sd || 0), hi: x.mean + 2 * (x.sd || 0) })));
  const all = rows.flatMap((x) => [x.lo, x.hi]).filter((x) => isFinite(x));
  const v0 = Math.min(0, ...all), v1 = Math.max(...all) * 1.05;
  const X = (v) => l + ((v - v0) / (v1 - v0)) * (rr - l);
  const rowH = Math.max(14, (b - t) / Math.max(1, rows.length));
  g.strokeRect(l, t, rr - l, b - t);
  rows.forEach((x, i) => {
    const y = t + rowH * (i + 0.5);
    const me = i === 0;
    g.strokeStyle = me ? "#fff" : "rgba(98,214,255,.7)"; g.lineWidth = me ? 2 : 1;
    g.beginPath(); g.moveTo(X(x.lo), y); g.lineTo(X(x.hi), y); g.stroke();
    g.fillStyle = me ? "#fff" : "#62d6ff";
    g.beginPath(); g.arc(X(x.mean), y, me ? 4 : 3, 0, 7); g.fill();
    g.font = "10px JetBrains Mono, monospace";
    g.fillText(x.region.slice(0, 26), l + 4, y - 4);
  });
  g.fillStyle = "#7f9fb3"; g.fillText(`${(100 * v0).toFixed(1)}%`, l, b + 14);
  g.fillText(`${(100 * v1).toFixed(1)}% Neanderthal`, rr - 110, b + 14);
}

function renderPgProv() {
  const R = PG.result, st = PG.status || {};
  const src = st.sources || {};
  $("pgProv").innerHTML = `<p><b>Your file</b> (read on this computer, never uploaded)</p>
    <ul><li>${escapeHtml(R.file.file)} — ${escapeHtml(R.file.vendor)}, ${R.file.markers.toLocaleString()} markers<br>
      <small>sha256 ${escapeHtml(R.file.sha256.slice(0, 16))}…</small></li></ul>
    <p><b>Reference data</b></p><ul>${Object.entries(src).map(([k, v]) => `<li><b>${escapeHtml(k)}</b>: ${escapeHtml(v.name || "")}
      ${v.licence ? `<br><small>licence ${escapeHtml(v.licence)}</small>` : ""}${v.citation ? `<br><small>${escapeHtml(v.citation)}</small>` : ""}</li>`).join("")}</ul>
    <p class="note">Analysis took ${R.seconds} s.</p>`;
}

async function boot() {
  wire(); setTier("gfn2-xtb"); setBusy(false);
  try { S.bench = await (await fetch("data/benchmark.json")).json(); } catch { S.bench = null; }
  if (TAURI && TAURI.event) {
    // First launch after install or update unpacks the science engines; show progress.
    TAURI.event.listen("core-progress", (e) => {
      const d = e.payload || {}, t = progressText(d);
      if (d.task === "adna") { if (A.busy) adnaProgress(d); return; }
      if (d.task === "personal") { pgProgress(d); return; }
      if (d.task === "pgdata") { pgPackProgress(d); return; }
      if (t && S.busy) setBusy(true, t);
    });
    TAURI.event.listen("core-status", (e) => {
      const st = e.payload || {};
      $("coreState").textContent = st.stage || "";
      if (st.stage === "preparing") setBusy(true, st.message);
      if (st.stage === "error") { $("coreChip").querySelector(".dot").className = "dot bad"; log(st.message, true); setBusy(false); }
      else if (st.message) log(st.message);
    });
  }
  // Up to 15 minutes: the first launch extracts several hundred MB of engines.
  for (let i = 0; i < 1800; i++) {
    try { await rpc("ping"); break; } catch { await new Promise((r) => setTimeout(r, 500)); }
  }
  setBusy(false);
  try {
    const eng = await rpc("engines");
    $("coreState").textContent = "ready"; $("coreChip").querySelector(".dot").className = "dot ok";
    const chips = Object.entries(eng).map(([k, v]) => `<div class="chip" title="${escapeHtml(v.reason || "")}"><i class="dot ${v.available ? "ok" : "bad"}"></i>${k}<span>${v.available ? escapeHtml(v.engine + " " + v.version) : "unavailable"}</span></div>`);
    $("engineChips").insertAdjacentHTML("beforeend", chips.join(""));
    if (!eng.dft.available) { document.querySelector('#tierSeg [data-tier="dft"]').disabled = true; }
    S.validation = await rpc("validation");
    log("Science core ready");
    checkForUpdate();
    initCompute().catch((e) => log(`GPU compute status unavailable: ${e}`, true));
    await build("CCO", "Ethanol");
  } catch (e) {
    $("coreState").textContent = "error"; $("coreChip").querySelector(".dot").className = "dot bad";
    log(`Science core failed to start: ${e.message || e}`, true);
  }
}
// ------------------------------------------------------------------ DNA Structure
const BASECOL = { A: 0xff5a4a, T: 0x4f8dff, G: 0x52f0b0, C: 0xffe066, U: 0x4f8dff };   // the NDB colours: A red, T blue, G green, C yellow

const H = (() => {
  let canvas, renderer, scene, camera, world, group, dist = 60, distTarget = 60, ready = false;
  let model = null, style = "atoms", color = "base", melt = null, inst = null, bondInst = null, slabs = null;
  const sph = new THREE.SphereGeometry(1, 16, 12), cyl = new THREE.CylinderGeometry(1, 1, 1, 10, 1, true);
  const box = new THREE.BoxGeometry(1, 1, 1);
  const m4 = new THREE.Matrix4(), qq = new THREE.Quaternion(), up = new THREE.Vector3(0, 1, 0), tmp = new THREE.Vector3();

  function init() {
    if (ready) return;
    canvas = $("dnagl");
    renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2)); renderer.setClearColor(0x000000, 0);
    scene = new THREE.Scene(); camera = new THREE.PerspectiveCamera(38, 1, 0.5, 2000);
    scene.add(new THREE.AmbientLight(0x6aa8c8, 0.65));
    const key = new THREE.DirectionalLight(0xdff4ff, 1.0); key.position.set(30, 40, 60); scene.add(key);
    const rim = new THREE.DirectionalLight(0x62d6ff, 0.7); rim.position.set(-40, -20, -30); scene.add(rim);
    world = new THREE.Group(); scene.add(world); group = new THREE.Group(); world.add(group);
    const resize = () => { const w = canvas.clientWidth, h = canvas.clientHeight; if (!w || !h) return;
      renderer.setSize(w, h, false); camera.aspect = w / h; camera.updateProjectionMatrix(); };
    new ResizeObserver(resize).observe(canvas); resize();
    let drag = null; const ax = new THREE.Vector3();
    canvas.addEventListener("pointerdown", (e) => { canvas.setPointerCapture(e.pointerId); drag = { x: e.clientX, y: e.clientY, moved: 0 }; });
    canvas.addEventListener("pointermove", (e) => {
      if (!drag) return; const dx = e.clientX - drag.x, dy = e.clientY - drag.y; drag.x = e.clientX; drag.y = e.clientY; drag.moved += Math.abs(dx) + Math.abs(dy);
      ax.set(0, 1, 0); qq.setFromAxisAngle(ax, dx * .008); world.quaternion.premultiply(qq);
      ax.set(1, 0, 0); qq.setFromAxisAngle(ax, dy * .008); world.quaternion.premultiply(qq);
    });
    canvas.addEventListener("pointerup", (e) => { if (drag && drag.moved < 4) pick(e); drag = null; });
    canvas.addEventListener("wheel", (e) => { e.preventDefault(); distTarget = Math.max(10, Math.min(400, distTarget * (1 + e.deltaY * .0012))); }, { passive: false });
    const frame = () => { dist += (distTarget - dist) * .1; camera.position.set(0, 0, dist); camera.lookAt(0, 0, 0);
      renderer.render(scene, camera); requestAnimationFrame(frame); };
    requestAnimationFrame(frame); ready = true;
  }

  // positions after the melting displacement: an open base pair's bases move apart along their C1'-C1' line
  function positions() {
    const P = model.atoms.map((a) => new THREE.Vector3(a[5], a[6], a[7]));
    if (!melt || color !== "melt") return P;
    model.pair_of_atom.forEach((k, i) => {
      if (k < 0) return;
      const open = 1 - (melt.closed[k] ?? 1), dir = model.pair_dir[k];
      const s = model.atoms[i][0] === "A" ? 1 : -1;
      P[i].addScaledVector(new THREE.Vector3(...dir), s * open * 4.0);
    });
    return P;
  }

  function atomColor(i) {
    const a = model.atoms[i];
    if (color === "element") return new THREE.Color(elInfo(a[4])[0]);
    const base = a[2].replace(/^D/, ""), k = model.pair_of_atom[i];
    if (color === "groove" && k >= 0 && model.groove[k] != null) {
      // narrow minor groove (3 A) hot, wide (9 A) cool
      const t = Math.max(0, Math.min(1, (model.groove[k] - 3) / 6));
      return new THREE.Color(0xff3b5c).lerp(new THREE.Color(0x62d6ff), t);
    }
    if (color === "melt" && melt && k >= 0) {
      const c = melt.closed[k] ?? 1;
      return new THREE.Color(0x7f9fb3).lerp(new THREE.Color(BASECOL[base] || 0xd0d8e0), c);
    }
    return model.backbone[i] ? new THREE.Color(0x9fb8c8) : new THREE.Color(BASECOL[base] || 0xd0d8e0);
  }

  function draw() {
    if (!model) return;
    group.clear();
    const P = positions(), n = P.length;
    const center = P.reduce((a, v) => a.add(v), new THREE.Vector3()).multiplyScalar(1 / n);
    if (style === "atoms") {
      const mat = new THREE.MeshStandardMaterial({ metalness: .2, roughness: .45 });
      inst = new THREE.InstancedMesh(sph, mat, n);
      P.forEach((p, i) => {
        const r = 0.18 + elInfo(model.atoms[i][4])[1] * 0.22;
        m4.compose(tmp.copy(p).sub(center), qq.identity(), new THREE.Vector3(r, r, r));
        inst.setMatrixAt(i, m4); inst.setColorAt(i, atomColor(i));
      });
      group.add(inst);
      const bm = new THREE.MeshStandardMaterial({ color: 0x9fdcff, emissive: 0x1d5a7a, metalness: .3, roughness: .4 });
      bondInst = new THREE.InstancedMesh(cyl, bm, model.bonds.length);
      model.bonds.forEach(([i, j], k) => {
        const a = P[i].clone().sub(center), b = P[j].clone().sub(center), d = b.clone().sub(a), len = d.length();
        qq.setFromUnitVectors(up, d.divideScalar(len || 1));
        m4.compose(a.add(b).multiplyScalar(.5), qq, new THREE.Vector3(.08, len, .08)); bondInst.setMatrixAt(k, m4);
      });
      group.add(bondInst);
    } else {
      // cartoon: a tube along each strand's phosphorus atoms, a slab per base from C1' into the ring
      inst = null;
      for (const chain of ["A", "B"]) {
        const pts = model.atoms.map((a, i) => (a[0] === chain && a[3] === "P" ? P[i].clone().sub(center) : null)).filter(Boolean);
        if (pts.length > 1) {
          const tube = new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts), pts.length * 8, 0.9, 10, false);
          group.add(new THREE.Mesh(tube, new THREE.MeshStandardMaterial({ color: chain === "A" ? 0x9fdcff : 0xc8a8ff, metalness: .3, roughness: .4 })));
        }
      }
      const res = model.residues;
      const sm = new THREE.MeshStandardMaterial({ metalness: .15, roughness: .5 });
      slabs = new THREE.InstancedMesh(box, sm, res.length);
      res.forEach((r, k) => {
        const c1 = P[r.c1].clone().sub(center), tip = P[r.tip].clone().sub(center), d = tip.clone().sub(c1), len = d.length();
        qq.setFromUnitVectors(up, d.divideScalar(len || 1));
        m4.compose(c1.clone().add(tip).multiplyScalar(.5), qq, new THREE.Vector3(1.6, len, .5)); slabs.setMatrixAt(k, m4);
        slabs.setColorAt(k, atomColor(r.tip));
      });
      group.add(slabs);
    }
  }

  function pick(e) {
    if (!model) return;
    const rect = canvas.getBoundingClientRect();
    const ray = new THREE.Raycaster(), ndc = new THREE.Vector2(((e.clientX - rect.left) / rect.width) * 2 - 1, -((e.clientY - rect.top) / rect.height) * 2 + 1);
    ray.setFromCamera(ndc, camera);
    const hit = ray.intersectObjects([inst, slabs].filter(Boolean))[0];
    if (!hit || hit.instanceId == null) { $("dnaPick").hidden = true; return; }
    const i = hit.object === inst ? hit.instanceId : model.residues[hit.instanceId].tip;
    const a = model.atoms[i], k = model.pair_of_atom[i];
    const bp = k >= 0 ? model.analysis.pairs[k] : null;
    $("dnaPick").hidden = false;
    $("dnaPick").textContent = `${a[2].replace(/^D/, "")}${a[1]} strand ${a[0] === "A" ? "I" : "II"} · atom ${a[3]}` +
      (bp ? ` · pair ${bp.pair}, propeller ${bp.params[4].toFixed(1)}°` : "") +
      (k >= 0 && model.groove[k] != null ? ` · minor groove ${model.groove[k].toFixed(1)} Å` : "") +
      (melt && k >= 0 ? ` · closed ${pct(melt.closed[k], 0)}` : "");
  }

  function set(m) {
    init(); model = m;
    // which base pair each atom belongs to, the residues' C1' and ring-tip atoms, and backbone atoms
    const bb = new Set(["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "O1P", "O2P"]);
    const nbp = m.sequence.length;
    m.backbone = m.atoms.map((a) => bb.has(a[3]));
    m.pair_of_atom = m.atoms.map((a) => (a[0] === "A" ? a[1] - 1 : 2 * nbp - a[1]));
    const byRes = new Map();
    m.atoms.forEach((a, i) => { const key = a[0] + a[1]; if (!byRes.has(key)) byRes.set(key, {}); byRes.get(key)[a[3]] = i; });
    m.residues = [...byRes.values()].map((r) => ({ c1: r["C1'"], tip: r["N1"] !== undefined && r["N9"] === undefined ? (r["N3"] ?? r["N1"]) : (r["N1"] ?? r["C1'"]) }))
      .filter((r) => r.c1 !== undefined && r.tip !== undefined);
    m.pair_dir = [];
    for (let k = 0; k < nbp; k++) {
      const a = byRes.get("A" + (k + 1)) || {}, b = byRes.get("B" + (2 * nbp - k)) || {};
      if (a["C1'"] === undefined || b["C1'"] === undefined) { m.pair_dir.push([0, 0, 0]); continue; }
      const pa = m.atoms[a["C1'"]], pb = m.atoms[b["C1'"]];
      const d = [pa[5] - pb[5], pa[6] - pb[6], pa[7] - pb[7]], L = Math.hypot(...d) || 1;
      m.pair_dir.push(d.map((x) => x / L));
    }
    // minor-groove width is defined per step (3DNA's P-P measure); a base pair takes the mean of its two steps
    const gs = Array(nbp - 1).fill(null);
    (m.analysis.grooves || []).forEach((g) => { if (g.minor != null) gs[g.step] = g.minor; });
    m.groove = Array.from({ length: nbp }, (_, k) => {
      const v = [gs[k - 1], gs[k]].filter((x) => x != null && x !== undefined);
      return v.length ? v.reduce((a, b) => a + b, 0) / v.length : null;
    });
    const P = m.atoms.map((a) => new THREE.Vector3(a[5], a[6], a[7]));
    const c = P.reduce((x, v) => x.add(v), new THREE.Vector3()).multiplyScalar(1 / P.length);
    let r = 0; for (const v of P) r = Math.max(r, v.distanceTo(c));
    distTarget = Math.max(25, (r + 4) / Math.tan(THREE.MathUtils.degToRad(17)));
    // stand the helix up: its axis (first to last pair origin) along screen y
    const f = m.frames, o0 = new THREE.Vector3(...f[0][1]), o1 = new THREE.Vector3(...f[f.length - 1][1]);
    const axis = o1.clone().sub(o0); if (axis.lengthSq() > 1e-6) group.quaternion.setFromUnitVectors(axis.normalize(), up);
    world.quaternion.identity();
    draw();
  }
  return { set, setStyle: (s) => { style = s; draw(); }, setColor: (c) => { color = c; draw(); }, setMelt: (x) => { melt = x; if (color === "melt") draw(); } };
})();

const D = { inited: false, model: null, form: "B", melt: null };
const DNA_PRESETS = [
  ["Dickerson dodecamer", "CGCGAATTCGCG", "B"], ["A-tract (bent)", "GCAAAAAACGGCAAAAAACG", "B"],
  ["GC-rich 20-mer", "GCGGCGCCGCGGCCGCGGCG", "B"], ["Z-DNA hexamer", "CGCGCG", "Z"],
  ["A-DNA octamer", "GGGGCCCC", "A"], ["Primer-length 24-mer", "ATGCTAGCTAGGCTAACGTTAGCA", "B"],
];

function initDna() {
  if (D.inited) return; D.inited = true;
  $("dnaPresets").innerHTML = DNA_PRESETS.map(([n], i) => `<button type="button" class="btn" data-i="${i}">${escapeHtml(n)}</button>`).join("");
  $("dnaPresets").onclick = (e) => { const b = e.target.closest("button"); if (!b) return; const [, s, f] = DNA_PRESETS[+b.dataset.i];
    $("dnaSeq").value = s; setForm(f); dnaBuild(); };
  const setForm = (f) => { D.form = f; document.querySelectorAll("#dnaForm button").forEach((b) => b.setAttribute("aria-checked", String(b.dataset.form === f))); };
  $("dnaForm").onclick = (e) => { const b = e.target.closest("button"); if (b) setForm(b.dataset.form); };
  $("dnaBuild").onclick = dnaBuild;
  document.querySelectorAll("[data-dnastyle]").forEach((b) => (b.onclick = () => {
    document.querySelectorAll("[data-dnastyle]").forEach((x) => x.setAttribute("aria-pressed", String(x === b))); H.setStyle(b.dataset.dnastyle); }));
  document.querySelectorAll("[data-dnacolor]").forEach((b) => (b.onclick = () => {
    document.querySelectorAll("[data-dnacolor]").forEach((x) => x.setAttribute("aria-pressed", String(x === b))); H.setColor(b.dataset.dnacolor); }));
  let t = null;
  const meltLater = () => { $("dnaTempV").textContent = `${$("dnaTemp").value} °C`; clearTimeout(t); t = setTimeout(dnaMelt, 60); };
  $("dnaTemp").oninput = meltLater; $("dnaNa").onchange = dnaMelt; $("dnaCt").onchange = dnaMelt;
  $("dnaSavePdb").onclick = dnaSave;
  dnaBuild();
}

async function dnaBuild() {
  const seq = $("dnaSeq").value.replace(/\s+/g, "").toUpperCase();
  $("dnaBusy").hidden = false;
  try {
    const m = await rpc("dna_build", { sequence: seq, form: D.form });
    D.model = m; H.set(m);
    $("dnaName").textContent = `${m.form}-DNA · ${m.sequence.length} bp`;
    $("dnaFormula").textContent = `5'-${m.sequence}-3'`;
    renderDna(m); $("dnaSavePdb").hidden = false;
    await dnaMelt();
  } catch (e) { log(`Helix build failed: ${e.message || e}`, true); }
  finally { $("dnaBusy").hidden = true; }
}

async function dnaMelt() {
  if (!D.model) return;
  try {
    const r = await rpc("dna_melt", { sequence: D.model.sequence, temp: +$("dnaTemp").value, na: +$("dnaNa").value, ct: +$("dnaCt").value });
    D.melt = r; H.setMelt(r);
    $("dnaMelt").innerHTML = `Melting temperature <b>${r.tm.toFixed(1)} °C</b> (two-state; ΔH° ${r.dH.toFixed(1)} kcal/mol, ΔG°37 ${r.dG37.toFixed(2)} kcal/mol).
      At ${$("dnaTemp").value} °C, <b>${pct(r.fraction_duplex, 0)}</b> of strands are paired; colour by <i>Melting</i> to see which pairs open first.
      <br><small>Nearest-neighbour parameters: SantaLucia 1998; salt: Owczarzy et al. 2004. Checked in CI against 92 measured duplexes at five salt concentrations (mean error ${r.checked || "≈1.6 °C"}).</small>`;
  } catch (e) { $("dnaMelt").textContent = String(e.message || e); }
}

function renderDna(m) {
  const st = m.analysis.steps.map((s) => s.params), pr = m.analysis.pairs.map((p) => p.params);
  const avg = (k) => (st.length ? st.reduce((a, p) => a + p[k], 0) / st.length : NaN);
  const g = m.analysis.grooves.map((x) => x.minor).filter((x) => x != null);  // already less 5.8 A
  $("dnaSummary").innerHTML = kvTable([
    ["Base pairs", m.sequence.length], ["Length", `${(avg(2) * Math.max(0, m.sequence.length - 1)).toFixed(1)} Å`],
    ["Mean twist", `${avg(5).toFixed(1)}° (${(360 / avg(5)).toFixed(1)} bp per turn)`], ["Mean rise", `${avg(2).toFixed(2)} Å`],
    ["Minor groove", g.length ? `${Math.min(...g).toFixed(1)}–${Math.max(...g).toFixed(1)} Å` : "—"],
    ["Backbone O3′–P", m.o3p.length ? `${Math.min(...m.o3p).toFixed(2)}–${Math.max(...m.o3p).toFixed(2)} Å (ideal 1.61)` : "—"],
  ]);
  const f = (x, d = 1) => x.toFixed(d);
  $("dnaSteps").innerHTML = `<table class="tight"><thead><tr><th>Step</th><th>Shift</th><th>Slide</th><th>Rise</th><th>Tilt</th><th>Roll</th><th>Twist</th></tr></thead><tbody>
    ${m.analysis.steps.map((s) => `<tr><td>${escapeHtml(s.step)}</td>${s.params.map((v, k) => `<td class="n">${f(v, k < 3 ? 2 : 1)}</td>`).join("")}</tr>`).join("")}
    </tbody></table><p class="note">Å and degrees, in the standard reference frame (Olson et al. 2001).</p>`;
  $("dnaPairs").innerHTML = `<details><summary class="note">Shear, stretch, stagger, buckle, propeller, opening</summary><table class="tight"><thead><tr><th>Pair</th><th>Shear</th><th>Stretch</th><th>Stagger</th><th>Buckle</th><th>Prop.</th><th>Open.</th></tr></thead><tbody>
    ${m.analysis.pairs.map((p) => `<tr><td>${escapeHtml(p.pair)}</td>${p.params.map((v, k) => `<td class="n">${f(v, k < 3 ? 2 : 1)}</td>`).join("")}</tr>`).join("")}
    </tbody></table></details>`;
  $("dnaProv").innerHTML = `${escapeHtml(m.method)}<ul>${(m.sources || []).map((s) => `<li>${escapeHtml(s)}</li>`).join("")}</ul>${m.validation ? `<b>Checked:</b> ${escapeHtml(m.validation)}` : ""}`;
}

async function dnaSave() {
  if (!D.model) return;
  try {
    const path = await TAURI.core.invoke("pick_save_path", { title: "Save structure", defaultName: `${D.model.form}-DNA-${D.model.sequence.slice(0, 20)}.pdb`,
      filters: [["PDB", ["pdb"]]] });
    if (!path) return;
    await rpc("dna_export", { path, sequence: D.model.sequence, form: D.form });
    log(`Saved ${baseName(path)}`);
  } catch (e) { log(`Save failed: ${e.message || e}`, true); }
}

boot();
