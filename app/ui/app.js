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
  if (d.task === "adna") return `Ancient DNA: ${d.stage}`;
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
  $("adnaCov").oninput = (e) => ($("adnaCovVal").textContent = `${e.target.value}×`);
  $("adnaCont").oninput = (e) => ($("adnaContVal").textContent = `${e.target.value}%`);
  $("adnaLen").oninput = (e) => ($("adnaLenVal").textContent = `${e.target.value} bp`);
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
    for (const v of ["molecule", "validation", "adna"]) $(`view-${v}`).hidden = b.dataset.view !== v;
    if (b.dataset.view === "validation") renderValidation();
    if (b.dataset.view === "adna") initAdna();
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
const A = { catalog: null, result: null, busy: false };

async function initAdna() {
  const linuxReady = S.compute && S.compute.stage === "ready";
  $("adnaRun").disabled = !linuxReady || A.busy;
  $("adnaEmpty").textContent = linuxReady
    ? "Choose a sample and references, then run."
    : "The Ancient DNA Lab runs in the Linux engines pack. Install or update it from the top bar first.";
  if (!linuxReady || A.catalog) return;
  try {
    A.catalog = await rpc("adna_catalog", {});
    const g = A.catalog.references.genomes;
    $("adnaEndo").innerHTML = g.map((x) => `<option value="${x.accession}" ${x.accession === "NC_007596" ? "selected" : ""}>${escapeHtml(x.organism)} — ${x.version}</option>`).join("");
    $("adnaRefs").innerHTML = g.map((x) => `<label><input type="checkbox" value="${x.accession}" ${["NC_007596", "NC_000934", "NC_005129"].includes(x.accession) ? "checked" : ""}>
      <span>${escapeHtml(x.organism)}<small>${x.version} · ${x.length.toLocaleString()} bp</small></span></label>`).join("");
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
  };
}

async function runAdna() {
  const p = adnaParams();
  if (!p.references.length) { log("Choose at least one reference", true); return; }
  A.busy = true; $("adnaRun").disabled = true; $("adnaBusy").hidden = false; $("adnaBusyText").textContent = "Starting…";
  const job = addJob(`Ancient DNA · ${p.coverage}× · ${p.references.length} reference(s)`);
  try {
    A.result = await rpc("adna_run", p);
    renderAdna(); job.done(`${A.result.references.length} reference(s)`);
    const best = A.result.references.reduce((a, b) => (a.accuracy > b.accuracy ? a : b));
    log(`Ancient DNA: best reconstruction ${(100 * best.accuracy).toFixed(3)}% correct against ${best.reference.organism}`);
  } catch (e) { job.fail(String(e)); log(`Ancient DNA run failed: ${e.message || e}`, true); }
  finally { A.busy = false; $("adnaRun").disabled = false; $("adnaBusy").hidden = true; }
}

const pct = (x, d = 2) => (x == null ? "—" : `${(100 * x).toFixed(d)}%`);
function renderAdna() {
  const R = A.result; if (!R) return;
  $("adnaEmpty").hidden = true; $("adnaResults").hidden = false;
  const rows = R.references.map((r) => {
    const c = r.calibration;
    return `<tr><td>${escapeHtml(r.reference.organism)}<br><small>${r.reference.version}</small></td>
      <td class="n">${pct(r.truth_vs_reference_divergence)}</td><td class="n">${pct(r.estimated_divergence_theta)}</td>
      <td class="n">${(r.mapped_reads ?? 0).toLocaleString()}</td><td class="n">${r.mean_depth.toFixed(1)}×</td>
      <td class="n">${pct(r.covered_fraction, 1)}</td><td class="n">${pct(r.accuracy, 3)}</td>
      <td class="n">${pct(c.predicted_error_rate, 3)}</td><td class="n">${pct(c.observed_error_rate, 3)}</td>
      <td class="n">${c.expected_calibration_error.toFixed(4)}</td></tr>`;
  }).join("");
  $("adnaTable").innerHTML = `<table><thead><tr><th>Reference</th><th>Truth vs reference</th><th>Estimated divergence</th>
    <th>Mapped reads</th><th>Depth</th><th>Covered</th><th>Consensus correct</th><th>Predicted error</th><th>Observed error</th><th>Calibration error</th></tr></thead><tbody>${rows}</tbody></table>
    <p class="note">"Predicted error" is what the posteriors claim; "observed" is measured against the true genome. Calibration error near 0 means the stated confidence can be trusted. Positions where the truth has an indel relative to the reference are not scored.</p>`;
  const sim = R.simulation.damage_briggs_v_l_d_s;
  const names = [["Nu", "ν nick frequency", sim[0]], ["Lambda", "λ overhang parameter", sim[1]], ["DeltaD", "δd double-stranded", sim[2]], ["DeltaS", "δs single-stranded", sim[3]]];
  const drows = names.map(([k, lab, v]) => `<tr><td>${lab}</td><td class="n">${v}</td>${R.references.map((r) => {
    const e = r.mapdamage_bayesian && r.mapdamage_bayesian[k];
    return `<td class="n">${e && e.Mean != null ? e.Mean.toFixed(4) : "—"}</td>`; }).join("")}</tr>`).join("");
  $("adnaDamageTable").innerHTML = `<table><thead><tr><th>Briggs parameter</th><th>Simulated</th>${R.references.map((r) => `<th>${escapeHtml(r.reference.organism.split(" (")[0])}</th>`).join("")}</tr></thead><tbody>${drows}</tbody></table>`;
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
  });
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
  $("adnaProv").innerHTML = `<p><b>Simulation</b>: ${escapeHtml(JSON.stringify(R.simulation))}</p>
    <p><b>Tools</b></p><ul>${tools}</ul><p><b>Reference genomes</b> (${escapeHtml(R.reference_catalog.source)}, retrieved ${escapeHtml(R.reference_catalog.retrieved_utc)})</p><ul>${refs}</ul>
    <p><b>Steps</b></p><ol>${steps}</ol>`;
}

async function boot() {
  wire(); setTier("gfn2-xtb"); setBusy(false);
  try { S.bench = await (await fetch("data/benchmark.json")).json(); } catch { S.bench = null; }
  if (TAURI && TAURI.event) {
    // First launch after install or update unpacks the science engines; show progress.
    TAURI.event.listen("core-progress", (e) => {
      const d = e.payload || {}, t = progressText(d);
      if (d.task === "adna") { if (A.busy) $("adnaBusyText").textContent = t; return; }
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
boot();
