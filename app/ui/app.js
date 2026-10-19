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
};

function methodKey() {
  if (S.tier !== "dft") return S.tier;
  return `dft:${$("funcSel").value}/${$("basisSel").value}`;
}
function methodLabel() {
  if (S.tier === "mmff94") return "MMFF94";
  if (S.tier === "gfn2-xtb") return "GFN2-xTB";
  return `${$("funcSel").selectedOptions[0].textContent}/${$("basisSel").value}`;
}
function engineOptions() {
  if (S.tier !== "dft") return {};
  return { functional: $("funcSel").value, basis: $("basisSel").value };
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
async function runOptimize() {
  const label = methodLabel(), key = methodKey();
  const job = addJob(`Optimize · ${label}`);
  setBusy(true, `Optimizing with ${label}…`);
  try {
    const r = await rpc("optimize", { molecule: currentMolecule(), tier: S.tier, options: engineOptions() });
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
  const pub = S.tier === "dft" ? c.dft.dipole_regularized_rms_percent : c["gfn2-xtb"].dipole_regularized_rms_percent;
  let measured = "";
  const b = S.bench && S.bench.tiers && S.bench.tiers[S.tier] && S.bench.tiers[S.tier].dipoles;
  if (b && b.regularized_rms_percent != null) measured = ` · Athanor benchmark vs experiment: ${b.regularized_rms_percent.toFixed(1)}% (n=${b.n})`;
  const basisNote = S.tier === "dft" ? ` (${c.dft.method}/${pub.basis})` : "";
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
  if (S.tier === "dft" && c && $("funcSel").value === c.method && $("basisSel").value === c.frequency_scale_factor.basis) {
    const sf = c.frequency_scale_factor.value; freqs = freqs.map((f) => f * sf);
    note = `scaled ×${sf} (CCCBDB, ${c.method}/${c.frequency_scale_factor.basis}, σ ${c.frequency_scale_factor.sigma})`;
  } else if (S.tier === "dft") note = "harmonic, unscaled: no literature scale factor for this exact level; harmonic frequencies typically exceed experimental fundamentals by a few percent";
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
    ["gfn2-xtb", "GFN2-xTB", "tblite", ["Energies and frequencies cross-checked against the independent xtb program on every build.",
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
  $("dftOpts").hidden = t !== "dft";
  $("tierNote").textContent = TIER_NOTES[t];
}

function wire() {
  $("presets").innerHTML = PRESETS.map(([n, s]) => `<button type="button" class="btn" data-s="${s}" data-n="${n}">${n}</button>`).join("");
  $("presets").onclick = (e) => { const b = e.target.closest("button"); if (b && !S.busy) build(b.dataset.s, b.dataset.n); };
  $("smilesForm").onsubmit = (e) => { e.preventDefault(); const s = $("smilesInput").value.trim(); if (s && !S.busy) build(s); };
  $("tierSeg").onclick = (e) => { const b = e.target.closest("button"); if (b) setTier(b.dataset.tier); };
  $("btnOpt").onclick = () => runOptimize().catch(() => {});
  $("btnFreq").onclick = runFrequencies;
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
    $("view-molecule").hidden = b.dataset.view !== "molecule";
    $("view-validation").hidden = b.dataset.view !== "validation";
    if (b.dataset.view === "validation") renderValidation();
  };
}

async function boot() {
  wire(); setTier("gfn2-xtb"); setBusy(false);
  try { S.bench = await (await fetch("data/benchmark.json")).json(); } catch { S.bench = null; }
  if (TAURI && TAURI.event) {
    // First launch after install or update unpacks the science engines; show progress.
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
    await build("CCO", "Ethanol");
  } catch (e) {
    $("coreState").textContent = "error"; $("coreChip").querySelector(".dot").className = "dot bad";
    log(`Science core failed to start: ${e.message || e}`, true);
  }
}
boot();
