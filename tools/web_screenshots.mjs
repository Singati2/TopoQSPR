// Drive headless Chrome over CDP to take screenshots of the TopoQSPR web interface.
import { spawn } from "node:child_process";
import { writeFileSync } from "node:fs";

const OUT = "/home/mpcrlab/topoqspr_web_data/screens";
const PORT = 9333;
const chrome = spawn("google-chrome", ["--headless=new", "--disable-gpu", "--no-sandbox", "--hide-scrollbars",
  `--remote-debugging-port=${PORT}`, "--window-size=1300,1000", "--user-data-dir=/tmp/topoqspr-shoot-profile", "about:blank"],
  { stdio: "ignore" });
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

let ws, id = 0;
const pending = new Map();
function send(method, params = {}) {
  return new Promise((resolve, reject) => {
    const i = ++id;
    pending.set(i, { resolve, reject });
    ws.send(JSON.stringify({ id: i, method, params }));
  });
}
async function evalJS(expr) {
  const r = await send("Runtime.evaluate", { expression: expr, awaitPromise: true, returnByValue: true });
  if (r.exceptionDetails) throw new Error(JSON.stringify(r.exceptionDetails).slice(0, 400));
  return r.result.value;
}
async function waitFor(expr, timeoutMs = 120000) {
  const t0 = Date.now();
  while (Date.now() - t0 < timeoutMs) {
    if (await evalJS(expr)) return true;
    await sleep(500);
  }
  throw new Error("timeout waiting for " + expr);
}
async function shot(name, selector) {
  const box = await evalJS(`(() => { const e = document.querySelector(${JSON.stringify(selector)});
    e.scrollIntoView(); const r = e.getBoundingClientRect();
    return {x: r.left + scrollX, y: r.top + scrollY, w: r.width, h: r.height}; })()`);
  const r = await send("Page.captureScreenshot", { format: "png", captureBeyondViewport: true,
    clip: { x: Math.max(0, box.x - 12), y: Math.max(0, box.y - 12), width: box.w + 24, height: Math.min(box.h + 24, 3000), scale: 1 } });
  writeFileSync(`${OUT}/${name}.png`, Buffer.from(r.data, "base64"));
  console.log("saved", name);
}

for (let i = 0; i < 40; i++) {
  try { const v = await (await fetch(`http://127.0.0.1:${PORT}/json/version`)).json(); if (v) break; } catch { await sleep(250); }
}
const tab = await (await fetch(`http://127.0.0.1:${PORT}/json/new?http://localhost:8000/`, { method: "PUT" })).json();
ws = new WebSocket(tab.webSocketDebuggerUrl);
ws.onmessage = (m) => { const d = JSON.parse(m.data); if (d.id && pending.has(d.id)) { const p = pending.get(d.id); pending.delete(d.id); d.error ? p.reject(new Error(d.error.message)) : p.resolve(d.result); } };
await new Promise((r) => (ws.onopen = r));
await send("Emulation.setEmulatedMedia", { features: [{ name: "prefers-color-scheme", value: "light" }] });
await send("Emulation.setDeviceMetricsOverride", { width: 1300, height: 1000, deviceScaleFactor: 1, mobile: false });
await send("Page.enable");
await send("Page.navigate", { url: "http://localhost:8000/" });
await waitFor("document.readyState === 'complete' && document.querySelectorAll('#catalog-table tr').length > 5");

const step = process.argv[2] || "all";
try {
  if (step === "all" || step === "indices") {
    await shot("02_disease_step", "#step-disease");
    await evalJS(`document.getElementById('indices-text').value = 'CC(=O)Oc1ccccc1C(=O)O aspirin\\nCC(C)Cc1ccc(cc1)C(C)C(=O)O ibuprofen\\nCn1cnc2c1c(=O)n(C)c(=O)n2C caffeine\\nCN(C)C(=N)N=C(N)N.Cl metformin_HCl'; true`);
    await evalJS(`document.getElementById('indices-view-button').click(); true`);
    await waitFor("document.querySelectorAll('#indices-table tbody tr').length > 20");
    await sleep(2500); // depictions
    await evalJS(`(() => { const f = document.getElementById('indices-filter'); f.value = 'M1'; f.dispatchEvent(new Event('input')); return true; })()`);
    await sleep(800);
    await shot("03_indices_drugs_filtered_M1", "#step-indices");
    await evalJS(`(() => { const f = document.getElementById('indices-filter'); f.value = ''; f.dispatchEvent(new Event('input')); return true; })()`);
    await evalJS(`(() => { const s = document.getElementById('indices-format'); s.value = 'adjacency_matrix'; s.dispatchEvent(new Event('change'));
      document.getElementById('indices-text').value = 'A = [[0,1,0,0,0,1],[1,0,1,0,0,0],[0,1,0,1,0,0],[0,0,1,0,1,0],[0,0,0,1,0,1],[1,0,0,0,1,0]]'; return true; })()`);
    await evalJS(`document.getElementById('indices-view-button').click(); true`);
    await waitFor("document.querySelector('#indices-table caption') && document.querySelector('#indices-table caption').textContent.indexOf('1 of 1') === 0");
    await evalJS(`(() => { const f = document.getElementById('indices-filter'); f.value = 'ENT_'; f.dispatchEvent(new Event('input')); return true; })()`);
    await sleep(800);
    await shot("04_indices_adjacency_C6_entropies", "#step-indices");
    await evalJS(`(() => { const f = document.getElementById('catalog-filter'); f.value = 'reverse'; f.dispatchEvent(new Event('input')); return true; })()`);
    await sleep(500);
    await shot("05_catalogue_reverse_degree", "#catalog");
  }
  if (step === "all" || step === "cohort") {
    await evalJS(`document.getElementById('heading-q').value = 'migraine'; document.getElementById('heading-form').requestSubmit(); true`);
    await waitFor("document.querySelectorAll('#heading-results button').length > 0", 180000);
    await shot("06_mesh_suggestions", "#step-disease");
    await evalJS(`(() => { const b = [...document.querySelectorAll('#heading-results button')].find(x => !x.disabled && /Migraine Disorders/.test(x.closest('li') ? x.closest('li').textContent : x.textContent)); (b || document.querySelector('#heading-results button')).click(); return true; })()`);
    await evalJS(`document.getElementById('build-cohort').click(); true`);
    await waitFor("!document.getElementById('step-cohort').hidden && document.querySelectorAll('#included-table tbody tr').length > 0", 300000);
    await sleep(800);
    await shot("07_cohort_manifest_review", "#step-cohort");
    await evalJS(`(() => { const c = document.querySelector('#computed-targets input[type=checkbox]'); if (c) { c.checked = true; c.dispatchEvent(new Event('change', {bubbles: true})); }
      const k = document.getElementById('confirm'); k.checked = true; k.dispatchEvent(new Event('change', {bubbles: true})); return true; })()`);
    await shot("08_targets_and_confirm", "#step-run");
    await evalJS(`document.getElementById('run-button').click(); true`);
    await waitFor("!document.getElementById('results').hidden && document.querySelectorAll('#ladder tr, #ladder table tr').length > 2", 1800000);
    await sleep(1000);
    await shot("09_results_ladder", "#step-results");
  }
} catch (e) {
  console.error("ERROR", e.message);
  await send("Page.captureScreenshot", { format: "png" }).then((r) => writeFileSync(`${OUT}/zz_error_state.png`, Buffer.from(r.data, "base64")));
} finally {
  chrome.kill();
  process.exit(0);
}
