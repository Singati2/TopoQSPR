/* TopoQSPR v1.0 web front-end.
 * Plain JavaScript, no dependencies, no network calls except this server's /api.
 * Every string that comes from the server or the user is inserted with textContent
 * (DOM building through el()); this file never assigns HTML markup from strings.
 */
"use strict";

(function () {
  const UUID_RE = /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/;
  const state = { meta: null, selected: new Map(), cohort: null, jobId: null, pollTimer: null };

  function $(id) { return document.getElementById(id); }

  /** Build an element: el("td", {className: "num"}, ["text", childNode]). */
  function el(tag, attrs, children) {
    const node = document.createElement(tag);
    if (attrs) {
      for (const [k, v] of Object.entries(attrs)) {
        if (v === null || v === undefined || v === false) continue;
        if (k === "className") node.className = v;
        else if (k === "text") node.textContent = String(v);
        else if (k.startsWith("on") && typeof v === "function") node.addEventListener(k.slice(2), v);
        else node.setAttribute(k, v === true ? "" : String(v));
      }
    }
    for (const c of children || []) {
      if (c === null || c === undefined) continue;
      node.appendChild(typeof c === "string" || typeof c === "number" ? document.createTextNode(String(c)) : c);
    }
    return node;
  }

  function clear(node) { while (node.firstChild) node.removeChild(node.firstChild); }

  function setStatus(id, text, kind) {
    const n = $(id);
    n.textContent = text || "";
    n.className = "status" + (kind ? " " + kind : "");
  }

  async function api(path, options) {
    let resp;
    try {
      resp = await fetch(path, Object.assign({ headers: { "Accept": "application/json" } }, options || {}));
    } catch (e) {
      throw new Error("The server could not be reached. Check your connection and try again.");
    }
    if (!resp.ok) {
      let msg = "Request failed (HTTP " + resp.status + ")";
      try {
        const body = await resp.json();
        if (body && body.detail) msg = typeof body.detail === "string" ? body.detail : JSON.stringify(body.detail);
      } catch (e) { /* not JSON */ }
      const err = new Error(msg);
      err.status = resp.status;
      throw err;
    }
    return resp;
  }

  async function postJSON(path, payload) {
    return api(path, { method: "POST", headers: { "Content-Type": "application/json", "Accept": "application/json" },
                       body: JSON.stringify(payload) });
  }

  function fmt(x, digits) {
    if (x === null || x === undefined || Number.isNaN(x)) return "n/a";
    return Number(x).toFixed(digits === undefined ? 3 : digits);
  }
  function fmtP(x) {
    if (x === null || x === undefined) return "n/a";
    return x < 0.001 ? "<0.001" : Number(x).toFixed(3);
  }
  function fmtCI(lo, hi) {
    if (lo === null || lo === undefined || hi === null || hi === undefined) return "n/a";
    return "[" + fmt(lo) + ", " + fmt(hi) + "]";
  }

  function table(node, headers, rows, caption) {
    clear(node);
    if (caption) node.appendChild(el("caption", { text: caption }));
    node.appendChild(el("thead", null, [el("tr", null, headers.map((h) => el("th", { scope: "col", text: h })))]));
    node.appendChild(el("tbody", null, rows.map((r) => el("tr", null, r.map((c) =>
      c instanceof Node ? el("td", null, [c]) : el("td", { text: c === null || c === undefined ? "" : c }))))));
  }

  function show(id) { $(id).hidden = false; }

  // ------------------------------------------------------------------ meta
  async function loadMeta() {
    const meta = await (await api("/api/meta")).json();
    state.meta = meta;
    $("version").textContent = "v" + meta.version;
    $("mode-note").textContent = meta.mode.note + " Web limits: up to " + meta.limits.max_molecules +
      " molecules per cohort, one job at a time.";
    $("chembl-licence").textContent = meta.chembl.licence;
    const ev = $("evidence-types");
    clear(ev);
    for (const t of meta.evidence_types) {
      const id = "ev-" + t;
      ev.appendChild(el("span", { className: "check" }, [
        el("input", { type: "checkbox", id: id, value: t, name: "evidence" }),
        el("label", { for: id, text: t })]));
    }
    const comp = $("computed-targets");
    clear(comp);
    for (const t of meta.computed_targets) {
      const id = "ct-" + t.name;
      comp.appendChild(el("div", { className: "check" }, [
        el("input", { type: "checkbox", id: id, value: t.name, name: "computed", checked: t.name === "MolLogP" || t.name === "TPSA" }),
        el("label", { for: id }, [el("strong", { text: t.name }), " - " + t.description + " (" + t.units + ", class " + t.target_class + ")"])]));
    }
    const meas = $("measured-targets");
    clear(meas);
    const cites = $("measured-citations");
    clear(cites);
    const seen = new Set();
    for (const t of meta.measured_targets) {
      const id = "mt-" + t.name;
      meas.appendChild(el("div", { className: "check" }, [
        el("input", { type: "checkbox", id: id, value: t.name, name: "measured" }),
        el("label", { for: id }, [el("strong", { text: t.name }), " - " + t.endpoint + " (" + t.units + "); " + t.source + "; licence " + t.licence])]));
      if (!seen.has(t.source)) {
        seen.add(t.source);
        cites.appendChild(el("li", { text: t.source + ". " + t.citation + ". doi:" + t.doi + ". Licence: " + t.licence + "." }));
      }
    }
    if (meta.mesh_allowlist) {
      $("heading-help").textContent = "At least two characters. This server is restricted to the MeSH headings " +
        meta.mesh_allowlist.join(", ") + ".";
    }
  }

  // ------------------------------------------------------------------ headings
  function renderSelected() {
    const ul = $("selected-headings");
    clear(ul);
    for (const [id, h] of state.selected) {
      ul.appendChild(el("li", null, [
        el("span", { text: h + " (" + id + ")" }),
        el("button", { type: "button", className: "small", "aria-label": "Remove " + h,
                       onclick: () => { state.selected.delete(id); renderSelected(); } }, ["Remove"])]));
    }
    $("selected-empty").hidden = state.selected.size > 0;
    $("build-cohort").disabled = state.selected.size === 0;
  }

  async function searchHeadings(ev) {
    ev.preventDefault();
    const q = $("heading-q").value.trim();
    const list = $("heading-results");
    clear(list);
    if (q.length < 2) { setStatus("heading-status", "Type at least two characters.", "error"); return; }
    setStatus("heading-status", "Searching the heading index (the first search after a new ChEMBL release can take a minute)...");
    try {
      const data = await (await api("/api/headings?q=" + encodeURIComponent(q) + "&limit=10")).json();
      setStatus("heading-status", data.headings.length + " suggestion(s) from " + data.release + ".");
      if (!data.headings.length) list.appendChild(el("li", { text: "No heading matched. Try another name or a MeSH term." }));
      for (const h of data.headings) {
        const label = h.mesh_heading + " (" + h.mesh_id + "), " + h.n_rows + " indication rows" +
          (h.matched_via && h.matched_via !== h.mesh_heading ? "; matched via " + h.matched_via : "");
        const btn = el("button", { type: "button", className: "small", disabled: !h.allowed,
                                   "aria-label": "Add " + h.mesh_heading,
                                   onclick: () => { state.selected.set(h.mesh_id, h.mesh_heading); renderSelected(); } },
                       [h.allowed ? "Add" : "Not enabled"]);
        list.appendChild(el("li", null, [btn, " ", el("span", { text: label })]));
      }
    } catch (e) {
      setStatus("heading-status", e.message, "error");
    }
  }

  // ------------------------------------------------------------------ cohort
  async function buildCohort() {
    const ids = Array.from(state.selected.keys());
    if (!ids.length) return;
    const phase = $("min-phase").value;
    const evidence = Array.from(document.querySelectorAll("input[name=evidence]:checked")).map((x) => x.value);
    const payload = {
      mesh_ids: ids,
      min_phase: phase === "none" ? null : Number(phase),
      evidence_sources: evidence.length ? evidence : null,
      exclude_sole_evidence: $("excl-dailymed").checked ? ["DailyMed"] : null,
      exclude_diagnostics: $("excl-diag").checked,
    };
    $("build-cohort").disabled = true;
    setStatus("cohort-status", "Building the cohort from ChEMBL (cached per release; a new cohort can take up to a minute)...");
    try {
      const c = await (await postJSON("/api/cohort", payload)).json();
      state.cohort = c;
      renderCohort(c);
      setStatus("cohort-status", "Cohort ready: " + c.n_included + " included parents. Review it below.", "ok");
      show("step-cohort");
      show("step-run");
      $("confirm").checked = false;
      updateRunButton();
      $("h-cohort").scrollIntoView({ behavior: "smooth" });
    } catch (e) {
      setStatus("cohort-status", e.message, "error");
    } finally {
      $("build-cohort").disabled = state.selected.size === 0;
    }
  }

  function dlRow(dl, k, v) { dl.appendChild(el("dt", { text: k })); dl.appendChild(el("dd", { text: v })); }

  function renderCohort(c) {
    const dl = $("cohort-summary");
    clear(dl);
    dlRow(dl, "Disease", c.disease_label);
    dlRow(dl, "MeSH IDs", c.mesh_ids.join(", "));
    dlRow(dl, "ChEMBL release", c.chembl_release + (c.chembl_release_date ? " (" + c.chembl_release_date + ")" : ""));
    const p = c.parameters;
    dlRow(dl, "Rules", "min phase " + (p.min_phase === null ? "none" : p.min_phase) +
      "; evidence " + (p.evidence_sources ? p.evidence_sources.join("/") : "any") +
      (p.exclude_sole_evidence ? "; excluding sole " + p.exclude_sole_evidence.join("/") : "") +
      (p.exclude_diagnostics ? "; diagnostics excluded" : "") + "; small molecules; InChIKey de-duplication");
    dlRow(dl, "Included / manifest rows", c.n_included + " / " + c.counts.manifest_rows);
    dlRow(dl, "Manifest SHA-256", c.manifest_sha256);
    dlRow(dl, "Data licence", c.data_licence);
    if (!c.runnable) dlRow(dl, "Cannot run on the web", c.not_runnable_reason);
    $("cohort-download").href = "/api/cohorts/" + encodeURIComponent(c.cohort_token) + "/manifest.csv";
    $("chembl-release").textContent = c.chembl_release;
    table($("reason-table"), ["Reason", "Rows"], Object.entries(c.counts.by_reason).map(([k, v]) => [k, String(v)]));
    $("n-included").textContent = String(c.included.length);
    $("n-excluded").textContent = String(c.excluded.length);
    table($("included-table"), ["Parent", "Name", "Phase", "Evidence", "SMILES"],
      c.included.map((r) => [r.parent_chembl_id, r.pref_name || "", fmt(r.phase_for_indication, 1), r.indication_ref_types || "", r.canonical_smiles]));
    table($("excluded-table"), ["Molecule", "Parent", "Name", "Phase", "Reason"],
      c.excluded.map((r) => [r.molecule_chembl_id, r.parent_chembl_id, r.pref_name || "", fmt(r.phase_for_indication, 1), r.reason]));
  }

  // ------------------------------------------------------------------ run
  function chosen(name) {
    return Array.from(document.querySelectorAll("input[name=" + name + "]:checked")).map((x) => x.value);
  }

  function updateRunButton() {
    const ok = state.cohort && state.cohort.runnable && $("confirm").checked &&
      (chosen("computed").length + chosen("measured").length) > 0;
    $("run-button").disabled = !ok;
  }

  async function startRun() {
    if (!state.cohort) return;
    $("run-button").disabled = true;
    setStatus("run-status", "Submitting...");
    try {
      const job = await (await postJSON("/api/run", {
        cohort_token: state.cohort.cohort_token, mode: "exploratory", confirmed: $("confirm").checked,
        targets: { computed: chosen("computed"), measured: chosen("measured") },
      })).json();
      setStatus("run-status", "Job submitted.", "ok");
      openJob(job.job_id);
    } catch (e) {
      setStatus("run-status", e.message, "error");
      updateRunButton();
    }
  }

  // ------------------------------------------------------------------ jobs
  function openJob(id) {
    if (!UUID_RE.test(id)) return;
    state.jobId = id;
    if (location.hash !== "#job=" + id) history.replaceState(null, "", "#job=" + id);
    const link = $("job-link");
    link.href = "#job=" + id;
    link.textContent = location.origin + location.pathname + "#job=" + id;
    $("results").hidden = true;
    show("step-results");
    $("h-results").scrollIntoView({ behavior: "smooth" });
    poll(0);
  }

  function poll(delay) {
    if (state.pollTimer) clearTimeout(state.pollTimer);
    state.pollTimer = setTimeout(checkJob, delay);
  }

  async function checkJob() {
    const id = state.jobId;
    if (!id) return;
    try {
      const job = await (await api("/api/jobs/" + id)).json();
      if (id !== state.jobId) return;
      const prog = $("job-progress");
      prog.value = job.progress || 0;
      prog.textContent = Math.round((job.progress || 0) * 100) + "%";
      let stage;
      if (job.status === "queued") stage = "Queued" + (job.queue_position ? " (position " + job.queue_position + ")" : "");
      else if (job.status === "running") stage = "Running: " + (job.message || job.stage || "");
      else if (job.status === "done") stage = "Finished";
      else stage = "Failed";
      $("job-stage").textContent = stage + " - " + Math.round((job.progress || 0) * 100) + "%";
      if (job.status === "done") {
        setStatus("job-status", "Finished at " + job.finished_utc + ".", "ok");
        await loadResults(id);
        return;
      }
      if (job.status === "failed") {
        setStatus("job-status", "The job failed: " + (job.error || "unknown reason"), "error");
        return;
      }
      setStatus("job-status", job.cohort ? "Cohort: " + job.cohort.disease_label + " (" + job.cohort.n_included + " parents)." : "");
      poll(job.status === "queued" ? 4000 : 2000);
    } catch (e) {
      setStatus("job-status", e.message, "error");
      if (e.status !== 404) poll(8000);
    }
  }

  async function loadResults(id) {
    const data = await (await api("/api/jobs/" + id + "/results")).json();
    const dl = $("results-summary");
    clear(dl);
    const job = data.job;
    if (job.cohort) {
      dlRow(dl, "Cohort", job.cohort.disease_label + " (" + job.cohort.mesh_ids.join(", ") + ")");
      dlRow(dl, "ChEMBL release", job.cohort.chembl_release);
    }
    dlRow(dl, "Molecules", (data.counts.standardised === null ? "n/a" : data.counts.standardised) +
      " standardised of " + data.counts.candidates);
    const pr = data.protocol;
    dlRow(dl, "Protocol", "exploratory; learner " + pr.learner + " over " + pr.families.join(", ") +
      "; " + pr.outer_splits + "-fold outer x " + pr.repeats + " repeat; " + pr.inner_splits + "-fold inner; AD " + pr.ad);
    dlRow(dl, "Code", (data.code && data.code.commit ? data.code.commit.slice(0, 12) + (data.code.dirty ? " (modified)" : "") : "not recorded") +
      "; topoqspr " + data.topoqspr_version);
    dlRow(dl, "Config SHA-256", data.config_sha256);
    $("dl-zip").href = "/api/jobs/" + id + "/results.zip";
    $("dl-manifest").href = "/api/jobs/" + id + "/manifest.csv";
    $("dl-record").href = "/api/jobs/" + id + "/run_record.json";
    const box = $("ladder");
    clear(box);
    for (const t of data.targets) {
      const sec = el("section", { className: "target" }, [
        el("h3", { text: t.name + " (" + (t.units || "") + ")" }),
        el("p", { className: "help", text: (t.description || "") + ". Class " + t.target_class + " (" + (t.target_class_label || "") +
          "). Values for " + (t.n_values === null || t.n_values === undefined ? "?" : t.n_values) + " molecules." })]);
      if (t.skipped_reason) {
        sec.appendChild(el("p", { className: "status error", text: "Not evaluated: " + t.skipped_reason }));
        box.appendChild(sec);
        continue;
      }
      const headers = ["Feature set", "n", "p", "R² (fold mean)", "R² (pooled)", "MAE", "ΔR² vs size", "95% CI", "p (ΔR²)",
                       "ΔMAE vs size", "95% CI", "Folds valid/failed", "Underpowered", "AD outside", "Note"];
      const rows = t.rows.map((r) => [
        r.set + (r.is_reference ? " (reference)" : ""), fmt(r.n, 0), fmt(r.p, 0), fmt(r.r2_fold_mean), fmt(r.r2_pooled_mean),
        fmt(r.mae_fold_mean), r.is_reference ? "-" : fmt(r.d_r2_mean), r.is_reference ? "-" : fmtCI(r.d_r2_ci_low, r.d_r2_ci_high),
        r.is_reference ? "-" : fmtP(r.d_r2_p), r.is_reference ? "-" : fmt(r.d_mae_mean),
        r.is_reference ? "-" : fmtCI(r.d_mae_ci_low, r.d_mae_ci_high),
        fmt(r.n_valid_folds, 0) + "/" + fmt(r.n_failed_folds, 0), r.underpowered ? "yes" : "no",
        r.ad_share_outside === null ? "n/a" : Math.round(r.ad_share_outside * 100) + "%",
        [r.reason, r.underpowered_reason, r.ad_reason].filter(Boolean).join("; ")]);
      const tbl = el("table");
      table(tbl, headers, rows, "Baseline ladder for " + t.name + " (learner " + t.learner + ")");
      sec.appendChild(el("div", { className: "table-wrap" }, [tbl]));
      box.appendChild(sec);
    }
    $("results").hidden = false;
  }

  // ------------------------------------------------------------------ indices
  async function loadFile() {
    const f = $("indices-file").files[0];
    if (!f) return;
    if (f.size > (state.meta ? state.meta.limits.max_body_bytes : 200000)) {
      setStatus("indices-status", "The file is too large for the web limit.", "error");
      return;
    }
    $("indices-text").value = await f.text();
    setStatus("indices-status", "Loaded " + f.name + ".");
  }

  function indicesPayload(view) {
    return { text: $("indices-text").value, format: $("indices-format").value,
             entropies: $("indices-entropies").checked, view: view };
  }

  async function computeIndices() {
    const text = $("indices-text").value;
    if (!text.trim()) { setStatus("indices-status", "Paste structures or choose a file first.", "error"); return; }
    $("indices-button").disabled = true;
    setStatus("indices-status", "Computing...");
    try {
      const resp = await postJSON("/api/indices", indicesPayload(false));
      const blob = await resp.blob();
      const url = URL.createObjectURL(blob);
      const a = el("a", { href: url, download: "topoqspr-indices.zip" });
      document.body.appendChild(a);
      a.click();
      a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 10000);
      setStatus("indices-status", resp.headers.get("X-Molecules-OK") + " of " + resp.headers.get("X-Molecules") +
        " structures accepted; " + resp.headers.get("X-Undefined-Values") + " undefined values (reasons in the zip).", "ok");
    } catch (e) {
      setStatus("indices-status", e.message, "error");
    } finally {
      $("indices-button").disabled = false;
    }
  }

  const META_COLS = ["id", "input_smiles", "status", "reason", "smiles", "inchikey", "kind", "n_vertices", "n_edges"];

  function renderIndexView() {
    const v = state.indexView;
    if (!v) return;
    const ids = v.rows.map((r) => r[0]);
    const q = $("indices-filter").value.trim().toLowerCase();
    const rows = [];
    v.columns.forEach((col, ci) => {
      if (ci === 0) return;
      if (q && META_COLS.indexOf(col) < 0 && col.toLowerCase().indexOf(q) < 0) return;
      rows.push([col].concat(v.rows.map((r) => {
        const x = r[ci];
        return typeof x === "number" ? (Number.isInteger(x) ? String(x) : x.toPrecision(8)) : (x === null ? "undefined" : x);
      })));
    });
    table($("indices-table"), ["index"].concat(ids), rows,
      v.n_ok + " of " + v.n + " structures accepted; " + v.failures.length + " undefined values or rejections");
  }

  async function viewIndices() {
    const text = $("indices-text").value;
    if (!text.trim()) { setStatus("indices-status", "Paste structures or choose a file first.", "error"); return; }
    $("indices-view-button").disabled = true;
    setStatus("indices-status", "Computing...");
    try {
      const resp = await postJSON("/api/indices", indicesPayload(true));
      state.indexView = await resp.json();
      renderIndexView();
      setStatus("indices-status", "Done.", "ok");
    } catch (e) {
      setStatus("indices-status", e.message, "error");
    } finally {
      $("indices-view-button").disabled = false;
    }
  }

  function renderCatalog() {
    const items = state.catalog || [];
    const q = $("catalog-filter").value.trim().toLowerCase();
    const rows = items.filter((i) => !q || (i.name + " " + i.family + " " + i.reference).toLowerCase().indexOf(q) >= 0)
      .map((i) => [i.name, i.family, i.formula_latex, i.reference]);
    table($("catalog-table"), ["index", "family", "formula (LaTeX)", "primary reference"], rows,
      rows.length + " of " + items.length + " entries");
  }

  async function loadCatalog() {
    try {
      const resp = await api("/api/catalog");
      state.catalog = (await resp.json()).items;
      renderCatalog();
    } catch (e) { /* the catalogue is optional; the rest of the page works without it */ }
  }

  // ------------------------------------------------------------------ init
  function routeHash() {
    const m = /^#job=([0-9a-f-]{36})$/.exec(location.hash);
    if (m && m[1] !== state.jobId) openJob(m[1]);
  }

  document.addEventListener("DOMContentLoaded", async () => {
    $("heading-form").addEventListener("submit", searchHeadings);
    $("build-cohort").addEventListener("click", buildCohort);
    $("run-button").addEventListener("click", startRun);
    $("confirm").addEventListener("change", updateRunButton);
    $("computed-targets").addEventListener("change", updateRunButton);
    $("measured-targets").addEventListener("change", updateRunButton);
    $("indices-file").addEventListener("change", loadFile);
    $("indices-button").addEventListener("click", computeIndices);
    $("indices-view-button").addEventListener("click", viewIndices);
    $("indices-filter").addEventListener("input", renderIndexView);
    $("catalog-filter").addEventListener("input", renderCatalog);
    loadCatalog();
    window.addEventListener("hashchange", routeHash);
    renderSelected();
    try {
      await loadMeta();
    } catch (e) {
      setStatus("heading-status", "Could not load the server settings: " + e.message, "error");
    }
    routeHash();
  });
})();
