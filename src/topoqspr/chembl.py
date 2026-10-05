"""Thin, cached REST clients for ChEMBL and NLM MeSH.

ChEMBL web services: https://www.ebi.ac.uk/chembl/api/data (Davies et al., Nucleic Acids
Res. 43 (2015) W612-W620, doi:10.1093/nar/gkv352; Zdrazil et al., Nucleic Acids Res. 52
(2024) D1180-D1192, doi:10.1093/nar/gkad1004).
NLM MeSH RDF / lookup API: https://id.nlm.nih.gov/mesh/ (MeSH linked data, NLM).

Design rules
------------
* Every ChEMBL cache entry is keyed by ``(chembl_db_version, endpoint, params)``; entries
  live under ``<cache_dir>/<chembl_db_version>/`` so two releases can never mix.  The
  release is read from ``/status`` once per client (or, offline, from the cached
  ``status.json``).
* An online client writes new entries to a private staging directory
  ``<cache_dir>/.staging-*/`` (visible only to that client) and promotes them into
  ``<cache_dir>/<release>/`` only when :meth:`ChemblClient.verify_release` confirms that
  the live release is still the one read at start.  If the release changed mid-run the
  staged entries are deleted, so data served by release N+1 is never filed under N.
  Entries of a client that never calls ``verify_release`` are never promoted.
* Cache files are gzip-compressed canonical JSON (sorted keys, ``mtime=0``), so the same
  response always produces the same bytes; each file stores its own key and is verified on
  load (hash-collision / tampering guard).
* Transient failures (connection errors, timeouts, HTTP 429/5xx) are retried with
  exponential backoff ``backoff * 2**attempt`` seconds (``Retry-After`` honoured); other
  HTTP errors raise immediately.  Nothing is ever silently replaced by an empty result.
* ``offline=True`` never touches the network: a cache miss raises :class:`CacheMiss`.
  The test-suite runs offline from recorded fixtures.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

import requests

__all__ = [
    "CHEMBL_BASE_URL",
    "MESH_BASE_URL",
    "CHEMBL_DATA_LICENCE",
    "CacheMiss",
    "ChemblClient",
    "MeshClient",
    "chembl_sort_key",
    "default_cache_dir",
]

CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
MESH_BASE_URL = "https://id.nlm.nih.gov/mesh"

#: Licence of ChEMBL data and of files derived from it (recorded in cohort
#: sidecars and run records; redistributed ChEMBL-derived files must carry it).
CHEMBL_DATA_LICENCE = ("CC BY-SA 3.0 (https://creativecommons.org/licenses/by-sa/3.0/); "
                       "ChEMBL, EMBL-EBI")

#: Fields kept from a ChEMBL ``molecule`` record (the molfile is dropped to keep the cache small).
MOLECULE_FIELDS = (
    "molecule_chembl_id",
    "pref_name",
    "molecule_type",
    "structure_type",
    "inorganic_flag",
    "max_phase",
    "atc_classifications",
    "molecule_hierarchy",
    "molecule_structures",
)
_STRUCTURE_FIELDS = ("canonical_smiles", "standard_inchi_key")
_RETRY_STATUS = {429, 500, 502, 503, 504}


class CacheMiss(LookupError):
    """Raised by an offline client when a request is not in the cache."""


def default_cache_dir() -> Path:
    """``$TOPOQSPR_CACHE`` if set, else ``$XDG_CACHE_HOME/topoqspr`` (``~/.cache/topoqspr``)."""
    env = os.environ.get("TOPOQSPR_CACHE")
    if env:
        return Path(env)
    base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    return Path(base) / "topoqspr"


def _canonical(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _write_gz_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _canonical(obj).encode("utf-8")
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    with open(tmp, "wb") as fh:
        # mtime=0 and no filename -> byte-identical files for identical content
        with gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0) as gz:
            gz.write(raw)
    os.replace(tmp, path)


def _read_gz_json(path: Path) -> Any:
    with gzip.open(path, "rb") as gz:
        return json.loads(gz.read().decode("utf-8"))


def chembl_sort_key(chembl_id: str) -> tuple[int, str]:
    """Sort key: numeric part of a ChEMBL ID, then the full string (CHEMBL2 < CHEMBL10)."""
    digits = "".join(ch for ch in str(chembl_id) if ch.isdigit())
    return (int(digits) if digits else -1, str(chembl_id))


class _HttpJson:
    """Shared GET-with-retries logic."""

    def __init__(
        self,
        base_url: str,
        session: requests.Session | None,
        timeout: float,
        max_retries: int,
        backoff: float,
        sleep: Callable[[float], None],
        offline: bool,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.session = session if session is not None else requests.Session()
        self.timeout = timeout
        self.max_retries = max_retries
        self.backoff = backoff
        self._sleep = sleep
        self.offline = offline
        self.n_http_requests = 0

    def _get_json(self, path: str, params: Mapping[str, Any] | None = None) -> Any:
        if self.offline:
            raise CacheMiss(f"offline client cannot fetch {path} {dict(params or {})}")
        url = f"{self.base_url}/{path.lstrip('/')}"
        last_exc: Exception | None = None
        for attempt in range(self.max_retries + 1):
            wait = self.backoff * (2**attempt)
            try:
                self.n_http_requests += 1
                resp = self.session.get(
                    url,
                    params=dict(params or {}),
                    timeout=self.timeout,
                    headers={"Accept": "application/json"},
                )
            except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as exc:
                last_exc = exc
            else:
                if resp.status_code == 200:
                    try:
                        return resp.json()
                    except ValueError as exc:  # e.g. an HTML maintenance page served with 200
                        last_exc = exc
                        if attempt < self.max_retries:
                            self._sleep(wait)
                        continue
                if resp.status_code not in _RETRY_STATUS:
                    resp.raise_for_status()
                    raise requests.HTTPError(f"HTTP {resp.status_code} for {resp.url}")
                last_exc = requests.HTTPError(f"HTTP {resp.status_code} for {resp.url}")
                retry_after = resp.headers.get("Retry-After")
                if retry_after and retry_after.isdigit():
                    wait = max(wait, float(retry_after))
            if attempt < self.max_retries:
                self._sleep(wait)
        raise RuntimeError(f"GET {url} failed after {self.max_retries + 1} attempts") from last_exc


class ChemblClient(_HttpJson):
    """Cached, paginating client for the ChEMBL REST API.

    Parameters
    ----------
    cache_dir:
        Cache root (default :func:`default_cache_dir` ``/chembl``).  ``None`` disables
        caching only when ``offline`` is False.
    offline:
        Serve from cache only; raise :class:`CacheMiss` otherwise.
    release:
        Offline only: which cached release to use (default: the release in
        ``<cache_dir>/status.json``).
    page_limit:
        Page size for paginated endpoints (ChEMBL maximum is 1000).
    """

    def __init__(
        self,
        cache_dir: str | os.PathLike | None = None,
        *,
        offline: bool = False,
        release: str | None = None,
        base_url: str = CHEMBL_BASE_URL,
        session: requests.Session | None = None,
        timeout: float = 60.0,
        max_retries: int = 5,
        backoff: float = 1.0,
        page_limit: int = 1000,
        molecule_batch: int = 50,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        super().__init__(base_url, session, timeout, max_retries, backoff, sleep, offline)
        self.cache_dir = Path(cache_dir) if cache_dir is not None else default_cache_dir() / "chembl"
        self.page_limit = int(page_limit)
        self.molecule_batch = int(molecule_batch)
        self._requested_release = release
        self._status: dict | None = None
        self._staging: Path | None = None
        self._staged: list[Path] = []

    def __enter__(self) -> "ChemblClient":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        """Verify-and-promote on a clean exit; discard staged entries on an exception."""
        if exc_type is not None:
            self.discard_staged()
        else:
            self.verify_release()

    # ------------------------------------------------------------------ status / release
    def status(self) -> dict:
        """``/status`` payload (fetched once per client; offline: the cached copy)."""
        if self._status is not None:
            return self._status
        if self.offline:
            path = (
                self.cache_dir / self._requested_release / "status.json"
                if self._requested_release
                else self.cache_dir / "status.json"
            )
            if not path.exists():
                raise CacheMiss(f"no cached ChEMBL status at {path}")
            status = json.loads(path.read_text())
        else:
            status = self._get_json("status.json")
            if self._requested_release and status.get("chembl_db_version") != self._requested_release:
                raise RuntimeError(
                    f"requested ChEMBL release {self._requested_release!r} but the live API "
                    f"serves {status.get('chembl_db_version')!r}"
                )
            # keep only release-identifying fields: the row counts drift with no meaning here
            status = {
                k: status[k] for k in ("chembl_db_version", "chembl_release_date") if k in status
            }
            for path in (self.cache_dir / "status.json", self.cache_dir / status["chembl_db_version"] / "status.json"):
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(status, sort_keys=True, indent=1) + "\n")
        if "chembl_db_version" not in status:
            raise RuntimeError(f"ChEMBL status lacks chembl_db_version: {status}")
        self._status = status
        return status

    @property
    def release(self) -> str:
        """ChEMBL release string, e.g. ``'ChEMBL_37'``."""
        return self.status()["chembl_db_version"]

    # ------------------------------------------------------------------ cache
    def _cache_path(self, endpoint: str, params: Mapping[str, Any]) -> tuple[Path, dict]:
        key = {"release": self.release, "endpoint": endpoint, "params": dict(sorted(params.items()))}
        digest = hashlib.sha256(_canonical(key).encode()).hexdigest()[:24]
        safe_ep = endpoint.replace("/", "_")
        return self.cache_dir / self.release / safe_ep / f"{digest}.json.gz", key

    def _staged_path(self, path: Path) -> Path:
        if self._staging is None:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            self._staging = Path(tempfile.mkdtemp(prefix=".staging-", dir=self.cache_dir))
        return self._staging / path.relative_to(self.cache_dir)

    def _load(self, path: Path, key: dict) -> tuple[bool, Any]:
        """``(True, data)`` from the promoted cache or this client's staging area."""
        candidates = [path]
        if self._staging is not None:
            candidates.append(self._staging / path.relative_to(self.cache_dir))
        for cand in candidates:
            if cand.exists():
                blob = _read_gz_json(cand)
                if blob.get("key") != key:
                    raise RuntimeError(f"cache key mismatch in {cand}")
                return True, blob["data"]
        return False, None

    def _store(self, path: Path, key: dict, data: Any) -> None:
        """Write to staging; promoted by :meth:`verify_release`."""
        staged = self._staged_path(path)
        _write_gz_json(staged, {"key": key, "data": data})
        self._staged.append(path)

    def _cached(self, endpoint: str, params: Mapping[str, Any], producer: Callable[[], Any]) -> Any:
        params = {k: (str(v) if not isinstance(v, str) else v) for k, v in params.items()}
        path, key = self._cache_path(endpoint, params)
        hit, data = self._load(path, key)
        if hit:
            return data
        if self.offline:
            raise CacheMiss(f"{endpoint} {params} not cached for {self.release} under {self.cache_dir}")
        data = producer()
        self._store(path, key, data)
        return data

    def discard_staged(self) -> int:
        """Delete this client's unpromoted cache entries; returns how many were dropped."""
        n = len(self._staged)
        if self._staging is not None:
            shutil.rmtree(self._staging, ignore_errors=True)
        self._staging, self._staged = None, []
        return n

    def _promote(self) -> None:
        for path in sorted(set(self._staged)):
            staged = self._staging / path.relative_to(self.cache_dir)
            path.parent.mkdir(parents=True, exist_ok=True)
            os.replace(staged, path)
        self.discard_staged()

    def verify_release(self) -> None:
        """Re-read ``/status``; raise if the live release changed since this client started.

        On success the entries staged by this client are promoted into
        ``<cache_dir>/<release>/``; on a change they are deleted first, so nothing fetched
        after a switch is ever filed under the old release.  No-op offline.
        """
        if self.offline:
            return
        live = self._get_json("status.json").get("chembl_db_version")
        if live != self.release:
            n = self.discard_staged()
            raise RuntimeError(
                f"ChEMBL release changed during the run: {self.release} -> {live} "
                f"({n} unverified cache entries discarded)"
            )
        self._promote()

    # ------------------------------------------------------------------ pagination
    def _paginate(self, endpoint: str, params: Mapping[str, Any]) -> list[dict]:
        records: list[dict] = []
        offset = 0
        total = None
        collection = None
        while True:
            page = self._get_json(
                f"{endpoint}.json", {**params, "limit": self.page_limit, "offset": offset}
            )
            if collection is None:
                collection = next(k for k in page if k != "page_meta")
            batch = page[collection]
            meta = page["page_meta"]
            if total is None:
                total = int(meta["total_count"])
            elif int(meta["total_count"]) != total:
                raise RuntimeError(f"{endpoint}: total_count changed during pagination")
            records.extend(batch)
            offset += len(batch)
            if not meta.get("next") or not batch:
                break
        if len(records) != total:
            raise RuntimeError(f"{endpoint}: got {len(records)} records, expected {total}")
        return records

    def fetch_all(self, endpoint: str, params: Mapping[str, Any]) -> list[dict]:
        """All records of a filtered list endpoint, fully paginated and cached.

        ``params`` are ChEMBL filter arguments (e.g. ``{"mesh_id__in": "D003924",
        "order_by": "drugind_id"}``).  Include an ``order_by`` on a unique field so the page
        boundaries are stable.
        """
        return self._cached(endpoint, params, lambda: self._paginate(endpoint, params))

    # ------------------------------------------------------------------ endpoints
    def drug_indications(self, mesh_ids: Iterable[str]) -> list[dict]:
        """All ``drug_indication`` rows whose ``mesh_id`` is in ``mesh_ids`` (sorted by drugind_id)."""
        ids = sorted({str(m).strip() for m in mesh_ids})
        out: list[dict] = []
        for mid in ids:  # one cache entry per descriptor so cohorts can share entries
            out.extend(self.fetch_all("drug_indication", {"mesh_id": mid, "order_by": "drugind_id"}))
        return sorted(out, key=lambda r: int(r["drugind_id"]))

    def drug_indications_icontains(self, text: str) -> list[dict]:
        """``drug_indication`` rows with ``mesh_heading__icontains=text`` (v0.x retrieval)."""
        return self.fetch_all(
            "drug_indication", {"mesh_heading__icontains": text, "order_by": "drugind_id"}
        )

    def mesh_headings(self) -> list[dict]:
        """Distinct ``(mesh_id, mesh_heading)`` pairs in ``drug_indication`` with row counts.

        Derived from a full scan of the endpoint (cached per release); sorted by mesh_id,
        then heading.
        """

        def produce() -> list[dict]:
            rows = self._paginate(
                "drug_indication", {"only": "drugind_id,mesh_id,mesh_heading", "order_by": "drugind_id"}
            )
            counts: dict[tuple[str, str], int] = {}
            for r in rows:
                k = (r.get("mesh_id") or "", r.get("mesh_heading") or "")
                counts[k] = counts.get(k, 0) + 1
            return [
                {"mesh_id": m, "mesh_heading": h, "n_rows": n}
                for (m, h), n in sorted(counts.items())
            ]

        return self._cached("drug_indication/_mesh_headings", {}, produce)

    @staticmethod
    def _project_molecule(rec: dict) -> dict:
        out = {k: rec.get(k) for k in MOLECULE_FIELDS}
        structs = rec.get("molecule_structures")
        out["molecule_structures"] = (
            {k: structs.get(k) for k in _STRUCTURE_FIELDS} if isinstance(structs, dict) else None
        )
        return out

    def molecules(self, chembl_ids: Iterable[str]) -> dict[str, dict | None]:
        """Molecule records by ChEMBL ID; ``None`` for IDs the API does not return.

        Each ID is cached individually (endpoint ``molecule``, params
        ``{"molecule_chembl_id": id}``), so overlapping cohorts reuse entries; missing IDs
        are fetched in sorted batches of ``molecule_batch``.
        """
        ids = sorted({str(i) for i in chembl_ids}, key=chembl_sort_key)
        out: dict[str, dict | None] = {}
        missing: list[str] = []
        for cid in ids:
            path, key = self._cache_path("molecule", {"molecule_chembl_id": cid})
            hit, data = self._load(path, key)
            if hit:
                out[cid] = data
            else:
                missing.append(cid)
        if missing and self.offline:
            raise CacheMiss(f"{len(missing)} molecules not cached for {self.release}: {missing[:5]}...")
        for i in range(0, len(missing), self.molecule_batch):
            batch = missing[i : i + self.molecule_batch]
            recs = self._paginate(
                "molecule",
                {
                    "molecule_chembl_id__in": ",".join(batch),
                    "only": ",".join(MOLECULE_FIELDS),
                    "order_by": "molecule_chembl_id",
                },
            )
            got = {r["molecule_chembl_id"]: self._project_molecule(r) for r in recs}
            for cid in batch:
                rec = got.get(cid)  # None is cached too: "the release has no such record"
                path, key = self._cache_path("molecule", {"molecule_chembl_id": cid})
                self._store(path, key, rec)
                out[cid] = rec
        return {cid: out[cid] for cid in ids}

    def molecules_by_atc(self, atc_prefix: str) -> list[dict]:
        """Molecules with an ATC level-5 code starting with ``atc_prefix`` (projected records)."""
        recs = self.fetch_all(
            "molecule",
            {
                "atc_classifications__level5__startswith": atc_prefix,
                "only": ",".join(MOLECULE_FIELDS),
                "order_by": "molecule_chembl_id",
            },
        )
        return [self._project_molecule(r) for r in recs]

    def atc_classes(self, atc_prefix: str) -> list[dict]:
        """``atc_class`` rows (level1-5 codes and descriptions) under ``atc_prefix``."""
        return self.fetch_all(
            "atc_class", {"level5__startswith": atc_prefix, "order_by": "level5"}
        )


class MeshClient(_HttpJson):
    """Optional NLM MeSH client: entry terms (synonyms) and tree descendants.

    Cache entries live under ``<cache_dir>/mesh/`` keyed by (endpoint, params).  MeSH has
    yearly versions; NLM serves the current year, which is recorded by the caller via
    :meth:`descendants` results being stored in the cohort JSON sidecar.
    """

    def __init__(
        self,
        cache_dir: str | os.PathLike | None = None,
        *,
        offline: bool = False,
        base_url: str = MESH_BASE_URL,
        session: requests.Session | None = None,
        timeout: float = 60.0,
        max_retries: int = 4,
        backoff: float = 1.0,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        super().__init__(base_url, session, timeout, max_retries, backoff, sleep, offline)
        self.cache_dir = Path(cache_dir) if cache_dir is not None else default_cache_dir() / "mesh"

    def _cached(self, endpoint: str, params: Mapping[str, Any], producer: Callable[[], Any]) -> Any:
        key = {"endpoint": endpoint, "params": dict(sorted(params.items()))}
        digest = hashlib.sha256(_canonical(key).encode()).hexdigest()[:24]
        path = self.cache_dir / endpoint.replace("/", "_") / f"{digest}.json.gz"
        if path.exists():
            blob = _read_gz_json(path)
            if blob.get("key") != key:
                raise RuntimeError(f"cache key mismatch in {path}")
            return blob["data"]
        if self.offline:
            raise CacheMiss(f"MeSH {endpoint} {params} not cached under {self.cache_dir}")
        data = producer()
        _write_gz_json(path, {"key": key, "data": data})
        return data

    def _sparql(self, query: str, page: int = 1000) -> list[dict]:
        """Run a SELECT query (which must have an ORDER BY), paging with limit/offset."""
        out: list[dict] = []
        offset = 0
        while True:
            res = self._get_json(
                "sparql", {"query": query, "format": "JSON", "limit": page, "offset": offset}
            )
            rows = [{k: v["value"] for k, v in b.items()} for b in res["results"]["bindings"]]
            out.extend(rows)
            if len(rows) < page:
                return out
            offset += page

    @staticmethod
    def _did(uri: str) -> str:
        return uri.rsplit("/", 1)[-1]

    def entry_terms(self, descriptor_id: str) -> list[str]:
        """All term labels (preferred + entry terms) of a descriptor, sorted."""

        def produce() -> list[str]:
            d = self._get_json("lookup/details", {"descriptor": descriptor_id})
            return sorted({t["label"] for t in d.get("terms", [])})

        return self._cached("lookup/details", {"descriptor": descriptor_id}, produce)

    def descriptors_for_text(self, text: str, limit: int = 50, match: str = "contains") -> dict[str, list[str]]:
        """Descriptors having a term that matches ``text`` -> sorted matched term labels.

        ``match`` is the NLM lookup mode: ``contains`` (default) or ``exact``
        (case-insensitive whole-label match).
        """
        if match not in ("contains", "exact"):
            raise ValueError(f"match must be 'contains' or 'exact', not {match!r}")
        text = " ".join(text.split())

        def produce() -> dict[str, list[str]]:
            terms = self._get_json(
                "lookup/term", {"label": text, "match": match, "limit": limit}
            )
            if not terms:
                return {}
            values = " ".join(f"<{t['resource']}>" for t in terms)
            rows = self._sparql(
                "PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#> "
                f"SELECT DISTINCT ?t ?d WHERE {{ VALUES ?t {{ {values} }} "
                "?d meshv:concept|meshv:preferredConcept ?c . "
                "?c meshv:term|meshv:preferredTerm ?t . } ORDER BY ?t ?d"
            )
            label = {t["resource"]: t["label"] for t in terms}
            out: dict[str, set] = {}
            for r in rows:
                out.setdefault(self._did(r["d"]), set()).add(label[r["t"]])
            return {k: sorted(v) for k, v in sorted(out.items())}

        return self._cached("lookup/term", {"label": text, "match": match, "limit": limit}, produce)

    def descendants(self, descriptor_id: str) -> list[dict]:
        """All descriptors below ``descriptor_id`` in any of its tree positions.

        SPARQL: ``?d meshv:treeNumber/meshv:parentTreeNumber+ ?pt`` with ``?pt`` a tree
        number of the root.  MeSH is a polyhierarchy, so descendants can include terms
        that are only partly about the root concept (e.g. *Carcinoma, Lobular* under
        *Breast Neoplasms*); the caller must show and confirm them.
        """

        def produce() -> list[dict]:
            rows = self._sparql(
                "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> "
                "PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#> "
                "PREFIX mesh: <http://id.nlm.nih.gov/mesh/> "
                f"SELECT DISTINCT ?d ?label WHERE {{ mesh:{descriptor_id} meshv:treeNumber ?pt . "
                "?t meshv:parentTreeNumber+ ?pt . ?d meshv:treeNumber ?t . ?d rdfs:label ?label } "
                "ORDER BY ?d"
            )
            return [{"mesh_id": self._did(r["d"]), "mesh_heading": r["label"]} for r in rows]

        return self._cached("sparql/descendants", {"descriptor": descriptor_id}, produce)
