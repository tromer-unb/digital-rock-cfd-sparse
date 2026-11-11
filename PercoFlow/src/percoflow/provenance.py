"""Run manifests: what was run, on what, with which software, by whom and when.

Every run directory contains ``manifest.json`` with

* the complete configuration (``config.RunConfig``) actually used,
* SHA-256 hashes of every input file and of every output file written by the run,
* package version and git commit of the installed source (when available), Python and library versions,
* the OpenFOAM version used by the cross-check (or the explicit reason why it did not run),
* host name, start and end time (UTC), wall time, command line.

``verify_manifest(run_dir)`` recomputes the output hashes, so a result can be checked for tampering or corruption later.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from . import __version__

LIBRARIES = ("numpy", "scipy", "skfem", "skimage", "gmsh", "meshio", "matplotlib", "PIL", "yaml")


def sha256(path, chunk=1 << 20) -> str:
    """SHA-256 hex digest of a file, read in 1 MiB blocks."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _git_commit():
    """Commit of the installed source and whether the working tree had uncommitted changes (None outside git)."""
    here = Path(__file__).resolve().parent
    try:
        out = subprocess.run(["git", "-C", str(here), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=5)
        dirty = subprocess.run(["git", "-C", str(here), "status", "--porcelain"], capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            return {"commit": out.stdout.strip(), "dirty": bool(dirty.stdout.strip())}
    except (OSError, subprocess.SubprocessError):
        pass
    return {"commit": None, "dirty": None}


def environment() -> dict:
    """Software and machine identity recorded in every manifest."""
    versions = {}
    for name in LIBRARIES:
        try:
            mod = __import__(name)
            versions[name] = getattr(mod, "__version__", getattr(mod, "GMSH_API_VERSION", "unknown"))
        except Exception:      # noqa: BLE001 - a missing optional library is recorded, not fatal
            versions[name] = None
    return {"percoflow": __version__, "git": _git_commit(), "python": sys.version.split()[0], "platform": platform.platform(),
            "host": socket.gethostname(), "libraries": versions,
            "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}}


class Manifest:
    """Collects provenance during a run and writes it at the end (``with Manifest(run_dir, ...) as m:``)."""

    def __init__(self, run_dir, kind: str, config: dict, inputs=()):
        self.run_dir = Path(run_dir)
        self.data = {"kind": kind, "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
                     "command": " ".join(sys.argv), "config": config, "environment": environment(),
                     "inputs": {str(p): sha256(p) for p in inputs if Path(p).is_file()}, "notes": []}
        self._t0 = time.time()

    def note(self, text: str):
        """Add a free-text remark (e.g. "OpenFOAM cross-check DISABLED by the user") to the manifest."""
        self.data["notes"].append(text)

    def __enter__(self):
        self.run_dir.mkdir(parents=True, exist_ok=True)
        return self

    def __exit__(self, exc_type, exc, tb):
        self.data["status"] = "failed" if exc_type else self.data.get("status_override", "completed")
        if exc_type:
            self.data["error"] = f"{exc_type.__name__}: {exc}"
        self.data["finished_utc"] = datetime.now(timezone.utc).isoformat()
        self.data["wall_time_s"] = time.time() - self._t0
        outputs = {}
        for p in sorted(self.run_dir.rglob("*")):
            if p.is_file() and p.name != "manifest.json":
                outputs[str(p.relative_to(self.run_dir))] = sha256(p)
        self.data["outputs"] = outputs
        (self.run_dir / "manifest.json").write_text(json.dumps(self.data, indent=1, default=str))
        return False


def verify_manifest(run_dir) -> dict:
    """Recompute the hashes of the outputs listed in ``manifest.json``; returns {"ok": bool, "mismatch": [...]}."""
    run_dir = Path(run_dir)
    data = json.loads((run_dir / "manifest.json").read_text())
    bad = [rel for rel, h in data.get("outputs", {}).items() if not (run_dir / rel).is_file() or sha256(run_dir / rel) != h]
    return {"ok": not bad, "mismatch": bad, "n_files": len(data.get("outputs", {}))}
