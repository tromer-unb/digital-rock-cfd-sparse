"""Read-only access to the author's pipeline (the code under audit). Nothing in the legacy file is modified."""
import hashlib, importlib.util, sys
from pathlib import Path

LEGACY_MD5 = "57afdb1f"          # prefix of the md5 of the audited file (full md5 in docs/AUDIT-TRAIL.md and MANIFEST.sha256)


def pipeline_path(root=None) -> Path:
    here = Path(__file__).resolve()
    cands = [Path(root) / "digital_rock_cfd_pipeline_sparse.py"] if root else []
    cands += [here.parents[2] / "legacy" / "digital_rock_cfd_pipeline_sparse.py"]
    for c in cands:
        if c.exists():
            return c
    raise FileNotFoundError("digital_rock_cfd_pipeline_sparse.py not found (looked in %s)" % cands)


def load_pipeline(root=None, check_md5=True):
    p = pipeline_path(root)
    if check_md5:
        h = hashlib.md5(p.read_bytes()).hexdigest()
        if not h.startswith(LEGACY_MD5):
            print(f"WARNING: legacy pipeline md5 {h[:8]} differs from the audited file ({LEGACY_MD5})", file=sys.stderr)
    spec = importlib.util.spec_from_file_location("drc_legacy", p)
    m = importlib.util.module_from_spec(spec); sys.modules["drc_legacy"] = m; spec.loader.exec_module(m)
    return m
