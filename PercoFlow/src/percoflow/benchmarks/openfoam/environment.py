"""Locate an OpenFOAM installation and run its utilities.

Search order (first hit wins):
    1. the explicit ``env_script`` argument (OpenFOAMOptions.env_script),
    2. the environment variable PERCOFLOW_OPENFOAM_ENV,
    3. ``simpleFoam`` already on PATH (e.g. Debian/Ubuntu package, or an environment sourced beforehand),
    4. standard locations: /usr/lib/openfoam/openfoam*/etc/bashrc, /opt/openfoam*/etc/bashrc,
       /usr/share/openfoam/etc/bashrc, ~/OpenFOAM/OpenFOAM-*/etc/bashrc, ~/opt/openfoam-sysroot/env_openfoam.sh.

If nothing is found ``OpenFOAMNotFound`` is raised with installation instructions (docs/OPENFOAM.md). The workflow does not
silently skip the cross-check: switching it off requires ``--no-openfoam`` and is recorded in the run manifest.
"""
from __future__ import annotations

import glob
import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

INSTALL_HINT = (
    "OpenFOAM was not found. The OpenFOAM cross-check is part of every percoflow run.\n"
    "  Ubuntu/Debian:  sudo apt install openfoam            (reference version: OpenFOAM v1912)\n"
    "  no root:        bash scripts/install_openfoam_userspace.sh   (extracts the Ubuntu packages under ~/opt)\n"
    "  other systems:  https://www.openfoam.com/download\n"
    "Then point percoflow to it with  export PERCOFLOW_OPENFOAM_ENV=/path/to/etc/bashrc  (or --openfoam-env).\n"
    "To run without the cross-check (recorded as such in manifest.json) use --no-openfoam."
)


class OpenFOAMNotFound(RuntimeError):
    """No usable OpenFOAM installation was found (the message contains installation instructions)."""


@dataclass
class OpenFOAM:
    """A located OpenFOAM installation: the script that sets it up (if any) and its version string."""
    env_script: str | None       # None when the utilities are already on PATH
    version: str

    def command(self, args, cwd, log, nprocs: int = 1, timeout: float | None = None):
        """Run an OpenFOAM utility in ``cwd``; stdout and stderr go to ``cwd/log``. Raises on a non-zero exit code."""
        cmd = " ".join(shlex.quote(a) for a in args)
        if nprocs > 1:
            cmd = f"mpirun --oversubscribe -np {int(nprocs)} {cmd} -parallel"
        if self.env_script:
            cmd = f"source {shlex.quote(self.env_script)} >/dev/null 2>&1; {cmd}"
        with open(Path(cwd) / log, "w") as fh:
            r = subprocess.run(["bash", "-c", cmd], cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, timeout=timeout)
        if r.returncode != 0:
            raise RuntimeError(f"OpenFOAM command failed ({args[0]}), see {Path(cwd) / log}")

    def describe(self) -> dict:
        return {"env_script": self.env_script, "version": self.version}


def _probe(script: str | None) -> str | None:
    """Return the OpenFOAM version if simpleFoam can be found with this environment, else None."""
    pre = f"source {shlex.quote(script)} >/dev/null 2>&1; " if script else ""
    cmd = pre + 'command -v simpleFoam >/dev/null && echo "${WM_PROJECT_VERSION:-unknown}"'
    try:
        r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    out = r.stdout.strip().splitlines()
    return out[-1] if r.returncode == 0 and out else None


def find_openfoam(env_script: str | None = None) -> OpenFOAM:
    """Locate OpenFOAM following the search order of the module docstring; raise OpenFOAMNotFound otherwise."""
    candidates = []
    if env_script:
        candidates.append(os.path.expanduser(env_script))
    if os.environ.get("PERCOFLOW_OPENFOAM_ENV"):
        candidates.append(os.path.expanduser(os.environ["PERCOFLOW_OPENFOAM_ENV"]))
    if env_script or os.environ.get("PERCOFLOW_OPENFOAM_ENV"):
        for c in candidates:
            v = _probe(c)
            if v:
                return OpenFOAM(c, v)
        raise OpenFOAMNotFound(f"the given OpenFOAM environment {candidates} does not provide simpleFoam.\n" + INSTALL_HINT)
    if shutil.which("simpleFoam"):
        v = _probe(None)
        if v:
            return OpenFOAM(None, v)
    patterns = ["/usr/lib/openfoam/openfoam*/etc/bashrc", "/opt/openfoam*/etc/bashrc", "/usr/share/openfoam/etc/bashrc",
                os.path.expanduser("~/OpenFOAM/OpenFOAM-*/etc/bashrc"), os.path.expanduser("~/opt/openfoam-sysroot/env_openfoam.sh")]
    for pat in patterns:
        for c in sorted(glob.glob(pat), reverse=True):
            v = _probe(c)
            if v:
                return OpenFOAM(c, v)
    raise OpenFOAMNotFound(INSTALL_HINT)
