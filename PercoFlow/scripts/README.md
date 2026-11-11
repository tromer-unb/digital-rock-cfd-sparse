# scripts/

| script | what it does | when to use it |
|---|---|---|
| `install_openfoam_userspace.sh [PREFIX]` | Downloads the Ubuntu 24.04 OpenFOAM v1912 packages with `apt-get download`, which needs no root, and extracts them under `PREFIX` (default `~/opt/openfoam-v1912`). It writes `PREFIX/env_openfoam.sh`, records the SHA-256 of every package in `PREFIX/debs/SHA256SUMS` and runs a smoke test. | No root access. With root, use `sudo apt install openfoam`. |
| `make_checksums.sh` | Regenerates `SHA256SUMS.txt`: the SHA-256 of every file tracked by git, except the checksum file itself. | Before every commit that is meant to be released (see docs/AUDIT_AND_TRACEABILITY.md). |

After installing OpenFOAM in user space, run `export PERCOFLOW_OPENFOAM_ENV=PREFIX/env_openfoam.sh` and then
`percoflow openfoam-check`.
