#!/usr/bin/env bash
# Install OpenFOAM v1912 WITHOUT root privileges by extracting the Ubuntu 24.04 (noble) packages into a user prefix.
# This reproduces the reference installation used to validate PercoFlow 1.0.0 (Ubuntu package openfoam 1912.200626-2build3).
#
#   bash scripts/install_openfoam_userspace.sh [PREFIX]        (default PREFIX = ~/opt/openfoam-v1912)
#   source PREFIX/env_openfoam.sh && simpleFoam -help
#   export PERCOFLOW_OPENFOAM_ENV=PREFIX/env_openfoam.sh      (so that percoflow finds it)
#
# With root on Ubuntu/Debian the simpler route is:  sudo apt install openfoam
# Requirements: apt-get (to download the .deb files), dpkg-deb, an OpenMPI runtime (libopenmpi3 / openmpi-bin) on the
# system for parallel runs. Serial runs do not need mpirun. Every downloaded package is hashed into PREFIX/debs/SHA256SUMS.
set -euo pipefail
PREFIX="${1:-$HOME/opt/openfoam-v1912}"
PKGS=(openfoam libopenfoam libscotch-7.0 libptscotch-7.0 libboost-system1.83.0)
mkdir -p "$PREFIX/debs"
cd "$PREFIX/debs"
echo "Downloading ${PKGS[*]} (apt-get download, no root needed)"
apt-get download "${PKGS[@]}"
sha256sum ./*.deb > SHA256SUMS
for d in ./*.deb; do dpkg-deb -x "$d" "$PREFIX"; done
cat > "$PREFIX/env_openfoam.sh" <<EOF
# OpenFOAM v1912 extracted in user space by percoflow/scripts/install_openfoam_userspace.sh - source this file
export OF_ROOT="$PREFIX"
export PATH="\$OF_ROOT/usr/bin:\$PATH"
export LD_LIBRARY_PATH="\$OF_ROOT/usr/lib:\$OF_ROOT/usr/lib/x86_64-linux-gnu:\$OF_ROOT/usr/lib/openmpi-system:\${LD_LIBRARY_PATH:-}"
export WM_PROJECT=OpenFOAM WM_PROJECT_VERSION=v1912
export WM_PROJECT_DIR="\$OF_ROOT/usr/share/openfoam"
export FOAM_ETC="\$OF_ROOT/usr/share/openfoam/etc"
export FOAM_LIBBIN="\$OF_ROOT/usr/lib" FOAM_APPBIN="\$OF_ROOT/usr/bin"
export WM_MPLIB=SYSTEMOPENMPI WM_COMPILER=Gcc WM_PRECISION_OPTION=DP WM_LABEL_SIZE=32
EOF
# smoke test
# shellcheck disable=SC1091
source "$PREFIX/env_openfoam.sh"
if simpleFoam -help >/dev/null 2>&1; then
  echo "OK: simpleFoam works. Now:  export PERCOFLOW_OPENFOAM_ENV=$PREFIX/env_openfoam.sh && percoflow openfoam-check"
else
  echo "simpleFoam does not start; check missing shared libraries with:  ldd $PREFIX/usr/bin/simpleFoam | grep 'not found'" >&2
  exit 1
fi
