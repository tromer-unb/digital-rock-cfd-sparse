#!/usr/bin/env bash
# Build the Supplementary Data archive of the paper from a commit of this repository.
# The archive holds exactly the files tracked by git at that commit (git archive), under the folder
# digital-rock-cfd-sparse/, so it has the same content as the GitHub repository at that commit.
# Usage:  bash tools/make_supplementary_zip.sh [commit-or-tag] [output-directory]     (defaults: HEAD, dist/)
# Writes <output>/digital-rock-cfd-sparse-<commit>.zip and <output>/digital-rock-cfd-sparse-<commit>.zip.sha256
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
REF="${1:-HEAD}"
OUT="${2:-dist}"
if [ -n "$(git status --porcelain)" ]; then
  echo "warning: the working tree has uncommitted changes; they are NOT in the archive (it is built from $REF)" >&2
fi
SHORT="$(git rev-parse --short=7 "$REF")"
NAME="digital-rock-cfd-sparse-${SHORT}"
mkdir -p "$OUT"
git archive --format=zip -9 --prefix=digital-rock-cfd-sparse/ -o "$OUT/$NAME.zip" "$REF"
(cd "$OUT" && sha256sum "$NAME.zip" > "$NAME.zip.sha256")
echo "wrote $OUT/$NAME.zip ($(du -h "$OUT/$NAME.zip" | cut -f1)), commit $(git rev-parse "$REF")"
cat "$OUT/$NAME.zip.sha256"
