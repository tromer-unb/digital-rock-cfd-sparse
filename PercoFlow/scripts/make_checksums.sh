#!/usr/bin/env bash
# Regenerate SHA256SUMS.txt: SHA-256 of every file tracked by git (the checksum file itself excluded).
# Check the repository later with:  sha256sum -c SHA256SUMS.txt
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
git ls-files -z | grep -zv '^SHA256SUMS.txt$' | sort -z | xargs -0 sha256sum > SHA256SUMS.txt
echo "wrote SHA256SUMS.txt ($(wc -l < SHA256SUMS.txt) files)"
