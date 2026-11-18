#!/bin/bash
set -euo pipefail
# Download 1000 CAFs from NERSC portal to Eagle on Polaris compute node.
# Proxy already exported by job_preamble (http_proxy=http://proxy.alcf.anl.gov:3128).
EAGLE_BASE="/lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/caf_MC"
SUBDIR="0000000"
DEST="$EAGLE_BASE/$SUBDIR"
BASE_URL="https://portal.nersc.gov/project/dune/data/2x2/simulation/productions/MiniRun6.5_1E19_RHC/MiniRun6.5_1E19_RHC.caf/CAF/$SUBDIR"
mkdir -p "$DEST"
echo "Downloading 1000 CAFs to $DEST from $BASE_URL"
echo "Proxy: $http_proxy"
# Portal files are MiniRun6.5_1E19_RHC.caf.<7-digit>.CAF.root, 0000000-0000999
# We fetch in a loop with resume and check; on Polaris, curl respects http_proxy automatically.
for i in $(seq -w 0 999); do
  fname="MiniRun6.5_1E19_RHC.caf.0000${i}.CAF.root"
  # Some portals use 0000000/ prefix already in BASE_URL; adjust if 404:
  url="$BASE_URL/$fname"
  dest="$DEST/$fname"
  if [ -f "$dest" ] && [ $(stat -c%s "$dest" 2>/dev/null || stat -f%z "$dest" 2>/dev/null) -gt 100000 ]; then
    echo "[$i/999] exists $fname"
    continue
  fi
  echo "[$i/999] curl $fname"
  if ! curl -L --fail --retry 2 --connect-timeout 15 -o "$dest.tmp" "$url"; then
    echo "  WARN: $url failed, trying alternate listing..."
    # Fallback: try bare CAF/ without SUBDIR duplication (portal layout varies)
    alt_url="https://portal.nersc.gov/project/dune/data/2x2/simulation/productions/MiniRun6.5_1E19_RHC/MiniRun6.5_1E19_RHC.caf/CAF/$fname"
    curl -L --fail --retry 2 -o "$dest.tmp" "$alt_url" || { echo "  ERROR: $fname not found at either URL"; rm -f "$dest.tmp"; continue; }
  fi
  mv "$dest.tmp" "$dest"
  # throttle every 100 to be nice to portal
  if [ $((10#$i % 100)) -eq 0 ]; then sleep 1; fi
done
echo "=== Done ==="
ls "$DEST" | wc -l
du -sh "$EAGLE_BASE"
ls "$DEST" | head -5
ls "$DEST" | tail -5
