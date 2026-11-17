#!/bin/bash
set -e
# Fetch 12 MiniRun6.5 RHC CAFs for local_full verification
# Update BASE_URL if portal moves; files are 0000000-0000011.CAF.root (~1.8 MB each)
BASE_URL="https://portal.nersc.gov/project/dune/data/nd_prototypes/MiniRun6.5_1E19_RHC.caf"
mkdir -p caf_MC
for i in $(seq -w 0 11); do
  fname="MiniRun6.5_1E19_RHC.caf.00000${i}.CAF.root"
  # Try portal, fallback to local_full original path hint
  echo "Fetching $fname ..."
  curl -L -o "caf_MC/$fname" "$BASE_URL.$i.CAF.root" || echo "  -> missing at $BASE_URL.$i.CAF.root (check portal.nersc.gov)"
done
ls -lh caf_MC/*.root 2>&1 | head -20
