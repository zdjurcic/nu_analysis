# caf_MC — 2×2 MC inputs

12 MiniRun6.5 RHC CAFs used for local verification (0000000-0000011, ~21 MB total, 2368 spills).
Not committed (see top-level `.gitignore`) — fetch with:

```bash
# From 2x2_mc/
./download_caf.sh
# or manual:
# xrdcp / curl from portal.nersc.gov as in original local_full/caf_MC/
ls -lh caf_MC/*.root
```

For CI, the verified PNGs + cutflow JSONs in `plots/` are committed, so plots render without inputs.
