# nu_analysis

DUNE ND 2×2 CAFs — reproducible Python analyses of MiniRun6.5 RHC (and future data).

## Layout

```
nu_analysis/
  2x2_mc/   — MC verified analysis (MiniRun6.5 RHC, 12-file slice, 21 MB)
            mirrors nominal_mc.cc : fiducial+track+MINERvA → hadron_mult / length / angular
            see 2x2_mc/README.md and 2x2_mc/analysis_summary.html
  2x2_data/ — placeholder for data (same selections, data/MC overlay)
```

## Quick start

```bash
# MC verified (tight MINERvA = nominal_mc.cc L2625-2785)
conda activate dune_root   # ROOT 6.40.04
python 2x2_mc/analyze_verified_minerva.py   # ~15s for 12 files → 2x2_mc/plots_verified/
open 2x2_mc/analysis_summary.html
```

Each `2x2_*` folder is self-contained with its own `deps/duneanaobj/`, `nominal_mc.cc`, and `README.md`.

## Provenance

MC work developed and verified 2026-11-16 in `Service_Enabled_Science/05_Agentic_Workflows/caf_analysis_code/local_full/`
(2368 spills, 868 fiducial, 349 track, loose 166 matched / tight 4 matched).
`2x2_mc/` is a clean copy of `local_full/` minus ROOT inputs (see `2x2_mc/caf_MC/README.md`).

## License

MIT — see LICENSE
