# 2×2 MC — Polaris Full (1000 files)

Polaris-native verified analysis: 50-file test → 1000-file production, all 13 plots, same tight cuts as `2x2_mc/` but built for Polaris ROOT 6.40.04.

## What is here (validated 2026-11-18)
- **Env:** `WORK=/home/zdjurcic/2x2/root_env` with `ROOT 6.40.04` + `matplotlib` (conda-forge, 676 MB) — see `POLARIS_GUIDE.html` §2
- **Lib:** `deps/duneanaobj/install/lib/libduneanaobj_StandardRecord.so` (1.1M, `rootcling + g++` for Polaris, no `cetmodules`)
- **Scripts:** `download_cafs.sh` (portal → Eagle 1000 CAFs, 1.8 GB), `run_50_full.py` (50 files, 13 plots), `run_1000_full.py` (1000 files, 13 plots), `run_2x2_mc_polaris.py`
- **Plots:** `plots_50_full/` (50 files: 9925 spills → 23 matched) + `plots_full_1000/` (1000 files: 198460 spills → 499 matched) — on Polaris at `/lus/eagle/.../plots_*_full/`, tarred to `$HOME/plots_*.tar.gz` and mirrored here
- **Guide:** `POLARIS_GUIDE.html` (self-run instructions, no assistant) + `PROMPTS.md` (11 prompts to recreate from scratch)

## Run yourself on Polaris (no assistant)
```bash
# Every new shell/job:
source $WORK/polaris_setup/env_polaris.sh
conda activate $WORK/root_env
python -u $HOME/polaris_full/run_50_full.py    # 50 files, 30-min debug, ~1 min → plots_50_full/ on Eagle
python -u $HOME/polaris_full/run_1000_full.py  # 1000 files, 60-min debug, ~5 min → plots_full_1000/ on Eagle
ls /lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/caf_MC/plots_50_full/
cat /lus/eagle/.../plots_50_full/cutflow_50_full.json
```
See `POLARIS_GUIDE.html` for full `qsub` / IRI `submit_job` examples and how to bring plots back via Globus (`home` → laptop).

## Provenance
- `2x2_mc/` = local 12-file verification (laptop, `ROOT 6.32`, `libduneanaobj_full.so`)
- `2x2_mc_polaris/` = Polaris 50/1000-file production (Polaris, `ROOT 6.40.04`, `libduneanaobj_StandardRecord.so` built via `rootcling`)
- Both use identical `fiducial+track+MINERvA_tight` (L2461, L2573/2598, L2625-2785) — see `local_full/README.md` for cut table.

## To add 2x2_data later
Copy `run_1000_full.py` to `2x2_data/`, set `mnvOffsetX=-10 Y=5` for data, point `EAGLE` to data CAFs.
