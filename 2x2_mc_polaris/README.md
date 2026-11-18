# polaris_full — 1000-file 2x2 MC on Polaris (Option 1)

This is the Polaris job deck for the full MiniRun6.5 stats. It pulls 1000 CAFs from the NERSC portal
on the compute node (via the http_proxy the preamble sets) into your Eagle, then runs the verified
analysis (same cuts as local_full/).

## What it does
1. `download_cafs.sh` — `curl` loop from https://portal.nersc.gov/project/dune/data/2x2/simulation/productions/MiniRun6.5_1E19_RHC/MiniRun6.5_1E19_RHC.caf/CAF/0000000/
   → `/lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/caf_MC/0000000/` (≈1.75 GB, 1000 × ~1.8 MB)
2. `run_2x2_mc_polaris.py` — same fiducial+track+MINERvA_tight as local_full/analyze_verified_minerva.py, but
   pointed at Eagle and writing plots to Eagle (`plots_full_1000/`), not local.

## Run locally first (stage 1 check)
```bash
# From 05_Agentic_Workflows/:
cat caf_analysis_code/polaris_full/download_cafs.sh
cat caf_analysis_code/polaris_full/run_2x2_mc_polaris.py
```

## Stage to Polaris (stage 2)
Eagle writes must be `/lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/<you>/` with filename — the workshop guard enforces this.
Only the *deck* is staged via Globus (3 files, <1 MB); the 1000 CAFs are fetched *on* Polaris, so you stay under MAX_STAGE_FILES=500 / 1 GiB.

```
local_endpoint → globus_ls(alcf, /lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/)
stage_to_alcf(local_path=.../polaris_full/download_cafs.sh, alcf_path=/lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/polaris_full/download_cafs.sh)
stage_to_alcf(.../run_2x2_mc_polaris.py, .../polaris_full/run_2x2_mc_polaris.py)
stage_to_alcf(.../deps, .../polaris_full/deps, recursive=True)  # libduneanaobj_full.so
transfer_status(task_id) → SUCCEEDED
globus_ls to verify
```

## Submit (stage 3)
```
get_system_status("polaris") → check queues (debug vs reservation)
submit_job(system="polaris", commands="bash /lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/polaris_full/download_cafs.sh && python /lus/eagle/.../run_2x2_mc_polaris.py",
           stdout_path="/home/zdjurcic/polaris_1000.out", nodes=1, queue="debug", account="Nu_Novel", walltime_sec=1800)
```

## Monitor / report (stages 4-5)
```
get_job_state("polaris", job_id) polled 10s→30s
read_file("/home/zdjurcic/polaris_1000.out") + globus_ls plots
```

Expected cutflow (1000 ≈ 83×12): track ~29k, loose matched ~13.8k (11.2k CC), tight ~330 (all CC).

## Changing Event Selection

Edit the same 6 cuts at the top of `run_50_full.py` / `run_1000_full.py` / `run_2x2_mc_polaris.py` (mirror of `local_full/analyze_verified*.py` — keep in sync via `../../nu_analysis/sync_local_full.sh to_nu`):

| Selection | Line | Example |
|-----------|------|---------|
| **Fiducial** `5< |x|<59, 5< |z|<59.5, |y|<57` | `if not (abs(vx)>5 and abs(vx)<59 ...):` | `5→10` tighter |
| **Track** `primary && |start-vtx|<5 && len>3` | `minTrkLength=3.0` + `diff>5` + `abs(pdg) in (2212,13,211,321)` | `3.0→5.0` cm |
| **Hadron mult** `p/π/K` | `if abs(pdg) in (2212,211,321):` | drop `321` for p/π only |
| **MINERvA tight** `maxDot>0.99, |Δ|<15, |Δatan|<0.06` | `passes_minerva_tight()` | `0.99→0.995`, `15→10` |
| **MC vs Data offset** | `mnvOffsetX=0; mnvOffsetY=0` | `-10,5` for data (L1989) |
| **Truth match** `|Δvtx|<5` | `if abs(vx-tvx)<5` | `5→2` cm |

After edit on laptop:
```bash
# Stage updated script to Polaris (30-min debug, Nu_Novel):
# stage_to_alcf(local_path=".../polaris_full/run_50_full.py", alcf_path="/home/zdjurcic/polaris_full/run_50_full.py")
# submit_job(commands="source $WORK/polaris_setup/env_polaris.sh && conda activate $WORK/root_env && python -u $HOME/polaris_full/run_50_full.py", stdout_path="/home/zdjurcic/run_50.out", nodes=1, queue="debug", account="Nu_Novel", walltime_sec=1800)
# Check: cat /lus/eagle/.../plots_50_full/cutflow_50_full.json
```

This README is the Polaris complement to `local_full/README.md` — both now have the same selection table. See `POLARIS_GUIDE.html` §4 and `PROMPTS.md` #9.
