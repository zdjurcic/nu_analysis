# Prompts to Recreate the Polaris 2×2 Analysis from Scratch

Copy-paste these in order to a fresh chat. Each assumes you are in `05_Agentic_Workflows/` and have `caf_analysis_code/local_full/` with `nominal_mc.cc` and `deps/duneanaobj/src`.

### 0. Bootstrap local verification (12 files)
```
I'm in 05_Agentic_Workflows/caf_analysis_code/local_full/ with nominal_mc.cc and deps/duneanaobj/src.
Create README.md that documents the 6 cuts from nominal_mc.cc:
fiducial (L2461), track (L2573/2598 primary+|start-vtx|<5 len>3), MINERvA (L2625 tight dot>0.99 trans<15 ang<0.06),
truth match (truthOverlap argmax + |Δvtx|<5), and the 3 physics groups:
hadron_mult (p/π/K), track_length (|end-start| μ/p/π), angular (cosθ=pz/|p| μ/p/π).
Point to caf_MC/ (12 MiniRun6.5 files 0000000-0000011, 21 MB) and deps/duneanaobj.
```

### 1. Reproduce the uproot blocker
```
In local_full/caf_MC/ with 12 files, run with ../.venv/bin/python:
  br = t['rec/common/common.ixn.dlp/common.ixn.dlp.part.dlp']
  br.array(entry_start=0, entry_stop=1)
Capture the NotImplementedError: memberwise serialization of AsArray in caf::SRVector3D.
Then run f.file.show_streamers('caf::SRRecoParticle') and show TStreamerInfo for SRVector3D v10 vs AsArray.
```

### 2. Prove ROOT can read it (local)
```
ROOT is at /opt/homebrew/bin/root (6.32.08) and dune_root env at /Users/zdjurcic/miniconda3/envs/dune_root/bin/root (6.40.04).
Fix deps/duneanaobj/LinkDef.h (missing #endif), then:
  rootcling -f G__StandardRecord.cxx -c -p -I src src/duneanaobj/StandardRecord/StandardRecord.h LinkDef.h
  g++ -fPIC -shared -o libduneanaobj.so G__StandardRecord.cxx src/.../*.cxx -I src $(root-config --cflags --libs)
Load it in ROOT and dump entry 1: rec.common.ixn.dlp[0].part.dlp size (expect 12), with pdg/primary/E/start/end/p.
```

### 3. Build verified loose analyzer (local, 12 files)
```
Write local_full/analyze_verified.py that runs with dune_root's python and uses ROOT.caf.StandardRecord via libduneanaobj.so.
Loop 12 files, apply fiducial 5<|x|<59 etc., track primary+|start-vtx|<5 len>3, MINERvA loose (exiting |z|>58 proxy),
truthOverlap argmax + vertex<5, and fill Enu, vertex X/Y/Z, mode, nprim, hadron_mult (p/π/K), length μ/p/π, angular μ/p/π truth vs reco.
Save to plots_verified/ + cutflow_verified.json (expect 868 fiducial → 349 track → 166 matched).
```

### 4. Add tight MINERvA (local)
```
Extend to analyze_verified_minerva.py: build libduneanaobj_full.so with LinkDef_full.h (+SRTrack, SRMINERvA, SRTMS...).
Implement passes_minerva_tight(ixn, rec.nd.minerva) replicating L2625-2785: for each exiting reco track, loop nd.minerva.ixn tracks,
lenMin>10, dot>prev && |ΔY|<15 && |ΔX|<15 && |Δatan|<0.06, pass maxDot>0.99 (mcOnly → mnvOffset 0,0).
Gate on it (requireMinerva=True). Expect loose 166 matched → tight 4 matched on 12 files.
```

### 5. Create Polaris deck (Option 1, HTTP portal)
```
My Polaris user is zdjurcic, project Nu_Novel, staging area /lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/caf_MC,
home /home/zdjurcic, portal https://portal.nersc.gov/project/dune/data/2x2/.../CAF/0000000/.
Create caf_analysis_code/polaris_full/ with download_cafs.sh (curl loop 0-999 via proxy to $EAGLE/0000000/) and run_2x2_mc_polaris.py
(same cuts as analyze_verified_minerva.py but EAGLE_CAF=/lus/eagle/.../0000000 and outdir=plots_full_1000/).
Show the 3-file deck before staging (skill stage 1).
```

### 6. Stage deck to Polaris (Globus)
```
Polaris is up (get_system_status). Stage the deck via Globus (MAX_STAGE_FILES=500/1 GiB, so only the deck, not the 1000 CAFs):
stage_to_alcf(local_path=.../polaris_full/download_cafs.sh, alcf_path=/home/zdjurcic/polaris_full/download_cafs.sh) etc. for run_*.py + deps/*.so,
poll transfer_status → SUCCEEDED, globus_ls verify, then get_system_status → submit_job.
```

### 7. Setup conda/ROOT on Polaris (one-time)
```
Polaris has no ROOT module. On Polaris $WORK=/home/zdjurcic/2x2:
  module load conda; export CONDA_PKGS_DIRS=$WORK/conda_pkgs
  conda config --set proxy_servers.http $http_proxy
  conda create -p $WORK/root_env -c conda-forge root python=3.11 matplotlib -y  # 676 MB, ROOT 6.40.04
  conda activate $WORK/root_env && python -c "import ROOT; print(ROOT.gROOT.GetVersion())"
  conda clean --all -y
Every new job: source $WORK/polaris_setup/env_polaris.sh && conda activate $WORK/root_env
```

### 8. Build duneanaobj for Polaris ROOT (no cetmodules)
```
On Polaris with $WORK/root_env active:
  module load cmake/3.27.7
  cd /tmp && cp $HOME/polaris_full/deps/duneanaobj/LinkDef_full.h ./LinkDef.h; cp -r $WORK/deps/duneanaobj/src ./src
  rootcling -f G__StandardRecord.cxx -c -p -I src src/duneanaobj/StandardRecord/StandardRecord.h LinkDef.h
  g++ -fPIC -shared -o libtest.so G__StandardRecord.cxx src/duneanaobj/StandardRecord/*.cxx -I src $(root-config --cflags --libs)
  ls -lh libtest.so  # → 1.1M proves Option 1 works without cmake/cetmodules
Then for production:
  mkdir -p $WORK/deps/duneanaobj/build $WORK/deps/duneanaobj/install
  cd $WORK/deps/duneanaobj/build && cmake $WORK/deps/duneanaobj/src -DCMAKE_INSTALL_PREFIX=$WORK/deps/duneanaobj/install -DCMAKE_CXX_STANDARD=17 && make -j4 && make install
  # → $WORK/deps/duneanaobj/install/lib/libduneanaobj_StandardRecord.so for Polaris ROOT 6.40.04
If cmake fails with cetmodules, use the rootcling direct path above — it's faster for Python.
```

### 9. Submit 50-file tight test (verify env)
```
Submit a 30-min debug job (walltime 1800, queue debug, account Nu_Novel, nodes 1):
  source $WORK/polaris_setup/env_polaris.sh && conda activate $WORK/root_env
  python -u $HOME/polaris_full/run_50_tight_fixed.py  # 50 files, tight cuts → plots_50_tight/ on Eagle
Expect: 50 files → 9925 spills → 3670 fiducial → 1482 track → 26 tight → 23 matched (22 CC). If lib mismatch (Failed to get attribute caf), rebuild lib for Polaris ROOT as in step 8.
```

### 10. Submit 1000-file full (production)
```
Same 30-min debug (or 60-min for queue), but run_1000_full.py with files[:1000]:
  source $WORK/polaris_setup/env_polaris.sh && conda activate $WORK/root_env
  python -u $HOME/polaris_full/run_1000_full.py  # → /lus/eagle/.../plots_full_1000/ + cutflow_full1000.json (13 PNGs)
Expect: 1000 files → 198460 spills → 72951 fiducial → 29267 track → 532 tight → 499 matched (488 CC/11 NC) — ~5 min compute.
Then tar on Polaris: tar -czf $HOME/plots_full_1000.tar.gz -C /lus/eagle/.../caf_MC plots_full_1000
and Globus transfer home → local laptop (home collection 6115bd2c-957a-4543-abff-5fae52992ff2 → local c83fd5dc-deeb-11f1-9c8b-02ffe3faed53).
```

Run these 11 prompts in order — you get the same 12-file local verified + 50-file Polaris verified + 1000-file production with all 13 plots, and a reproducible nu_analysis/2x2_mc_polaris/ for GitHub.
