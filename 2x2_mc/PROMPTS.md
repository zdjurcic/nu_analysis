# Prompts to Recreate local_full from Scratch (No Assistant Needed)

Use these 11 prompts in order in a fresh chat. Each assumes you are in `05_Agentic_Workflows/` and have `caf_analysis_code/local_full/` empty except for `caf_MC/` (12 MiniRun6.5 RHC files), `nominal_mc.cc`, and `deps/duneanaobj/src/` (from `https://github.com/DUNE/duneanaobj.git -b v03_06_01`).

### 0. Bootstrap the doc
```
I'm in 05_Agentic_Workflows/caf_analysis_code/local_full/ with nominal_mc.cc and deps/duneanaobj/src.
Create README.md that documents the 6 cuts from nominal_mc.cc:
fiducial (L2461 5< |x|<59, 5< |z|<59.5, |y|<57), track (L2573/2598 primary+|start-vtx|<5 len>3),
MINERvA (L2625 tight dot>0.99 trans<15 ang<0.06), truth match (truthOverlap argmax + |Δvtx|<5),
and the 3 physics groups: hadron_mult (p/π/K), track_length (|end-start| μ/p/π), angular (cosθ=pz/|p| μ/p/π).
Point to caf_MC/ (12 files 0000000-0000011, 21 MB) and deps/duneanaobj.
```

### 1. Reproduce the uproot blocker (confirm the need for ROOT)
```
In local_full/caf_MC/ with 12 files, run with ../.venv/bin/python:
  br = t['rec/common/common.ixn.dlp/common.ixn.dlp.part.dlp']
  br.array(entry_start=0, entry_stop=1)
Capture the full traceback: NotImplementedError: memberwise serialization of AsArray in caf::SRVector3D.
Then run f.file.show_streamers('caf::SRRecoParticle') and f.file.show_streamers('caf::SRVector3D') to show
SRVector3D v10 (x,y,z) vs AsArray memberwise.
```

### 2. Prove ROOT can read it (set up ROOT + duneanaobj)
```
ROOT is at /opt/homebrew/bin/root (6.32.08) and dune_root env at /Users/zdjurcic/miniconda3/envs/dune_root/bin/root (6.40.04).
Fix deps/duneanaobj/LinkDef.h (add missing #endif), then:
  rootcling -f G__StandardRecord.cxx -c -p -I src src/duneanaobj/StandardRecord/StandardRecord.h LinkDef.h
  g++ -fPIC -shared -o libduneanaobj.so G__StandardRecord.cxx src/duneanaobj/StandardRecord/*.cxx -I src $(root-config --cflags --libs)
  ls -lh libduneanaobj.so  # → 272K
Load it in ROOT and dump entry 1: rec.common.ixn.dlp[0].part.dlp size (expect 12), with pdg/primary/E/start/end/p.
```

### 3. Set up the conda env (one-time)
```
Create the dune_root env that has ROOT + matplotlib:
  conda create -n dune_root -c conda-forge root python=3.11 matplotlib -y  # or use existing /Users/zdjurcic/miniconda3/envs/dune_root
  conda activate dune_root
  python -c "import ROOT, matplotlib; print(ROOT.gROOT.GetVersion(), matplotlib.__version__)"  # → 6.40.04
Every new shell: conda activate dune_root
```

### 4. Write the verified loose analyzer (12 files, 6 plots)
```
Write local_full/analyze_verified.py that runs with dune_root's python and uses ROOT.caf.StandardRecord via libduneanaobj.so:
- Loop 12 files in caf_MC/, for each spill: truth loop mc.nu.{E, vtx, mode, iscc, nprim} → Enu/vtx/mode/nprim truth
- Reco loop common.ixn.dlp: fiducial 5<|x|<59 etc. (L2461), track primary+|start-vtx|<5 len>3 (L2573), MINERvA loose (exiting |z|>58 proxy)
- Truth match via truth/truthOverlap argmax + |Δvtx|<5 → Enu_truth_fid + CC/NC split
- Fill hadron_mult (p/π/K) loose, track_length and angular via |end-start| and pz/|p| for μ/p/π are stubs for now
- Plots: Enu_truth, vertex X/Y/Z, mode, nprim → plots_verified/ + cutflow_verified.json (expect 868 fiducial → 349 track → 166 matched)
Run: /Users/zdjurcic/miniconda3/envs/dune_root/bin/python local_full/analyze_verified.py (15s, 6 PNGs)
```

### 5. Add tight MINERvA (full dict)
```
Extend to analyze_verified_minerva.py: build libduneanaobj_full.so with LinkDef_full.h (+SRTrack, SRMINERvA, SRTMS, SRSAND, SRGAr, SRNDLAr, SRTruthBranch).
Implement passes_minerva_tight(ixn, rec.nd.minerva) replicating L2625-2785: for each exiting reco track (primary, |start-vtx|<5, len>3, |z|>58),
loop nd.minerva.ixn tracks, lenMin>10, dot>prev && |ΔY|<15 && |ΔX|<15 && |Δatan|<0.06, pass maxDot>0.99 (mcOnly → mnvOffset 0,0).
Gate on it (requireMinerva=True). Expect loose 166 matched → tight 4 matched on 12 files, and now fill hadron_mult/length/angular fully.
```

### 6. Verify entry 1 (hadron mult check)
```
Run the tight analyzer on entry 1 and print:
  Entry 1 nDLP=1, ndlp=12: list pdg/primary/len/costh for each part.dlp
Check: 12 parts → 2μ (pdg 13, len 55/35 cm) + 1p (2212, 55 cm) + 1π (211, 4.9 cm) primary, hadron_mult reco = 2 (μ excluded) or 4 with secondaries — matches your 2p+2π+1K → 4 expectation before primary filter.
```

### 7. Make the full 13-plot suite
```
In analyze_verified_minerva.py, after the tight match, loop nu.prim for truth hadron_mult/length/angular (same len>3, pdg 13/211/2212/321, sp→ep distance, p.pz/|p|).
Fill both truth and reco: hadron_mult (0-10), length μ/p/π (0-100 cm, 50 bins), angular μ/p/π (-1 to 1, 50 bins).
Add plots: hadron_mult_truth_vs_reco.png, length_mu/proton/pion_truth_vs_reco.png, angular_mu/proton/pion_truth_vs_reco.png.
Total should be 13 PNGs in plots_verified/ (6 from step 4 + 7 new).
```

### 8. Set/change selection (how to tune cuts)
```
Show the 6 cuts at the top of analyze_verified_minerva.py and how to change them:
- Fiducial: if not (abs(vx)>5 and abs(vx)<59 ...) → 5→10 for tighter
- Track: minTrkLength=3.0, diff>5, abs(pdg) in (2212,13,211,321) → 3→5 cm, drop 321 for p/π only
- Hadron mult: if abs(pdg) in (2212,211,321) → remove 321
- MINERvA tight: passes_minerva_tight() → maxDotDS>0.99, |Δ|<15, |Δatan|<0.06 → 0.99→0.995
- MC vs Data: mnvOffsetX=0 Y=0 → -10,5 for data (L1989)
- Truth match: abs(vx-tvx)<5 → 5→2 cm
After edit, rerun: dune_root/bin/python analyze_verified_minerva.py and check plots_*.png + cutflow_*.json update.
```

### 9. Make it reproducible (HTML + prompts)
```
Write local_full/analysis_summary.html that embeds the 6 cut definitions, the why-ROOT-needed (AsArray), the 13 plots (via <img src="plots/...">), and the cutflow table (2368 spills → 868 fiducial → 349 track → 166 loose/4 tight).
Also write this PROMPTS.md (the 11 prompts you are reading) so a fresh model can redo steps 0-8 without help.
```

### 10. Save for GitHub (2x2_mc)
```
Copy local_full/ (minus caf_MC/*.root, keep plots/*.png + cutflow*.json + *.py + deps/src + libs) to nu_analysis/2x2_mc/ as 1:1.
Add sync script: nu_analysis/sync_local_full.sh to_nu / to_local (rsync, excludes ROOTs).
Commit: git add 2x2_mc/ && git commit -m "Initial commit: 2x2_mc verified 12-file" && git push
Local stays runnable: dune_root/bin/python nu_analysis/2x2_mc/analyze_verified_minerva.py
```

Run these 11 prompts in order — you get the same 12-file verified analysis (6+7=13 plots, 868 fiducial, 349 track, 166/4 matched) with ROOT + duneanaobj, tunable cuts, and a GitHub-ready nu_analysis/2x2_mc/ without ever needing the assistant.

