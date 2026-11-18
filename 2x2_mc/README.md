# Local Full CAF Analysis (no Polaris) — Verified 2026-11-16

This is the fresh subdirectory for the complete 2x2 CAF analysis in python,
as discussed: find reco and truth information for each event passing
fiducial+track+MINERvA selection (as in nominal_mc.cc) and make physics
distributions.

## Structure (mirrors nominal_mc.cc event loop)
- `caf_MC/` — 12 files from portal.nersc.gov (0000000-0000011, 21 MB) for local verification
- `nominal_mc.cc` — reference C++: shows where branches are and how truth is matched to reco
- `analyze_full.py` — uproot port (fiducial+track proxy, no SRRecoParticle)
- `analyze_verified.py` — **verified** ROOT port (real `part.dlp` via libduneanaobj.so, loose MINERvA)
- `analyze_verified_minerva.py` — **tight MINERvA** (replicates L2625-2785 dot>0.99, trans<15, ang<0.06)
- `deps/duneanaobj/libduneanaobj.so` (272K) — minimal dict; `libduneanaobj_full.so` (729K) — full nd.minerva dict

### Cuts (nominal_mc.cc)
1. **Truth loop:** `mc.nu.{E, vtx.x/y/z, mode, iscc, pdg, targetPDG, nprim}` → no cuts → `true_nuE_all` etc.
2. **Fiducial:** `common.ixn.dlp.vtx` in `5< |x|<59, 5< |z|<59.5, |y|<57` (L2461)
3. **Track:** `part.dlp.{primary, pdg, start, end}` with `|start-vtx|<5`, `length>3` (L2573, L2598)
4. **MINERvA tight:** exiting `|z|>58` + `nd.minerva.ixn` dot>0.99, `|trans|<15, |ang|<0.06` (L2625-2785); loose proxy counts exiting only
5. **Truth match:** `common.ixn.dlp.truth`/`truthOverlap` argmax → `mc.nu[best]` with `|recoVtx-trueVtx|<5` and `targetPDG==Ar` (L2494, L2491)
6. **Physics:** `hadron_mult` (p/π/K primary), `track_length=|end-start|`, `angular cosθ=pz/|p|` for μ/p/π truth vs reco

## Why ROOT is needed (not uproot)
`common.ixn.dlp.part.dlp` is `vector<vector<SRRecoParticle>>` with `SRVector3D` stored memberwise (`AsArray`).
uproot 5.7 raises `NotImplementedError: memberwise serialization of AsArray in caf::SRVector3D` — the fix is a
ROOT dict (`libduneanaobj.so` built via `rootcling` from `StandardRecord.h`). The verified scripts use
`ROOT.caf.StandardRecord` directly, so `part.dlp` reads correctly (validated on entry 1: 12 parts, 2μ+1p+1π primary etc.).

## Verified outputs (12 files, 2368 spills)
- `plots/` — 30+ PNGs: `Enu_truth`, `vertex_X/Y/Z`, `mode`, `nprim`, `hadron_mult`, `length_mu/proton/pion`, `angular_mu/proton/pion`
  - `*_verified.png` — loose MINERvA (349 track, 166 matched, 135 CC/31 NC)
  - `*_verified_minerva.png` — tight MINERvA (4 matched, 4 CC, 0 NC) — matches nominal_mc.cc `requireMinerva=True`
  - `cutflow_verified.json` / `cutflow_verified_minerva_tight.json`
- `plots_verified/` — same, source for `plots/` copy

### Cutflow (12-file slice)
```
total_spills 2368, total_nu 9619, total_ixn_reco 3212
fiducial_reco 868, fiducial_truth 205
track (real) 349  [vs 579 proxy in analyze_full.py]
minerva loose 295/349, tight 4/349  [tight = nominal_mc.cc dot>0.99]
matched loose 166 (135 CC/31 NC), tight 4 (4 CC/0 NC)
hadron_mult reco entries 166 loose / 4 tight
length reco μ 196, p 87, π 35 (loose)
```

## How to run (verified, today)
```bash
# From 05_Agentic_Workflows/
# Loose (stats for plots): 
/Users/zdjurcic/miniconda3/envs/dune_root/bin/python caf_analysis_code/local_full/analyze_verified.py

# Tight MINERvA (exact nominal_mc.cc):
/Users/zdjurcic/miniconda3/envs/dune_root/bin/python caf_analysis_code/local_full/analyze_verified_minerva.py

# Outputs land in caf_analysis_code/local_full/plots_verified/ and are copied to plots/
# Needs deps/duneanaobj/libduneanaobj*.so already built (rootcling + g++ as above)
```

Run: `python analyze_full.py` (uproot proxy, ~10s, no MINERvA-tight) — kept for reference.
Verified runs are ~15s for 12 files via ROOT in `dune_root` (ROOT 6.40.04, conda).

## Changing Event Selection

Edit the 6 cuts at the top of `analyze_verified.py` / `analyze_verified_minerva.py` (and mirror to `polaris_full/run_*.py` via `sync_local_full.sh to_nu`):

| Selection | Line to edit | Example |
|-----------|--------------|---------|
| **Fiducial** `5< |x|<59, 5< |z|<59.5, |y|<57` | `if not (abs(vx)>5 and abs(vx)<59 and abs(vz)>5 and abs(vz)<59.5 and abs(vy)<57):` | `5→10` tighter, `59→55` looser |
| **Track** `primary && |start-vtx|<5 && len>3` | `minTrkLength=3.0` + `diff>5` + `abs(pdg) in (2212,13,211,321)` | `3.0→5.0` cm, `5→3` cm, drop `321` to exclude kaons |
| **Hadron mult** `p/π/K` | `if abs(pdg) in (2212,211,321): hadron_mult+=1` | remove `321` for p/π only |
| **MINERvA tight** `maxDot>0.99, |Δ|<15, |Δatan|<0.06` | `passes_minerva_tight()` → `maxDotDS>0.99` + `abs(extrapY)<15` + `angXdiff<0.06` | `0.99→0.995`, `15→10` |
| **MC vs Data offset** | `mnvOffsetX=0; mnvOffsetY=0` | `-10,5` for data (L1989) |
| **Truth match** `|Δvtx|<5` | `if abs(vx-tvx)<5 and ...` | `5→2` cm stricter |

After edit:
```bash
/Users/zdjurcic/miniconda3/envs/dune_root/bin/python caf_analysis_code/local_full/analyze_verified_minerva.py
open caf_analysis_code/local_full/plots_verified/hadron_mult_truth_vs_reco.png
# For Polaris: stage_to_alcf(.../run_50_full.py → /home/zdjurcic/polaris_full/...) then submit_job with 30-min debug
```

See `../polaris_full/POLARIS_GUIDE.html` §5 and `PROMPTS.md` #9 for prompts to ask a model to edit and resubmit.
