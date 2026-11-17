# 2x2_data — placeholder

Same fiducial+track+MINERvA selections as `2x2_mc/`, for data/MC overlay.

Structure will mirror `2x2_mc/`:
```
2x2_data/
  analyze_verified.py / analyze_verified_minerva.py  (data-aware: mnvOffsetX=-10 Y=5)
  caf_MC/  (data CAFs, not MC)
  plots/   (data vs MC)
```

To populate: copy `2x2_mc/analyze_verified*.py`, set `mnvOffsetX=-10, mnvOffsetY=5` (nominal_mc.cc L1989) and point `caf_dir` to data CAFs.
