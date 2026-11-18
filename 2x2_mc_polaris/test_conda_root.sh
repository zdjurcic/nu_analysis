#!/bin/bash
set -e
echo "=== conda test ==="
conda env list | head -20
echo "=== try root_env ==="
if [ -d /home/zdjurcic/root_env ]; then
  conda activate /home/zdjurcic/root_env
  echo "Activated $CONDA_PREFIX"
  python -c "import ROOT; print('ROOT OK', ROOT.gROOT.GetVersion())"
  python -c "import matplotlib; print('mpl', matplotlib.__version__)"
else
  echo "No /home/zdjurcic/root_env"
  ls -d /home/zdjurcic/*env* 2>&1 | head -20
  ls /lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/ 2>&1 | head -20
fi
echo "=== test 12 files ==="
python << 'PYEOF'
import ROOT
ROOT.gROOT.SetBatch(True)
# Load staged lib
import pathlib as pl
for cand in [pl.Path("/home/zdjurcic/polaris_full/deps/duneanaobj/libduneanaobj_full.so")]:
    if cand.exists():
        ROOT.gSystem.Load(str(cand.resolve()))
        print(f"Loaded {cand}")
for cand in [pl.Path("/home/zdjurcic/polaris_full/deps/duneanaobj/src")]:
    if cand.exists():
        ROOT.gInterpreter.AddIncludePath(str(cand))
f=ROOT.TFile.Open("/home/zdjurcic/caf_MC/MiniRun6.5_1E19_RHC.caf.0000000.CAF.root")
t=f.Get("cafTree")
rec=ROOT.caf.StandardRecord()
t.SetBranchAddress("rec", rec)
t.GetEntry(1)
print(f"Read OK nDLP={rec.common.ixn.dlp.size()} nMinerva={rec.nd.minerva.nixn}")
PYEOF
