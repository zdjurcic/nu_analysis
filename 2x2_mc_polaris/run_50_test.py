import pathlib, glob
EAGLE=pathlib.Path("/lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/caf_MC/0000000")
files=sorted(glob.glob(str(EAGLE/"*.root")))[:50]
print(f"50-test found {len(files)} files")
import ROOT
ROOT.gROOT.SetBatch(True)
import pathlib as pl
for cand in [pl.Path("/home/zdjurcic/polaris_full/deps/duneanaobj/libduneanaobj_full.so"), pl.Path("/home/zdjurcic/2x2/deps/duneanaobj/install/lib/libduneanaobj_StandardRecord.so")]:
    if cand.exists():
        ROOT.gSystem.Load(str(cand.resolve()))
        print(f"Loaded {cand}")
        break
for cand in [pl.Path("/home/zdjurcic/polaris_full/deps/duneanaobj/src"), pl.Path("/home/zdjurcic/2x2/deps/duneanaobj/src")]:
    if cand.exists():
        ROOT.gInterpreter.AddIncludePath(str(cand))
        print(f"Include {cand}")
        break
f=ROOT.TFile.Open(files[0])
t=f.Get("cafTree")
rec=ROOT.caf.StandardRecord()
t.SetBranchAddress("rec", rec)
t.GetEntry(1)
print(f"Read OK nDLP={rec.common.ixn.dlp.size()} nMinerva={rec.nd.minerva.nixn}")
print("50-file test imports OK")
