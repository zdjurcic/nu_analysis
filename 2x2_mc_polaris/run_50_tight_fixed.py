import pathlib, glob, json, math, traceback
import ROOT
ROOT.gROOT.SetBatch(True)
import pathlib as pl
for cand in [pl.Path("/home/zdjurcic/2x2/deps/duneanaobj/install/lib/libduneanaobj_StandardRecord.so"), pl.Path("/home/zdjurcic/2x2/deps/duneanaobj/install/lib/libduneanaobj_StandardRecord_dict.so"), pl.Path("/home/zdjurcic/polaris_full/deps/duneanaobj/libduneanaobj_full.so")]:
    if cand.exists():
        try:
            ROOT.gSystem.Load(str(cand.resolve()))
            print(f"Loaded {cand}", flush=True)
        except Exception as e:
            print(f"Load failed {cand}: {e}", flush=True)
        break
for cand in [pl.Path("/home/zdjurcic/2x2/deps/duneanaobj/src"), pl.Path("/home/zdjurcic/polaris_full/deps/duneanaobj/src")]:
    if cand.exists():
        ROOT.gInterpreter.AddIncludePath(str(cand))
        print(f"Include {cand}", flush=True)
        break
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
EAGLE=pathlib.Path("/lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/caf_MC/0000000")
files=sorted(glob.glob(str(EAGLE/"*.root")))[:50]
print(f"50-file tight analysis using {len(files)} files from {EAGLE}", flush=True)
outdir=pathlib.Path("/lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/caf_MC/plots_50_tight")
outdir.mkdir(parents=True, exist_ok=True)
print(f"Outdir {outdir} exists {outdir.exists()}", flush=True)
vals={k:[] for k in ["Enu_truth_all","Enu_truth_fid","Enu_truth_fid_cc","Enu_truth_fid_nc","hadron_mult_truth","hadron_mult_reco","length_mu_truth","length_mu_reco","length_proton_truth","length_proton_reco","length_pion_truth","length_pion_reco","angular_mu_truth","angular_mu_reco","angular_proton_truth","angular_proton_reco","angular_pion_truth","angular_pion_reco"]}
cutflow={"total_spills":0,"total_nu":0,"total_ixn_reco":0,"fiducial_reco":0,"track":0,"minerva_tight":0,"matched":0,"cc_fid":0,"nc_fid":0}
minTrkLength=3.0
mnvOffsetX=0; mnvOffsetY=0
def passes_minerva_tight(ixn, rec_nd_minerva):
    maxDotDS=-999
    vx,vy,vz=float(ixn.vtx.x),float(ixn.vtx.y),float(ixn.vtx.z)
    for p_idx in range(ixn.part.dlp.size()):
        p=ixn.part.dlp[p_idx]
        if not p.primary: continue
        pdg=int(p.pdg)
        if abs(pdg) not in (2212,13,211,321): continue
        diff=math.sqrt((p.start.x-vx)**2+(p.start.y-vy)**2+(p.start.z-vz)**2)
        if diff>5: continue
        dX=p.end.x-p.start.x; dY=p.end.y-p.start.y; dZ=p.end.z-p.start.z
        length=math.sqrt(dX*dX+dY*dY+dZ*dZ)
        if math.isnan(length) or math.isinf(length) or length<=minTrkLength: continue
        dirX=dX/length; dirY=dY/length; dirZ=dZ/length
        if dirZ<0: dirZ=-dirZ; dirX=-dirX; dirY=-dirY
        if not (abs(p.start.z)>58 or abs(p.end.z)>58): continue
        dotDS_best=-999
        for i in range(rec_nd_minerva.nixn):
            mi=rec_nd_minerva.ixn[i]
            for j in range(mi.ntracks):
                trk=mi.tracks[j]
                dir_z=float(trk.dir.z); end_z=float(trk.end.z); start_z=float(trk.start.z)
                if start_z>0 and ((p.start.z>60) or (p.end.z>60)):
                    dXMnv=trk.end.x-trk.start.x; dYMnv=trk.end.y-trk.start.y; dZMnv=trk.end.z-trk.start.z
                    lengthMinerva=math.sqrt(dXMnv*dXMnv+dYMnv*dYMnv+dZMnv*dZMnv)
                    if lengthMinerva<10: continue
                    dirXMin=dXMnv/lengthMinerva; dirYMin=dYMnv/lengthMinerva; dirZMin=dZMnv/lengthMinerva
                    dot=dirXMin*dirX+dirYMin*dirY+dirZMin*dirZ
                    end_pos_z=max(float(p.start.z),float(p.end.z)); end_pos_x=float(p.end.x) if float(p.end.z)>float(p.start.z) else float(p.start.x); end_pos_y=float(p.end.y) if float(p.end.z)>float(p.start.z) else float(p.start.y)
                    if abs(dirZ)<1e-6: continue
                    extrapdZ=start_z-end_pos_z; extrapY=dirY/dirZ*extrapdZ+end_pos_y; extrapX=dirX/dirZ*extrapdZ+end_pos_x
                    angXdiff=abs(math.atan(dirXMin/dirZMin)-math.atan(dirX/dirZ)) if abs(dirZMin)>1e-6 else 999
                    angYdiff=abs(math.atan(dirYMin/dirZMin)-math.atan(dirY/dirZ)) if abs(dirZMin)>1e-6 else 999
                    if dot>dotDS_best and abs(extrapY-mnvOffsetY)<15 and angXdiff<0.06 and angYdiff<0.06 and abs(extrapX-mnvOffsetX)<15:
                        dotDS_best=dot
        if dotDS_best>maxDotDS: maxDotDS=dotDS_best
    return maxDotDS>0.99

try:
    for fpath in files:
        f=ROOT.TFile.Open(fpath)
        if not f or f.IsZombie():
            print(f"SKIP zombie {fpath}", flush=True)
            continue
        t=f.Get("cafTree")
        rec=ROOT.caf.StandardRecord()
        t.SetBranchAddress("rec", rec)
        nentries=t.GetEntries()
        cutflow["total_spills"]+=nentries
        for entry in range(nentries):
            t.GetEntry(entry)
            for ni in range(rec.mc.nu.size()):
                nu=rec.mc.nu[ni]
                vals["Enu_truth_all"].append(float(nu.E))
                cutflow["total_nu"]+=1
            n_ixn=rec.common.ixn.dlp.size()
            cutflow["total_ixn_reco"]+=n_ixn
            for ixn_idx in range(n_ixn):
                ixn=rec.common.ixn.dlp[ixn_idx]
                vx,vy,vz=float(ixn.vtx.x),float(ixn.vtx.y),float(ixn.vtx.z)
                if not (abs(vx)>5 and abs(vx)<59 and abs(vz)>5 and abs(vz)<59.5 and abs(vy)<57): continue
                cutflow["fiducial_reco"]+=1
                passes_track=False; hadron_mult=0; lens_mu=[]; lens_p=[]; lens_pi=[]
                for p_idx in range(ixn.part.dlp.size()):
                    p=ixn.part.dlp[p_idx]
                    if not p.primary: continue
                    pdg=int(p.pdg)
                    if abs(pdg) not in (2212,13,211,321): continue
                    diff=math.sqrt((p.start.x-vx)**2+(p.start.y-vy)**2+(p.start.z-vz)**2)
                    if diff>5: continue
                    dX=p.end.x-p.start.x; dY=p.end.y-p.start.y; dZ=p.end.z-p.start.z
                    length=math.sqrt(dX*dX+dY*dY+dZ*dZ)
                    if math.isnan(length) or math.isinf(length) or length<=minTrkLength: continue
                    passes_track=True
                    if abs(pdg) in (2212,211,321): hadron_mult+=1
                    pmag=math.sqrt(p.p.x**2+p.p.y**2+p.p.z**2); costh=(p.p.z/pmag) if pmag>0 else -2
                    if abs(pdg)==13: lens_mu.append((length,costh))
                    elif abs(pdg)==2212: lens_p.append((length,costh))
                    elif abs(pdg)==211: lens_pi.append((length,costh))
                if not passes_track: continue
                cutflow["track"]+=1
                passes_tight=passes_minerva_tight(ixn, rec.nd.minerva)
                if passes_tight: cutflow["minerva_tight"]+=1
                else: continue
                if ixn.truth.size()>0:
                    best_idx=-1; best_overlap=-1
                    for tii in range(ixn.truth.size()):
                        ov=float(ixn.truthOverlap[tii])
                        if ov>best_overlap: best_overlap=ov; best_idx=int(ixn.truth[tii])
                    if 0 <= best_idx < rec.mc.nu.size():
                        nu=rec.mc.nu[best_idx]
                        tvx,tvy,tvz=float(nu.vtx.x),float(nu.vtx.y),float(nu.vtx.z)
                        if abs(vx-tvx)<5 and abs(vy-tvy)<5 and abs(vz-tvz)<5:
                            vals["Enu_truth_fid"].append(float(nu.E))
                            if bool(nu.iscc): vals["Enu_truth_fid_cc"].append(float(nu.E)); cutflow["cc_fid"]+=1
                            else: vals["Enu_truth_fid_nc"].append(float(nu.E)); cutflow["nc_fid"]+=1
                            cutflow["matched"]+=1
                            true_hadron=0
                            for k in range(nu.prim.size()):
                                pr=nu.prim[k]
                                apdg=abs(int(pr.pdg))
                                if apdg not in (13,211,2212,321): continue
                                sp=pr.start_pos; ep=pr.end_pos
                                if math.isnan(sp.z): continue
                                plen=math.sqrt((ep.x-sp.x)**2+(ep.y-sp.y)**2+(ep.z-sp.z)**2)
                                if plen<=minTrkLength: continue
                                if apdg in (2212,211,321): true_hadron+=1
                                pmag_t=math.sqrt(pr.p.px**2+pr.p.py**2+pr.p.pz**2) if pr.p else 0
                                costh_t=(pr.p.pz/pmag_t) if pmag_t>0 else -2
                                if apdg==13: vals["length_mu_truth"].append(plen); vals["angular_mu_truth"].append(costh_t)
                                elif apdg==2212: vals["length_proton_truth"].append(plen); vals["angular_proton_truth"].append(costh_t)
                                elif apdg==211: vals["length_pion_truth"].append(plen); vals["angular_pion_truth"].append(costh_t)
                            vals["hadron_mult_truth"].append(true_hadron); vals["hadron_mult_reco"].append(hadron_mult)
                            for l,c in lens_mu: vals["length_mu_reco"].append(l); vals["angular_mu_reco"].append(c)
                            for l,c in lens_p: vals["length_proton_reco"].append(l); vals["angular_proton_reco"].append(c)
                            for l,c in lens_pi: vals["length_pion_reco"].append(l); vals["angular_pion_reco"].append(c)
        f.Close()
        print(f"Done {pathlib.Path(fpath).name} cutflow {cutflow}", flush=True)
    print("Cutflow 50 tight:", cutflow, flush=True)
    import json
    with open(outdir/"cutflow_50_tight.json","w") as jf: json.dump(cutflow, jf, indent=2)
    print("Wrote cutflow", flush=True)
except Exception as e:
    print(f"ERROR {e}", flush=True)
    traceback.print_exc()
    raise
# Plots
try:
    bins_Enu=np.linspace(0,20,50)
    fig, ax=plt.subplots(figsize=(8,6))
    ax.hist(vals["Enu_truth_all"], bins=bins_Enu, histtype='step', color='gray', label=f'All truth {len(vals["Enu_truth_all"])}')
    ax.hist(vals["Enu_truth_fid"], bins=bins_Enu, histtype='step', color='blue', label=f'Fid+matched tight {len(vals["Enu_truth_fid"])}')
    ax.set_xlabel('True E_nu [GeV]'); ax.set_ylabel('Events'); ax.set_title('E_nu 50 files tight')
    ax.legend(); fig.tight_layout(); fig.savefig(outdir/"Enu_50_tight.png", dpi=150); plt.close(fig)
    print(f"Plots saved to {outdir}", flush=True)
except Exception as e:
    print(f"Plot error {e}", flush=True)
    traceback.print_exc()
