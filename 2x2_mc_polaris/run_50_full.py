#!/usr/bin/env python3
"""
Full 50-file verification — all 13 plots as in local_full/analyze_verified_minerva.py
fiducial+track+MINERvA_tight → Enu, vertex X/Y/Z, mode, nprim, hadron_mult, length μ/p/π, angular μ/p/π
Writes to /lus/eagle/.../plots_50_full/
"""
import pathlib, glob, json, math
import ROOT
ROOT.gROOT.SetBatch(True)
import pathlib as pl
for cand in [pl.Path("/home/zdjurcic/2x2/deps/duneanaobj/install/lib/libduneanaobj_StandardRecord.so"),
             pl.Path("/home/zdjurcic/2x2/deps/duneanaobj/install/lib/libduneanaobj_StandardRecord_dict.so"),
             pl.Path("/home/zdjurcic/polaris_full/deps/duneanaobj/libduneanaobj_full.so")]:
    if cand.exists():
        try:
            ROOT.gSystem.Load(str(cand.resolve()))
            print(f"Loaded {cand}")
            break
        except: pass
for cand in [pl.Path("/home/zdjurcic/2x2/deps/duneanaobj/src"), pl.Path("/home/zdjurcic/polaris_full/deps/duneanaobj/src")]:
    if cand.exists():
        ROOT.gInterpreter.AddIncludePath(str(cand))
        print(f"Include {cand}")
        break
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
EAGLE=pathlib.Path("/lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/caf_MC/0000000")
files=sorted(glob.glob(str(EAGLE/"*.root")))[:50]
print(f"Full 50-file analysis using {len(files)} files from {EAGLE}")
outdir=pathlib.Path("/lus/eagle/clone/g2/projects/Nu_Novel/alcf_training/zdjurcic/caf_MC/plots_50_full")
outdir.mkdir(parents=True, exist_ok=True)
vals={k:[] for k in ["Enu_truth_all","Enu_truth_fid","Enu_truth_fid_cc","Enu_truth_fid_nc",
                     "vtxX_truth_all","vtxY_truth_all","vtxZ_truth_all","vtxX_truth_fid","vtxY_truth_fid","vtxZ_truth_fid",
                     "vtxX_reco_all","vtxY_reco_all","vtxZ_reco_all","vtxX_reco_fid","vtxY_reco_fid","vtxZ_reco_fid",
                     "mode_truth_all","mode_truth_fid","nprim_truth_all","nprim_truth_fid","nprim_reco_fid",
                     "hadron_mult_truth","hadron_mult_reco",
                     "length_mu_truth","length_mu_reco","length_proton_truth","length_proton_reco","length_pion_truth","length_pion_reco",
                     "angular_mu_truth","angular_mu_reco","angular_proton_truth","angular_proton_reco","angular_pion_truth","angular_pion_reco"]}
cutflow={"total_spills":0,"total_nu":0,"total_ixn_reco":0,"fiducial_reco":0,"fiducial_truth":0,"cc_fid":0,"nc_fid":0,"track":0,"minerva_loose":0,"minerva_tight":0,"matched":0,"hadron_mult_reco_entries":0,"length_reco_mu":0,"length_reco_p":0,"length_reco_pi":0}
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

for fpath in files:
    f=ROOT.TFile.Open(fpath)
    if not f or f.IsZombie(): continue
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
            vals["vtxX_truth_all"].append(float(nu.vtx.x)); vals["vtxY_truth_all"].append(float(nu.vtx.y)); vals["vtxZ_truth_all"].append(float(nu.vtx.z))
            vals["mode_truth_all"].append(int(nu.mode)); vals["nprim_truth_all"].append(int(nu.nprim))
            cutflow["total_nu"]+=1
            vx,vy,vz=float(nu.vtx.x),float(nu.vtx.y),float(nu.vtx.z)
            if abs(vx)>5 and abs(vx)<59 and abs(vz)>5 and abs(vz)<59.5 and abs(vy)<57: cutflow["fiducial_truth"]+=1
        n_ixn=rec.common.ixn.dlp.size()
        cutflow["total_ixn_reco"]+=n_ixn
        for ixn_idx in range(n_ixn):
            ixn=rec.common.ixn.dlp[ixn_idx]
            vx,vy,vz=float(ixn.vtx.x),float(ixn.vtx.y),float(ixn.vtx.z)
            vals["vtxX_reco_all"].append(vx); vals["vtxY_reco_all"].append(vy); vals["vtxZ_reco_all"].append(vz)
            if not (abs(vx)>5 and abs(vx)<59 and abs(vz)>5 and abs(vz)<59.5 and abs(vy)<57): continue
            cutflow["fiducial_reco"]+=1
            vals["vtxX_reco_fid"].append(vx); vals["vtxY_reco_fid"].append(vy); vals["vtxZ_reco_fid"].append(vz)
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
            vals["nprim_reco_fid"].append(int(ixn.part.dlp.size()))
            has_exiting=any(abs(float(ixn.part.dlp[p_idx].start.z))>58 or abs(float(ixn.part.dlp[p_idx].end.z))>58 for p_idx in range(ixn.part.dlp.size()) if ixn.part.dlp[p_idx].primary)
            if has_exiting: cutflow["minerva_loose"]+=1
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
                        vals["vtxX_truth_fid"].append(tvx); vals["vtxY_truth_fid"].append(tvy); vals["vtxZ_truth_fid"].append(tvz)
                        vals["mode_truth_fid"].append(int(nu.mode)); vals["nprim_truth_fid"].append(int(nu.nprim))
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
                        vals["hadron_mult_truth"].append(true_hadron); vals["hadron_mult_reco"].append(hadron_mult); cutflow["hadron_mult_reco_entries"]+=1
                        for l,c in lens_mu: vals["length_mu_reco"].append(l); vals["angular_mu_reco"].append(c)
                        for l,c in lens_p: vals["length_proton_reco"].append(l); vals["angular_proton_reco"].append(c)
                        for l,c in lens_pi: vals["length_pion_reco"].append(l); vals["angular_pion_reco"].append(c)
                        cutflow["length_reco_mu"]+=len(lens_mu); cutflow["length_reco_p"]+=len(lens_p); cutflow["length_reco_pi"]+=len(lens_pi)
    f.Close()
print("Cutflow 50 full:", cutflow)
import json
with open(outdir/"cutflow_50_full.json","w") as jf: json.dump(cutflow, jf, indent=2)
bins_Enu=np.linspace(0,20,50); bins_vtx=np.linspace(-70,70,70)
fig, ax=plt.subplots(figsize=(8,6))
ax.hist(vals["Enu_truth_all"], bins=bins_Enu, histtype='step', color='gray', label=f'All truth {len(vals["Enu_truth_all"])}')
ax.hist(vals["Enu_truth_fid"], bins=bins_Enu, histtype='step', color='blue', label=f'Fid+matched tight {len(vals["Enu_truth_fid"])}')
ax.hist(vals["Enu_truth_fid_cc"], bins=bins_Enu, histtype='step', color='red', linestyle='--', label=f'CC {len(vals["Enu_truth_fid_cc"])}')
ax.hist(vals["Enu_truth_fid_nc"], bins=bins_Enu, histtype='step', color='green', linestyle='--', label=f'NC {len(vals["Enu_truth_fid_nc"])}')
ax.set_xlabel('True E_nu [GeV]'); ax.set_ylabel('Events'); ax.set_title('E_nu 50 files tight (full)')
ax.legend(); fig.tight_layout(); fig.savefig(outdir/"Enu_50_full.png", dpi=150); plt.close(fig)
for coord in ["X","Y","Z"]:
    fig, ax=plt.subplots(figsize=(6,4))
    ax.hist(vals[f"vtx{coord}_truth_all"], bins=bins_vtx, histtype='step', color='gray', label='All truth')
    ax.hist(vals[f"vtx{coord}_truth_fid"], bins=bins_vtx, histtype='step', color='blue', label='Fid truth matched')
    ax.hist(vals[f"vtx{coord}_reco_fid"], bins=bins_vtx, histtype='step', color='darkorange', label='Fid reco')
    ax.set_xlabel(f'vtx {coord} [cm]'); ax.set_title(f'Vertex {coord} 50 files tight'); ax.legend()
    fig.tight_layout(); fig.savefig(outdir/f"vertex_{coord}_50_full.png", dpi=150); plt.close(fig)
fig, ax=plt.subplots(figsize=(8,5))
modes=sorted(set(vals["mode_truth_all"]))
for m in modes:
    ax.bar(str(m), vals["mode_truth_all"].count(m), color='lightgray', edgecolor='black')
    ax.bar(str(m), vals["mode_truth_fid"].count(m), color='blue', alpha=0.7, edgecolor='black')
ax.set_xlabel('GENIE mode'); ax.set_ylabel('Events'); ax.set_title('GENIE mode 50 files tight')
fig.tight_layout(); fig.savefig(outdir/"mode_50_full.png", dpi=150); plt.close(fig)
fig, ax=plt.subplots(figsize=(8,5))
ax.hist(vals["nprim_truth_all"], bins=range(0,11), histtype='step', color='gray', label=f'All truth {len(vals["nprim_truth_all"])}', align='left')
ax.hist(vals["nprim_truth_fid"], bins=range(0,11), histtype='step', color='blue', label=f'Fid truth {len(vals["nprim_truth_fid"])}', align='left')
ax.hist(vals["nprim_reco_fid"], bins=range(0,20), histtype='step', color='darkorange', label=f'Fid reco {len(vals["nprim_reco_fid"])}', align='left')
ax.set_xlabel('nPrim'); ax.set_ylabel('Events'); ax.set_title('nPrim 50 files tight')
ax.legend(); fig.tight_layout(); fig.savefig(outdir/"nprim_50_full.png", dpi=150); plt.close(fig)
fig, ax=plt.subplots(figsize=(8,5))
ax.hist(vals["hadron_mult_truth"], bins=range(0,11), histtype='step', color='blue', label=f'Truth {len(vals["hadron_mult_truth"])}', align='left', linewidth=1.5)
ax.hist(vals["hadron_mult_reco"], bins=range(0,11), histtype='step', color='darkorange', label=f'Reco tight {len(vals["hadron_mult_reco"])}', align='left', linewidth=1.5)
ax.set_xlabel('Hadron mult (p/π/K)'); ax.set_ylabel('Events'); ax.set_title('Hadron mult 50 files tight (full)')
ax.legend(); fig.tight_layout(); fig.savefig(outdir/"hadron_mult_50_full.png", dpi=150); plt.close(fig)
for species, label in [("mu","Muon"),("proton","Proton"),("pion","Pion")]:
    fig, ax=plt.subplots(figsize=(8,5))
    t_vals=vals[f"length_{species}_truth"]; r_vals=vals[f"length_{species}_reco"]
    if not t_vals and not r_vals: continue
    bins=np.linspace(0,100,50)
    ax.hist(t_vals, bins=bins, histtype='step', color='blue', label=f'Truth {label} {len(t_vals)}')
    ax.hist(r_vals, bins=bins, histtype='step', color='darkorange', label=f'Reco {label} {len(r_vals)}')
    ax.set_xlabel('Track length [cm]'); ax.set_ylabel('Tracks'); ax.set_title(f'Track length {label} 50 files tight (full)')
    ax.legend(); fig.tight_layout(); fig.savefig(outdir/f"length_{species}_50_full.png", dpi=150); plt.close(fig)
for species, label in [("mu","Muon"),("proton","Proton"),("pion","Pion")]:
    fig, ax=plt.subplots(figsize=(8,5))
    t_vals=vals[f"angular_{species}_truth"]; r_vals=vals[f"angular_{species}_reco"]
    if not t_vals and not r_vals: continue
    bins=np.linspace(-1,1,50)
    ax.hist(t_vals, bins=bins, histtype='step', color='blue', label=f'Truth {label} {len(t_vals)}')
    ax.hist(r_vals, bins=bins, histtype='step', color='darkorange', label=f'Reco {label} {len(r_vals)}')
    ax.set_xlabel('cosθ = pz/|p|'); ax.set_ylabel('Tracks'); ax.set_title(f'Angular {label} 50 files tight (full)')
    ax.legend(); fig.tight_layout(); fig.savefig(outdir/f"angular_{species}_50_full.png", dpi=150); plt.close(fig)
print(f"All 13 plots saved to {outdir}: {sorted([p.name for p in outdir.glob('*.png')])}")
