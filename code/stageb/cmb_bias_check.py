"""Sanity check of the tracer biases against ACT DR6 CMB lensing: for every (tracer, slice) unit-bias map m
(lowz_catalogues.py), <m kappa_CMB>_ell = b x C_ell^{(W_lya slice) x (W_CMB slice)}, the (l, c) element of the
slice's Limber covariance. Pseudo-C_ell with the ACT mask squared on the kappa side (ACT README) and the tracer
mask on the map side, divided by <M_K^2 M_L>; fsky approximation of the mask coupling (the auto-spectrum
transfer was 0.8-0.9 in the band); the cross-spectrum has no shot noise. Writes report/stageb/cmb_bias_check.json.
Usage: python cmb_bias_check.py [--lowz DIR] [--nside 512] [--lmax 1000]
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import healpy as hp
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DATA, RAW
from config import Config
from lowz import slice_spectra, bias_band, Tracer, TRACERS
from lowz_catalogues import binned_cl


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--lowz',type=Path,default=DATA/'lowz_split'); ap.add_argument('--nside',type=int,default=512)
    ap.add_argument('--lmax',type=int,default=1000); ap.add_argument('--out',type=Path,default=CODE.parent/'report/stageb/cmb_bias_check.json')
    a=ap.parse_args(); cfg=Config(scale=1.,r_perp_min=3.,fit_rperp_min=3.); cref=cfg.chi_ref; lmax=a.lmax
    act=RAW/'act/baseline'
    kalm=hp.read_alm(str(act/'kappa_alm_data_act_dr6_lensing_v1_baseline.fits')); kalm=np.nan_to_num(kalm)   # the release alm hold NaN in unused modes
    kalm=hp.almxfl(kalm,np.r_[np.zeros(2),np.ones(lmax-1),np.zeros(hp.Alm.getlmax(len(kalm))-lmax)])
    kmap=hp.alm2map(kalm,a.nside,verbose=False)   # the alm run to ell = 3000; modes above lmax were zeroed above
    mact=hp.ud_grade(hp.read_map(str(act/'mask_act_dr6_lensing_v1_healpix_nside_4096_baseline.fits')),a.nside)
    mk=mact**2   # ACT README: M_K^2 for the quadratic-estimator map
    edges=np.arange(0,lmax+40,40.); pw=hp.pixwin(a.nside,lmax=lmax)
    summary=json.loads((a.lowz/'summary.json').read_text()); out={'nside':a.nside,'lmax':lmax,'tracers':{}}
    print(f"{'tracer':14s} {'b_auto':>8s} {'b_cmb':>8s} {'+-':>6s} {'ratio':>6s}  fsky_joint  chi2/dof")
    for sl in summary['slices']:
        Lth,C=slice_spectra(sl['zmin'],sl['zmax'],cref,lmax); Slc=np.interp(np.arange(lmax+1),Lth,C[:,1,2]); Scc=np.interp(np.arange(lmax+1),Lth,C[:,2,2])
        band=bias_band(sl['zmin'],sl['zmax'])
        for lab,t in sl['tracers'].items():
            m,ml=hp.read_map(str(a.lowz/f'unitbias_{lab}_nside{a.nside}.fits'),field=(0,1)); joint=ml*mk; w2=float(np.mean(ml*mk))
            if w2<=0: continue
            cross=hp.anafast(m*ml,kmap*mk,lmax=lmax,iter=0)/w2; ag=hp.anafast(m*ml,lmax=lmax,iter=0)/max(float(np.mean(ml**2)),1e-30); ak=hp.anafast(kmap*mk,lmax=lmax,iter=0)/max(float(np.mean(mk**2)),1e-30)
            L,cx,nm=binned_cl(cross,lmax,edges); _,gg,_=binned_cl(ag,lmax,edges); _,kk,_=binned_cl(ak,lmax,edges); _,T,_=binned_cl(Slc*pw,lmax,edges)
            use=(L>=band[0])&(L<=band[1])&(T>0)&(nm>0); var=(gg*kk+cx**2)/np.maximum(nm*w2,1); w=1/np.maximum(var,1e-40)
            b=float(np.sum((w*cx*T)[use])/np.sum((w*T*T)[use])); sb=float(np.sqrt(1/np.sum((w*T*T)[use]))); chi2=float(np.sum((w*(cx-b*T)**2)[use]))
            ba=t['bias_fit']['b']; out['tracers'][lab]={'b_auto':ba,'sigma_b_auto':t['bias_fit']['sigma_b'],'b_cmb_cross':b,'sigma_b_cmb_cross':sb,'ratio_cross_over_auto':b/ba,
                                                        'chi2':chi2,'dof':int(use.sum()-1),'fsky_joint':w2,'band':list(band),'L':L[use].tolist(),'cross':cx[use].tolist(),'theory_unit_bias':T[use].tolist(),
                                                        'b_per_annulus':(cx/np.maximum(T,1e-30))[use].tolist()}
            print(f"{lab:14s} {ba:8.3f} {b:8.3f} {sb:6.3f} {b/ba:6.3f}  {w2:.3f}       {chi2:.1f}/{int(use.sum()-1)}   per annulus: {np.round((cx/np.maximum(T,1e-30))[use][:6],2)}",flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(out,indent=1,default=float)+'\n'); print('saved',a.out)


if __name__=='__main__': main()
