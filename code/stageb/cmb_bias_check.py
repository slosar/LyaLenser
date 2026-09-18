"""Sanity check of the tracer biases against CMB lensing (ACT DR6 baseline, Planck PR4): for every (tracer, slice)
unit-bias map m (lowz_catalogues.py), <m kappa_CMB>_ell = b x C_ell^{(W_lya slice) x (W_CMB slice)}, the (l, c)
element of the slice's Limber covariance. NaMaster decoupled bandpowers with the tracer mask on the map side and
the CMB mask (ACT: squared, per the release README) on the kappa side; the theory is pushed through the same
bandpower windows; errors from the NaMaster Gaussian covariance with the measured auto-spectra. The cross-spectrum
has no shot noise. Writes report/stageb/cmb_bias_check.json (ACT) and cmb_bias_check_planck.json.
Usage: python cmb_bias_check.py [--lowz DIR] [--nside 512] [--lmax 1000] [--surveys ACT Planck]
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import healpy as hp
HERE=Path(__file__).resolve().parent; CODE=HERE.parent
for p in (CODE,CODE/'pipeline',HERE):
    if str(p) not in sys.path: sys.path.insert(0,str(p))
from paths import DATA
from config import Config
from lowz import slice_spectra, bias_band, ANNULUS
from nmt_spectra import Spectra, fit_amplitude
from cmb_maps import load_kappa


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--lowz',type=Path,default=DATA/'lowz_split'); ap.add_argument('--nside',type=int,default=512)
    ap.add_argument('--lmax',type=int,default=1000); ap.add_argument('--out',type=Path,default=CODE.parent/'report/stageb/cmb_bias_check.json')
    ap.add_argument('--surveys',nargs='*',default=['ACT','Planck'])
    a=ap.parse_args(); cfg=Config(scale=1.,r_perp_min=3.,fit_rperp_min=3.); cref=cfg.chi_ref; lmax=a.lmax
    summary=json.loads((a.lowz/'summary.json').read_text()); pw=hp.pixwin(a.nside,lmax=lmax)
    for survey in a.surveys:
        t0=time.perf_counter(); S=Spectra(lmax,width=int(ANNULUS)); kmap,mk=load_kappa(survey,a.nside,lmax); fk=S.field(mk,[kmap],key=f'{survey}_mask')
        out={'survey':survey,'nside':a.nside,'lmax':lmax,'estimator':'NaMaster decoupled bandpowers, Gaussian covariance','tracers':{}}
        print(f"[{survey}] {'tracer':14s} {'b_auto':>8s} {'b_cmb':>8s} {'+-':>6s} {'ratio':>6s}  fsky_joint  chi2/dof")
        for sl in summary['slices']:
            Lth,C=slice_spectra(sl['zmin'],sl['zmax'],cref,lmax); Slc=np.interp(np.arange(lmax+1),Lth,C[:,1,2])
            band=bias_band(sl['zmin'],sl['zmax'])
            for lab,t in sl['tracers'].items():
                m,ml=hp.read_map(str(a.lowz/f'unitbias_{lab}_nside{a.nside}.fits'),field=(0,1)); w2=float(np.mean(ml*mk))
                if w2<=0: continue
                fm=S.field(ml,[m],key=f'{lab}_mask')
                cx=S.cross(fm,fk)[0]; T=S.theory(fm,fk,[Slc*pw])[0]; gg=S.cross(fm,fm)[0]; kk=S.cross(fk,fk)[0]
                cov=S.gaussian_covariance(fm,fk,fm,fk,S.spectrum_model(gg,0.),S.spectrum_model(cx),S.spectrum_model(cx),S.spectrum_model(kk,0.))
                use=(S.ell_eff>=band[0])&(S.ell_eff<=band[1])&(T>0)
                b,sb,chi2,dof=fit_amplitude(cx,T,cov,use)
                ba=t['bias_fit']['b']; err=np.sqrt(np.diag(cov))
                out['tracers'][lab]={'b_auto':ba,'sigma_b_auto':t['bias_fit']['sigma_b'],'b_cmb_cross':b,'sigma_b_cmb_cross':sb,'ratio_cross_over_auto':b/ba,
                                     'chi2':chi2,'dof':dof,'fsky_joint':w2,'band':list(band),'L':S.ell_eff[use].tolist(),'cross':cx[use].tolist(),'cross_err':err[use].tolist(),
                                     'theory_unit_bias':T[use].tolist(),'b_per_annulus':(cx/np.maximum(T,1e-30))[use].tolist(),
                                     'L_all':S.ell_eff.tolist(),'cross_all':cx.tolist(),'cross_err_all':err.tolist(),'theory_unit_bias_all':T.tolist()}
                print(f"[{survey}] {lab:14s} {ba:8.3f} {b:8.3f} {sb:6.3f} {b/ba:6.3f}  {w2:.3f}       {chi2:.1f}/{dof}   per annulus: {np.round((cx/np.maximum(T,1e-30))[use][:6],2)}",flush=True)
        out['wall_s']=time.perf_counter()-t0
        path=a.out if survey=='ACT' else a.out.with_name(a.out.stem+'_'+survey.lower()+a.out.suffix)
        path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(out,indent=1,default=float)+'\n'); print('saved',path)


if __name__=='__main__': main()
