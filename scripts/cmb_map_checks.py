"""Localising the ACT / Planck difference in the deflection validation (iteration 10 follow-up).

Three tests on the combined convergence template T (kappa_combined_alm) and the two CMB convergence maps:
  (a) same sky: T x ACT and T x Planck on the identical footprint (tracer union mask, ACT mask > 0.99, Planck
      mask), with the same NaMaster bins and covariance, so that neither sky coverage nor sample variance of
      different regions can explain a difference in A_L;
  (b) ACT x Planck: the cross-spectrum of the two convergence maps on their overlap against the LCDM
      C_L^{kappa kappa} (Limber to recombination), an amplitude that tests the relative calibration of the two maps
      without any tracer, plus each map's auto-spectrum against theory + the released N_L (ACT);
  (c) the T x ACT and T x Planck bandpowers on their full overlaps, for the figure (report/figures/template_cmb_cross.pdf).
Writes results/cmb_map_checks.json.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
import numpy as np
import healpy as hp
import sys; from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # repository root: `lyalenser` imports without installation
from lyalenser.paths import DATA, RAW
from lyalenser.config import Config
from lyalenser.cosmo import chi as chi_of_z
from lyalenser.lensing import Z_CMB, kernel
from lyalenser.lowz import slice_spectra, ANNULUS
from lyalenser.nmt_spectra import Spectra, fit_amplitude
from lyalenser.cmb_maps import load_kappa, MASKED_ON_INPUT


def prediction(summary,cref,lmax,nside):
    ell=np.arange(lmax+1); pw=hp.pixwin(nside,lmax=lmax); pred=np.zeros(lmax+1)
    for s in summary['slices']:
        Lth,C=slice_spectra(s['zmin'],s['zmax'],cref,lmax); pred+=np.asarray(s['effective_weight'],float)[:lmax+1]*np.interp(ell,Lth,C[:,1,2])*pw
    return pred


def kk_theory(lmax):
    from lyalenser.lensing import limber
    ccmb=float(chi_of_z(Z_CMB)); L=np.unique(np.r_[2.,np.linspace(2,lmax,400)])
    c=limber(L,lambda x: kernel(x,ccmb),lambda x: kernel(x,ccmb),1.,ccmb,nchi=800,to_recombination=True)
    return np.interp(np.arange(lmax+1),L,c)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--lowz',type=Path,default=DATA/'lowz_v2'); ap.add_argument('--nside',type=int,default=512)
    ap.add_argument('--lmax',type=int,default=1000); ap.add_argument('--out',type=Path,default=ROOT/'results/cmb_map_checks.json')
    a=ap.parse_args(); cfg=Config(scale=1.,r_perp_min=3.,fit_rperp_min=3.); lmax=a.lmax; t0=time.perf_counter()
    S=Spectra(lmax,width=int(ANNULUS)); summary=json.loads((a.lowz/'summary.json').read_text()); cfg=cfg.copy(chi_ref=float(summary.get('chi_ref',cfg.chi_ref)))
    mT=hp.read_map(str(a.lowz/f'mask_combined_nside{a.nside}.fits')); alm=hp.read_alm(str(a.lowz/'kappa_combined_alm.fits'))
    T=hp.alm2map(hp.almxfl(alm,(np.arange(hp.Alm.getlmax(len(alm))+1)<=lmax).astype(float)),a.nside,verbose=False)
    pred=prediction(summary,cfg.chi_ref,lmax,a.nside); ckk=kk_theory(lmax)   # cfg.chi_ref follows the summary; ckk=kk_theory(lmax)
    maps={s:load_kappa(s,a.nside,lmax) for s in ('ACT','Planck')}
    use=lambda Tb: (S.ell_eff>=40)&(S.ell_eff<=500)&(Tb!=0)
    out={'nside':a.nside,'lmax':lmax,'L':S.ell_eff.tolist()}
    from lyalenser.template_prediction import TemplatePrediction
    TP=TemplatePrediction(a.lowz,cfg.chi_ref,lmax,a.nside)
    def cross(name,f1,f2,theory,keyT):
        if keyT=='T': theory=TP.cross(f1.get_mask()*f2.get_mask())      # class fractions inside this overlap
        cx=S.cross(f1,f2)[0]; th=S.theory(f1,f2,[theory])[0]; a11=S.cross(f1,f1)[0]; a22=S.cross(f2,f2)[0]
        cov=S.gaussian_covariance(f1,f2,f1,f2,S.spectrum_model(a11,0.),S.spectrum_model(cx),S.spectrum_model(cx),S.spectrum_model(a22,0.))
        u=use(th); A,sA,chi2,dof=fit_amplitude(cx,th,cov,u)
        print(f"  {name:34s} A = {A:.3f} +- {sA:.3f} (chi2 {chi2:.1f}/{dof}) ({time.perf_counter()-t0:.0f} s)",flush=True)
        return {'A':A,'sigma_A':sA,'chi2':chi2,'dof':dof,'cross':cx.tolist(),'err':np.sqrt(np.diag(cov)).tolist(),'theory':th.tolist(),'auto_1':a11.tolist(),'auto_2':a22.tolist(),'fsky_joint':float(np.mean(f1.get_mask()*f2.get_mask()))}
    # (c) full overlaps
    fT=S.field(mT,[T],key='T'); out['full_overlap']={}
    for s,(k,mk) in maps.items():
        fk=S.field(mk,[k],key=s,masked_on_input=MASKED_ON_INPUT[s]); out['full_overlap'][s]=cross(f'T x {s} (full overlap)',fT,fk,pred,'T')
    # (a) same sky
    mA=maps['ACT'][1]; joint=(mT>0.5)&(mA>0.99**2)&(maps['Planck'][1]>0.5); jm=joint.astype(float)
    print(f"joint footprint fsky {jm.mean():.4f} ({time.perf_counter()-t0:.0f} s)",flush=True)
    fTj=S.field(jm,[T],key='joint'); out['same_sky']={'fsky':float(jm.mean())}
    for s,(k,mk) in maps.items():
        fk=S.field(jm,[k],key='joint'); out['same_sky'][s]=cross(f'T x {s} (same sky)',fTj,fk,pred,'T')
    # (b) ACT x Planck and the autos
    fA=S.field(maps['ACT'][1],[maps['ACT'][0]],key='ACT',masked_on_input=True); fP=S.field(maps['Planck'][1],[maps['Planck'][0]],key='Planck')
    out['act_x_planck']=cross('ACT x Planck vs LCDM C_kk',fA,fP,ckk,'k')
    nl=np.loadtxt(str(RAW/'act/baseline/N_L_kk_act_dr6_lensing_v1_baseline.txt')); NL=np.interp(np.arange(lmax+1),nl[:,0],nl[:,1])
    aa=S.cross(fA,fA)[0]; th=S.theory(fA,fA,[ckk])[0]; thn=S.theory(fA,fA,[ckk+NL])[0]
    out['act_auto']={'auto':aa.tolist(),'theory_signal':th.tolist(),'theory_signal_plus_NL':thn.tolist(),'ratio_40_500':float(np.mean((aa/thn)[use(thn)]))}
    print(f"  ACT auto / (C_kk + N_L) over 40-500: {out['act_auto']['ratio_40_500']:.3f}",flush=True)
    pp=S.cross(fP,fP)[0]; thp=S.theory(fP,fP,[ckk])[0]; out['planck_auto']={'auto':pp.tolist(),'theory_signal':thp.tolist()}
    a.out.write_text(json.dumps(out,indent=1,default=float)+'\n'); print('saved',a.out)


if __name__=='__main__': main()
